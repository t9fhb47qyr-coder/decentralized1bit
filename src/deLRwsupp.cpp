#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include <vector>
#include "decentralizedCQR.h"

//' Decentralized Lasso Regression with Support Evaluation
//' @description
//' This function performs decentralized Lasso regression using proximal gradient descent with consensus penalties and evaluates support recovery.
//' @param X Design matrix (N x p), where N is the total number of samples and p is the number of features.
//' @param y Response vector (N x 1).
//' @param adjacency_matrix Adjacency matrix representing the network structure (m x m), where m is the number of nodes.
//' @param B_init Initial coefficient matrix (p x m).
//' @param betaT True coefficient vector for normalization (p x 1).
//' @param T_inner Maximum number of inner iterations.
//' @param tau_penalty_factor Factor to scale the consensus penalty.
//' @param nlambda Number of lambda values to consider for tuning.
//' @param lambda_factor Factor to determine the minimum lambda value.
//' @param lambda_max Maximum lambda value.
//' @param quiet If TRUE, suppresses output during iterations.
//' @param epsilon Convergence threshold for stopping criteria.
//' @return A list containing the estimated coefficients, error history, F1 score history, and results for each iteration.
//' @export
// [[Rcpp::export]]
Rcpp::List deLRwsupp(arma::mat &X, arma::vec &y,
                     arma::mat &adjacency_matrix,
                     arma::mat &B_init, arma::vec &betaT,
                     int T_inner = 20,
                     double tau_penalty_factor = 1 / 6.0,
                     int nlambda = 100,
                     double lambda_factor = 1e-4,
                     double lambda_max = 1.0,
                     bool quiet = true,
                     double epsilon = 1e-6) {
  // Parameters
  int m = adjacency_matrix.n_cols;
  int p = X.n_cols;
  int N = X.n_rows;
  int n = N / m;

  arma::vec s_vec = arma::eig_sym(arma::trans(X.rows(0, n - 1)) * X.rows(0, n - 1));
  double tau_penalty_beta = s_vec.back() / n * tau_penalty_factor;

  // Results
  arma::mat B_inner(p, m, arma::fill::zeros);
  arma::mat normalized_B(p, m, arma::fill::zeros);
  arma::mat B_inner_old(p, m, arma::fill::zeros);
  arma::mat P_beta(p, m, arma::fill::zeros);
  std::vector<double> error_history;
  error_history.reserve(T_inner);
  std::vector<double> f1_history;
  f1_history.reserve(T_inner);
  std::vector<arma::mat> result_history;
  result_history.reserve(T_inner);

  // Cache omega and rho
  arma::vec rho(m), omega(m);
  for (int j = 0; j < m; ++j) {
    arma::uvec idx = calN_j_cpp(n, j);
    s_vec = arma::eig_sym(arma::trans(X.rows(idx)) * X.rows(idx) / n) * tau_penalty_factor;
    rho(j) = s_vec.back();
    omega(j) = 1.0 / (2.0 * tau_penalty_beta * arma::accu(adjacency_matrix.row(j) != 0) + rho(j));
  }

  // Cache covariances
  arma::cube covariances(p, p, m, arma::fill::zeros);
  arma::mat Xty(p, m, arma::fill::zeros);
  for (int j = 0; j < m; ++j) {
    arma::uvec idx = calN_j_cpp(n, j);
    covariances.slice(j) = arma::trans(X.rows(idx)) * X.rows(idx);
    Xty.col(j) = arma::trans(X.rows(idx)) * y(idx);
  }

  // Tune parameter
  arma::vec bic_array(nlambda - 1, arma::fill::zeros);
  arma::vec lambda_array = arma::exp(arma::linspace(std::log(1.0), std::log(lambda_factor), nlambda));
  lambda_array.shed_row(0);
  lambda_array *= lambda_max;
  lambda_max = std::min(arma::max(arma::abs(X.t() * y / N)), lambda_max);

  if (!quiet) {
    Rcpp::Rcout << "lambda_max = " << lambda_max << std::endl;
  }

  arma::vec shat_array(nlambda - 1, arma::fill::zeros);

  for (int ilambda = 0; ilambda < lambda_array.n_elem; ++ilambda) {
    Rcpp::checkUserInterrupt();

    double lambda = lambda_array(ilambda);
    B_inner = B_init;
    P_beta.zeros();

    for (int t = 0; t < T_inner; ++t) {
      B_inner_old = B_inner;

      for (int j = 0; j < m; ++j) {
        arma::uvec idx = calN_j_cpp(n, j);
        arma::uvec neighbors_j = arma::find(adjacency_matrix.row(j));

        // Update P
        P_beta.col(j) += tau_penalty_beta *
          arma::sum(arma::repmat(B_inner_old.col(j), 1, neighbors_j.n_elem) - B_inner_old.cols(neighbors_j), 1);

        // Compute proximal update direction
        arma::vec update_vec = omega(j) * (rho(j) * B_inner_old.col(j)
          - 1.0 / (m * n) * (covariances.slice(j) * B_inner_old.col(j) - Xty.col(j))
          - P_beta.col(j)
          + tau_penalty_beta *
            arma::sum(arma::repmat(B_inner_old.col(j), 1, neighbors_j.n_elem) + B_inner_old.cols(neighbors_j), 1));

        // Update B
        B_inner.col(j) = soft_thresholding_cpp(update_vec, lambda * omega(j));
      }

      double change = arma::norm(B_inner - B_inner_old, "fro");
      if (epsilon > 0.0 && change <= epsilon) {
        break;
      }
    }

    // Calculate BIC
    int shat = 0;
    for (int j = 0; j < m; ++j) {
      arma::uvec idx = calN_j_cpp(n, j);
      bic_array(ilambda) += std::pow(arma::norm(y(idx) - X.rows(idx) * B_inner.col(j), "fro"), 2);
      shat = std::max(shat, int(arma::accu(B_inner.col(j) != 0)));
    }

    bic_array(ilambda) = bic_array(ilambda) / N + std::log(N) * std::log(std::log(N)) / N * std::max(shat, 1);
    shat_array(ilambda) = shat;
  }

  // Select best lambda
  arma::uword ilambda = arma::index_min(bic_array);
  double lambda = lambda_array(ilambda);

  if (!quiet) {
    Rcpp::Rcout << "lambda = " << lambda << std::endl;
  }

  B_inner = B_init;
  P_beta.zeros();

  for (int t = 0; t < T_inner; ++t) {
    B_inner_old = B_inner;

    for (int j = 0; j < m; ++j) {
      arma::uvec idx = calN_j_cpp(n, j);
      arma::uvec neighbors_j = arma::find(adjacency_matrix.row(j));

      // Update P
      P_beta.col(j) += tau_penalty_beta *
        arma::sum(arma::repmat(B_inner_old.col(j), 1, neighbors_j.n_elem) - B_inner_old.cols(neighbors_j), 1);

      // Compute proximal update direction
      arma::vec update_vec = omega(j) * (rho(j) * B_inner_old.col(j)
        - 1.0 / (m * n) * (covariances.slice(j) * B_inner_old.col(j) - Xty.col(j))
        - P_beta.col(j)
        + tau_penalty_beta *
          arma::sum(arma::repmat(B_inner_old.col(j), 1, neighbors_j.n_elem) + B_inner_old.cols(neighbors_j), 1));

      // Update B
      B_inner.col(j) = soft_thresholding_cpp(update_vec, lambda * omega(j));
    }

    for (int j = 0; j < m; ++j) {
      double col_norm = arma::norm(B_inner.col(j));
      if (col_norm > 0.0) {
        normalized_B.col(j) = B_inner.col(j) / col_norm;
      } else {
        normalized_B.col(j).zeros();
      }
    }

    double iteration_error = std::sqrt(arma::accu(arma::square(normalized_B - arma::repmat(betaT, 1, m))) / m);
    error_history.push_back(iteration_error);
    result_history.push_back(B_inner);

    const double support_tol = 1e-8;
    int tp = 0, fp = 0, fn = 0;
    for (arma::uword feature = 0; feature < betaT.n_elem; ++feature) {
      bool truth_active = std::abs(betaT(feature)) > support_tol;
      bool pred_active = arma::any(arma::abs(B_inner.row(feature)) > support_tol);

      if (truth_active && pred_active) {
        ++tp;
      } else if (!truth_active && pred_active) {
        ++fp;
      } else if (truth_active && !pred_active) {
        ++fn;
      }
    }

    double precision = (tp + fp) > 0 ? static_cast<double>(tp) / static_cast<double>(tp + fp) : 0.0;
    double recall = (tp + fn) > 0 ? static_cast<double>(tp) / static_cast<double>(tp + fn) : 0.0;
    double f1 = (precision + recall) > 0.0 ? 2.0 * precision * recall / (precision + recall) : 0.0;
    f1_history.push_back(f1);

    if (!quiet) {
      Rcpp::Rcout << iteration_error << "\t";
    }

    double change = arma::norm(B_inner - B_inner_old, "fro");
    if (epsilon > 0.0 && change <= epsilon) {
      break;
    }
  }

  Rcpp::NumericVector errors_inner = Rcpp::wrap(error_history);
  Rcpp::NumericVector f1_inner = Rcpp::wrap(f1_history);
  Rcpp::List result(result_history.size());
  for (std::size_t iter = 0; iter < result_history.size(); ++iter) {
    result[iter] = result_history[iter];
  }

  return Rcpp::List::create(
    Rcpp::Named("B") = B_inner,
    Rcpp::Named("history") = Rcpp::List::create(
      Rcpp::Named("errors_inner") = errors_inner,
      Rcpp::Named("F1_inner") = f1_inner),
    Rcpp::Named("result") = result);
}
