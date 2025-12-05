#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <algorithm>
#include <cmath>
#include <vector>
#include "decentralizedCQR.h"

//' Decentralized Subgradient Descent with Least Squares Loss
//' @description This function performs decentralized subgradient descent for least squares loss optimization.
//' @param X Design matrix (N x p), where N is the total number of samples and p is the number of features.
//' @param y Response vector (N x 1).
//' @param adjacency_matrix Adjacency matrix representing the network structure (m x m), where m is the number of nodes.
//' @param B_init Initial coefficient matrix (p x m).
//' @param betaT True coefficient vector for normalization (p x 1).
//' @param T Maximum number of iterations.
//' @param a Parameter for stepsize calculation.
//' @param b Parameter for stepsize calculation.
//' @param lambda Penalty parameter.
//' @param eps Convergence threshold for stopping criteria.
//' @param quiet If TRUE, suppresses output during iterations.
//' @return A list containing the estimated coefficients, error history, and inner error history.
//' @export
// [[Rcpp::export]]
Rcpp::List deSubGD_LS_single(arma::mat &X,
                          arma::vec &y,
                   arma::mat &adjacency_matrix,
                   arma::mat &B_init,
                   arma::vec &betaT,
                   int T = 20,
                   double a = 3.5, // for stepsize
                   double b = 0.51, // for stepsize
                   double lambda = 1e-2, // penalty parameter
                   double eps = 1e-4,
                   bool quiet = true) {

  const int m = adjacency_matrix.n_cols;
  const int p = X.n_cols;
  const int N = X.n_rows;

  if (m <= 0) {
    Rcpp::stop("Adjacency matrix must have positive dimension.");
  }
  if (N % m != 0) {
    Rcpp::stop("Number of observations must be divisible by the number of nodes.");
  }
  const int n = N / m;

  if (betaT.n_elem != static_cast<arma::uword>(p)) {
    Rcpp::stop("betaT must have the same length as the number of columns in X.");
  }

  arma::mat g_mat(p, m); // subgradient matrix
  arma::mat B(size(B_init)), B_old(size(B_init));
  arma::mat Phi(size(B_init));

  arma::vec errors(T, arma::fill::zeros);
  arma::vec errors_inner(T, arma::fill::zeros);
  std::vector<arma::mat> history_B;
  history_B.reserve(T);

  std::vector<arma::uvec> local_indices(m);
  std::vector<arma::uvec> neighbors(m);
  arma::uvec neighbor_counts(m, arma::fill::zeros);

  for (int j = 0; j < m; ++j) {
    local_indices[j] = calN_j_cpp(n, j);
    neighbors[j] = arma::find(adjacency_matrix.row(j));
    neighbor_counts(j) = neighbors[j].n_elem;
  }

  arma::mat C(m, m, arma::fill::zeros);
  for (int j = 0; j < m; ++j) {
    const arma::uvec &neighbors_j = neighbors[j];
    for (arma::uword idx : neighbors_j) {
      const arma::uword denom = std::max(neighbor_counts(j), neighbor_counts(idx));
      if (denom > 0) {
        C(idx, j) = 1.0 / static_cast<double>(denom);
      }
    }
  }

  for (int j = 0; j < m; ++j) {
    const arma::uvec &neighbors_j = neighbors[j];
    if (neighbors_j.is_empty()) {
      C(j, j) = 1.0;
    } else {
      double col_sum = 0.0;
      for (arma::uword idx : neighbors_j) {
        col_sum += C(idx, j);
      }
      C(j, j) = 1.0 - col_sum;
    }
  }

  const double invN = 1.0 / static_cast<double>(N);
  const double invM = 1.0 / static_cast<double>(m);
  const arma::mat beta_mat = arma::repmat(betaT, 1, m);

  B = B_init;
  int t = 0;
  arma::mat normalized_B(p, m, arma::fill::zeros);

  for (; t < T; ++t) {
    const double eta = a * std::pow(1.0 / (t + 1), b);
    const double penalty_scale = lambda * invM;

    for (int j = 0; j < m; ++j) {
      const arma::uvec &idx = local_indices[j];
      const arma::mat X_sub = X.rows(idx);
      const arma::vec y_sub = y.elem(idx);
      const arma::vec residual = y_sub - X_sub * B.col(j);
      const arma::vec indicator = arma::conv_to<arma::vec>::from(residual <= 0.0);

      g_mat.col(j) = invN * X_sub.t() * indicator + penalty_scale * arma::sign(B.col(j));
    }

    Phi = B - eta * g_mat;
    B_old = B;
    B = Phi * C;

    history_B.push_back(B);

    const double frob_norm = arma::norm(B - B_old, "fro");
    errors(t) = frob_norm / std::sqrt(static_cast<double>(m));

    normalized_B.zeros();
    for (int j = 0; j < m; ++j) {
      const double col_norm = arma::norm(B.col(j));
      if (col_norm > 0.0) {
        normalized_B.col(j) = B.col(j) / col_norm;
      }
    }
    errors_inner(t) = std::sqrt(arma::accu(arma::square(normalized_B - beta_mat)) * invM);

    if (t > 0 && errors(t - 1) > 0.0) {
      const double rel_change = std::abs(errors(t) - errors(t - 1)) / errors(t - 1);
      if (rel_change < eps) {
        ++t;
        break;
      }
    }
  }

  const arma::uword total_iter = history_B.size();
  errors = errors.head(total_iter);
  errors_inner = errors_inner.head(total_iter);

  Rcpp::List result(total_iter);
  for (arma::uword iter = 0; iter < total_iter; ++iter) {
    result[iter] = history_B[iter];
  }

  return Rcpp::List::create(
    Rcpp::Named("B") = B,
    Rcpp::Named("history") = Rcpp::List::create(
      Rcpp::Named("errors") = errors,
      Rcpp::Named("errors_inner") = errors_inner,
      Rcpp::Named("result") = result));

}


//' @export
// [[Rcpp::export]]
Rcpp::List deSubGD_LS(arma::mat &X,
                          arma::vec &y,
                          arma::mat &adjacency_matrix,
                          arma::mat &B_init,
                          arma::vec &betaT,
                          int T = 20,
                          double a = 3.5, // for stepsize
                          double b = 0.51, // for stepsize
                          int K = 19, // quantile level
                          int nlambda = 100L,
                          double lambda_factor = 1e-4,
                          double lambda_max = 1,
                          double eps = 1e-2,
                          bool quiet = true) {

  (void)K; // parameter kept for interface compatibility

  const int m = adjacency_matrix.n_cols;
  const int N = X.n_rows;

  if (m <= 0) {
    Rcpp::stop("Adjacency matrix must have positive dimension.");
  }
  if (N % m != 0) {
    Rcpp::stop("Number of observations must be divisible by the number of nodes.");
  }
  if (N <= 1) {
    Rcpp::stop("Need at least two observations to compute the information criterion.");
  }
  if (nlambda < 2) {
    Rcpp::stop("nlambda must be at least 2.");
  }
  if (lambda_factor <= 0.0 || lambda_factor >= 1.0) {
    Rcpp::stop("lambda_factor must lie in (0, 1).");
  }
  if (lambda_max <= 0.0) {
    Rcpp::stop("lambda_max must be positive.");
  }

  const int n = N / m;
  const double invN = 1.0 / static_cast<double>(N);

  const arma::vec gradient = X.t() * y * invN;
  const double data_lambda_max = arma::max(arma::abs(gradient));
  lambda_max = std::min(data_lambda_max, lambda_max);

  arma::vec lambda_scaling = arma::exp(arma::linspace<arma::vec>(0.0, std::log(lambda_factor), nlambda));
  lambda_scaling.shed_row(0);
  arma::vec lambda_array = lambda_max * lambda_scaling;
  arma::vec bic_array(lambda_array.n_elem, arma::fill::zeros);

  if (!quiet) {
    Rcpp::Rcout << "lambda_max = " << lambda_max << std::endl;
  }

  std::vector<arma::uvec> local_indices(m);
  for (int j = 0; j < m; ++j) {
    local_indices[j] = calN_j_cpp(n, j);
  }

  Rcpp::List L;
  arma::mat B;

  for (arma::uword ilambda = 0; ilambda < lambda_array.n_elem; ++ilambda) {
    Rcpp::checkUserInterrupt();
    const double lambda = lambda_array(ilambda);

    L = deSubGD_LS_single(X,
                          y,
                          adjacency_matrix,
                          B_init,
                          betaT,
                          T,
                          a,
                          b,
                          lambda,
                          eps,
                          quiet);
    B = Rcpp::as<arma::mat>(L["B"]);

    double bic_value = 0.0;
    int shat = 0;
    for (int j = 0; j < m; ++j) {
      const arma::uvec &idx = local_indices[j];
      const arma::vec residual = y.elem(idx) - X.rows(idx) * B.col(j);
      bic_value += arma::dot(residual, residual);

      const int active = static_cast<int>(arma::accu(B.col(j) != 0.0));
      shat = std::max(shat, active);
    }
    bic_array(ilambda) = bic_value * invN + std::log(N) * std::log(std::log(N)) * invN * std::max(shat, 1);
  }

  const arma::uword ilambda = arma::index_min(bic_array);
  const double lambda = lambda_array(ilambda);

  if (!quiet) {
    Rcpp::Rcout << "lambda = " << lambda << std::endl;
  }

  return deSubGD_LS_single(X,
                           y,
                           adjacency_matrix,
                           B_init,
                           betaT,
                           T,
                           a,
                           b,
                           lambda,
                           eps,
                           quiet);
}
