#include <iostream>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

//' @title Element-wise maximum between a vector and a scalar
//' @description
//' This function takes a vector and a scalar as input and returns a new vector
//' where each element is the maximum between the corresponding element in the input
//' vector and the scalar.
//' @param x A numeric vector
//' @param bound A numeric scalar
//' @return A numeric vector with element-wise maximum values
//' @examples
//' x <- c(1, 5, 3, 7, 2)
//' bound <- 4
//' pmax_arma(x, bound)
//' @export
// [[Rcpp::export]]
arma::vec pmax_arma(arma::vec x, double bound) {
  int length_x = x.n_elem;
  for (int i = 0; i < length_x; i++) {
    if (x(i) < bound) x(i) = bound;
  }
  return x;
}

//' @title Epanechnikov Kernel Density Estimation
//' @description
//' This function performs kernel density estimation using the Epanechnikov kernel.
//' @param E A matrix of data points where each row represents a data point
//' @param h A vector of bandwidths for each dimension
//' @return A vector of estimated density values for each dimension
//' @examples
//' set.seed(1)
//' E <- matrix(rnorm(100), nrow=50, ncol=2)
//' h <- c(0.5, 0.5)
//' kernel(E, h)
//' @export
// [[Rcpp::export]]
arma::mat kernel(arma::mat E, arma::vec h)
{
  arma::mat f(arma::size(E));
  E = E / arma::repmat(h.t(), E.n_rows, 1);
  // E.each_row() /= h.t();
  f = (-315 * arma::pow(E, 6) + 735 * arma::pow(E, 4) -
       525 * arma::pow(E, 2) + 105) /
      64 % (arma::abs(E) < 1);
  return arma::trans(arma::mean(f, 0)) / h;
}

//' @title Mode of a vector
//' @description
//' This function calculates the mode (most frequently occurring value) of a numeric vector.
//' @param v A numeric vector
//' @return The mode of the vector
//' @examples
//' v <- c(1, 2, 2, 3, 4, 4, 4, 5)
//' mod(v)
//' @export
// [[Rcpp::export]]
double mod(arma::vec v) {
  std::map<int, size_t> frequencyCount;
  using pair_type = decltype(frequencyCount)::value_type;

  for (auto i : v)
    frequencyCount[i]++;

  auto pr = std::max_element
    (
        std::begin(frequencyCount), std::end(frequencyCount),
        [] (const pair_type & p1, const pair_type & p2) {
          return p1.second < p2.second;
        }
    );
  return pr->first;
}

//' @title Calculate local indices for node j
//' @description
//' Calculate the local indices for node j in a distributed setting
//' @param n Number of samples per node
//' @param j Node index (0-based)
//' @return A vector of local indices for node j
//' @examples
//' n <- 5
//' j <- 2
//' calN_j_cpp(n, j)
//' @export
// [[Rcpp::export]]
arma::uvec calN_j_cpp(int n, int j) {
  return arma::regspace<arma::uvec>(n * (j), 1, n * (j + 1) - 1);
}


//' Soft-thresholding operator
//' @description Apply the soft-thresholding operator to a given vector.
//' @param x A vector
//' @param t Thresholding parameter
//' @return The result after applying the soft-thresholding operator
//' @examples
//' x <- c(-3, -1, 0, 1, 3)
//' t <- 1.5
//' soft_thresholding_cpp(x, t)
//' @export
// [[Rcpp::export]]
arma::vec soft_thresholding_cpp(arma::vec x, double t) {
  return pmax_arma(x - t, 0) - pmax_arma(-x - t, 0);
}