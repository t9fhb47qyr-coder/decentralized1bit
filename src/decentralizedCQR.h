#ifndef CQR_H
#define CQR_H
#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

arma::mat kernel(arma::mat E, arma::vec h);

double mod(arma::vec v);

arma::uvec calN_j_cpp(int n, int j);

arma::vec pmax_arma(arma::vec x, double bound);

arma::vec soft_thresholding_cpp(arma::vec x, double t);

#endif
