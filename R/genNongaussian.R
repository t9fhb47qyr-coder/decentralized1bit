#' @title Generate Non-Gaussian Data for Distributed Binary Classification
#' @description This function generates non-Gaussian data for distributed binary classification tasks.
#' It allows for different types of error distributions and various cases of data heterogeneity across multiple machines.
#' @param N Total sample size.
#' @param m Number of machines.
#' @param p Number of features.
#' @param s Number of non-zero coefficients in the true signal.
#' @param pc Probability of edge creation in the communication graph.
#' @param p_flip Proportion of labels to be flipped (noise level).
#' @param rho Correlation coefficient for the feature covariance matrix.
#' @param sigma2 Variance parameter for the feature covariance matrix.
#' @param type Type of error distribution. Options are "Cauchy", "Normal", and "T2".
#' @param sigma2_set Set of variance parameters for heterogeneous data generation.
#' @param rho_set Set of correlation coefficients for heterogeneous data generation.
#' @param p_flip_set Set of label flip proportions for heterogeneous data generation.
#' @param case Case of data heterogeneity:
#'             0 - Homogeneous data across machines,
#'             1 - Heterogeneous data with varying covariance parameters,
#'             2 - Heterogeneous data with varying error distributions,
#'             3 - Homogeneous data with varying label noise levels across machines.
#' @return A list containing:
#'         \item{X}{Feature matrix of size N x p.}
#'         \item{y}{Label vector of size N x 1.}
#'         \item{betaT}{True signal vector of size p x 1.}
#'         \item{graph}{Communication graph among machines.}
#' @details Each machine receives an equal share of observations, so the total
#' sample size \code{N} must be divisible by the number of machines \code{m}.
#' The sparsity level \code{s} is also required to be no larger than the number
#' of features \code{p}.
#' @importFrom igraph sample_gnp is_connected
#' @importFrom MASS mvrnorm
#' @importFrom copula claytonCopula rCopula
#' @export
genNongaussianData <- function(N = 2000, m = 10, p = 100, s = 5,
                    pc = 0.3, p_flip = 0.2, rho = 0.1, sigma2 = 1,
                    type = c("Cauchy", "Normal", "T2"),
                    sigma2_set = c(1, 3),
                    rho_set = c(0.1, 0.3),
                    p_flip_set = c(0.05, 0.1),
                    case = 1) {
  type <- match.arg(type)

  if (any(c(N, m, p, s) <= 0L))
    stop("N, m, p, and s must be positive integers.")
  if (N %% m != 0)
    stop("N must be divisible by m so that each machine receives the same number of samples.")
  if (s > p)
    stop("The sparsity level s cannot exceed the number of features p.")

  n <- N %/% m

  sample_beta <- function() {
    beta <- numeric(p)
    support <- sample.int(p, s)
    beta[support] <- sample(c(-1, 1), s, replace = TRUE)
    beta
  }

  generate_noise <- function(n_obs, noise_type) {
    switch(noise_type,
      Cauchy = stats::rcauchy(n_obs),
      Normal = stats::rnorm(n_obs),
      T2 = stats::rt(n_obs, df = 2),
      stop("Unsupported noise type: ", noise_type)
    )
  }

  make_covariance <- function(rho_val, sigma_sq) {
    stats::toeplitz(rho_val ^ (0:(p - 1))) * sigma_sq
  }

  add_noise_and_sign <- function(linear_predictor, noise_vec) {
    response <- sign(drop(linear_predictor) + noise_vec)
    ifelse(response == 0, 1, response)
  }

  flip_labels <- function(labels, flip_fraction) {
    if (flip_fraction <= 0)
      return(labels)
    n_flip <- floor(length(labels) * flip_fraction)
    if (n_flip > 0) {
      flip_idx <- sample(seq_along(labels), n_flip)
      labels[flip_idx] <- -labels[flip_idx]
    }
    labels
  }

  betaT <- sample_beta()

  if (case == 0) {
    clayton <- copula::claytonCopula(param = 2, dim = p)
    X <- copula::rCopula(N, clayton)
    X <- scale(X, center = TRUE, scale = FALSE)
    noise <- generate_noise(N, type)
    y <- add_noise_and_sign(X %*% betaT, noise)
    y <- flip_labels(y, p_flip)
  } else if (case == 1) {
    betaT <- sample_beta()
    X <- matrix(NA_real_, nrow = N, ncol = p)
    noise <- generate_noise(N, type)
    for (j in seq_len(m)) {
      idx <- ((j - 1) * n + 1):(j * n)
      sigma2_j <- sample(sigma2_set, 1)
      rho_j <- sample(rho_set, 1)
      Sigma_j <- make_covariance(rho_j, sigma2_j)
      X[idx, ] <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma_j)
    }
    y <- add_noise_and_sign(X %*% betaT, noise)
    y <- flip_labels(y, p_flip)
  } else if (case == 2) {
    betaT <- sample_beta()
    Sigma <- make_covariance(rho, sigma2)
    X <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = Sigma)
    noise_types <- sample(c("Cauchy", "Normal", "T2"), m, replace = TRUE)
    noise <- unlist(lapply(noise_types, function(noise_type) generate_noise(n, noise_type)))
    y <- add_noise_and_sign(X %*% betaT, noise)
    y <- flip_labels(y, p_flip)
  } else if (case == 3) {
    betaT <- sample_beta()
    Sigma <- make_covariance(rho, sigma2)
    X <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = Sigma)
    noise <- generate_noise(N, type)
    y <- add_noise_and_sign(X %*% betaT, noise)
    for (j in seq_len(m)) {
      idx <- ((j - 1) * n + 1):(j * n)
      flip_fraction <- sample(p_flip_set, 1)
      if (flip_fraction > 0) {
        n_flip <- floor(length(idx) * flip_fraction)
        if (n_flip > 0) {
          flip_idx <- sample(idx, n_flip)
          y[flip_idx] <- -y[flip_idx]
        }
      }
    }
  } else {
    stop("Unsupported heterogeneity case: ", case)
  }

  graph <- igraph::sample_gnp(m, pc)
  while (!igraph::is_connected(graph)) {
    graph <- igraph::sample_gnp(m, pc)
  }

  list(
    X = X,
    y = y,
    betaT = betaT,
    graph = graph
  )
}
