#' Run a single decentralized 1-bit simulation instance.
#'
#' @param n Local sample size per node.
#' @param m Number of nodes.
#' @param N Total sample size (`m * n` by default).
#' @param s Sparsity level of the true signal.
#' @param p Ambient dimension.
#' @param p_flip Probability of flipping the response sign.
#' @param sigma2 Noise variance used when generating features.
#' @param rho Autocorrelation parameter for the covariance matrix.
#' @param noise_type Distribution used to draw noise terms.
#' @param case Integer indicating the data-generating scenario.
#' @param pc Proportion of informative features for graph construction.
#' @param T_inner Number of inner iterations for optimization routines.
#' @param nlambda Number of candidate lambda values for `glmnet`.
#' @param quiet Passed to iterative estimators to silence progress output.
#' @param batch Monte Carlo repetition identifier.
#'
#' @return Data frame with columns `errors`, `f1s`, and `methods`.

library(decentralized1bit)
library(glmnet)
library(igraph)

if (!exists("SIM_METHOD_LEVELS", inherits = TRUE)) {
  SIM_METHOD_LEVELS <- c("Pooled", "Local", "Avg", "D-subGD", "deLR")
}

if (!exists("summarise_slurm_results", inherits = TRUE)) {
  source("Sims/sim_utils.R")
}

sim_single_simulation <- function(n = 100,
                                  m = 20,
                                  N = m * n,
                                  s = 10,
                                  p = 100,
                                  p_flip = 0.05,
                                  sigma2 = 1,
                                  rho = 0.1,
                                  noise_type = "Normal",
                                  case = 0,
                                  pc = 0.3,
                                  T_inner = 500,
                                  nlambda = 100,
                                  quiet = TRUE,
                                  batch = 1) {
  compute_mean_f1 <- function(coef_matrix, support) {
    column_f1 <- apply(
      coef_matrix,
      2,
      function(coefs) computeF1(which(abs(coefs) > 0), support)[[1]]
    )
    mean(unlist(column_f1))
  }

  # Ensure reproducibility across Slurm array runs.
  RNGkind("L'Ecuyer-CMRG")
  set.seed(42)
  random_seed <- fixRNGStream(.Random.seed, batch)
  .GlobalEnv$.Random.seed <- random_seed

  # Simulate features, responses, and the underlying graph structure.
  data <- genData(
    N = N,
    m = m,
    p = p,
    s = s,
    pc = pc,
    p_flip = p_flip,
    rho = rho,
    sigma2 = sigma2,
    type = noise_type,
    case = case
  )

  X <- data$X
  y <- data$y
  beta_true <- data$betaT / norm(data$betaT, "F")
  support_true <- which(abs(beta_true) > 0)
  adjacency_matrix <- as.matrix(as_adjacency_matrix(data$graph))

  # Pooled estimator using the entire dataset.
  beta_pooled <- bic.glmnet(y = y, X = X, nlambda = nlambda)
  error_pooled <- computeError(beta_pooled, beta_true)
  f1_pooled <- computeF1(which(abs(beta_pooled) > 0), support_true)[[1]]

  # Local estimators computed independently on each node.
  B_init <- matrix(0, nrow = p, ncol = m)
  for (node in seq_len(m)) {
    local_index <- calN_j(n, node)
    B_init[, node] <- bic.glmnet(X = X[local_index, ], y = y[local_index], nlambda = nlambda)
  }
  error_local <- computeError(B_init, beta_true)
  f1_local <- compute_mean_f1(B_init, support_true)

  # Decentralized averaging estimator.
  B_average <- rowMeans(B_init)
  error_average <- computeError(B_average, beta_true)
  f1_average <- computeF1(which(abs(B_average) > 0), support_true)[[1]]

  # DeLR estimator.
  B_lr <- deLR(
    X = X,
    y = y,
    adjacency_matrix = adjacency_matrix,
    B_init = B_init,
    betaT = beta_true,
    tau_penalty_factor = 10 / (3 * m),
    T_inner = T_inner,
    lambda_max = 0.01,
    nlambda = nlambda,
    quiet = quiet
  )$B
  error_lr <- computeError(B_lr, beta_true)
  f1_lr <- compute_mean_f1(B_lr, support_true)

  # Decentralized subgradient estimator.
  B_subgd <- deSubGD_LS(
    X,
    y,
    adjacency_matrix,
    B_init = B_init,
    beta_true,
    T = T_inner,
    nlambda = nlambda,
    lambda_max = 20
  )$B
  error_subgd <- computeError(B_subgd, beta_true)
  f1_subgd <- compute_mean_f1(B_subgd, support_true)

  methods <- factor(SIM_METHOD_LEVELS, levels = SIM_METHOD_LEVELS)

  data.frame(
    errors = c(error_pooled, error_local, error_average, error_subgd, error_lr),
    f1s = c(f1_pooled, f1_local, f1_average, f1_subgd, f1_lr),
    methods = methods
  )
}
