library(decentralized1bit)
library(glmnet)
library(igraph)
library(NonNorMvtDist)

source("Sims/sim_utils.R")
source("R/genNongaussian.R")

compute_mean_f1 <- function(coef_matrix, support_true) {
  column_f1 <- apply(
    coef_matrix,
    2,
    function(coefs) decentralized1bit::computeF1(which(abs(coefs) > 0), support_true)[[1]]
  )
  mean(unlist(column_f1))
}

run_estimators <- function(X,
                           y,
                           beta_true,
                           graph,
                           config) {
  beta_unit <- beta_true / base::norm(beta_true, "F")
  support_true <- which(abs(beta_unit) > 0)
  adjacency_matrix <- as.matrix(igraph::as_adjacency_matrix(graph))

  B_init <- matrix(0, nrow = ncol(X), ncol = config$m)
  for (node in seq_len(config$m)) {
    idx <- decentralized1bit::calN_j(config$n, node)
    B_init[, node] <- decentralized1bit::bic.glmnet(
      X = X[idx, ],
      y = y[idx],
      nlambda = config$nlambda
    )
  }

  error_local <- decentralized1bit::computeError(B_init, beta_unit)
  f1_local <- compute_mean_f1(B_init, support_true)

  B_average <- rowMeans(B_init)
  error_average <- decentralized1bit::computeError(B_average, beta_unit)
  f1_average <- decentralized1bit::computeF1(which(abs(B_average) > 0), support_true)[[1]]

  beta_pooled <- decentralized1bit::bic.glmnet(y = y, X = X, nlambda = config$nlambda)
  error_pooled <- decentralized1bit::computeError(beta_pooled, beta_unit)
  f1_pooled <- decentralized1bit::computeF1(which(abs(beta_pooled) > 0), support_true)[[1]]

  B_lr <- decentralized1bit::deLR(
    X = X,
    y = y,
    adjacency_matrix = adjacency_matrix,
    B_init = B_init,
    betaT = beta_unit,
    tau_penalty_factor = 10 / (3 * config$m),
    T_inner = config$T_inner,
    lambda_max = 0.01,
    nlambda = config$nlambda,
    quiet = config$quiet
  )$B
  error_lr <- decentralized1bit::computeError(B_lr, beta_unit)
  f1_lr <- compute_mean_f1(B_lr, support_true)

  B_subgd <- decentralized1bit::deSubGD_LS(
    X = X,
    y = y,
    adjacency_matrix = adjacency_matrix,
    B_init = B_init,
    betaT = beta_unit,
    T = config$T_inner,
    nlambda = config$nlambda,
    lambda_max = 20
  )$B
  error_subgd <- decentralized1bit::computeError(B_subgd, beta_unit)
  f1_subgd <- compute_mean_f1(B_subgd, support_true)

  data.frame(
    error_pooled = error_pooled,
    error_local = error_local,
    error_avg = error_average,
    error_subGD = error_subgd,
    error_deLR = error_lr,
    f1_pooled = f1_pooled,
    f1_local = f1_local,
    f1_avg = f1_average,
    f1_subGD = f1_subgd,
    f1_deLR = f1_lr
  )
}

generate_standard_data <- function(config) {
  decentralized1bit::genData(
    N = config$N,
    m = config$m,
    p = config$p,
    s = config$s,
    pc = config$pc,
    p_flip = config$p_flip,
    rho = config$rho,
    sigma2 = config$sigma2,
    type = config$noise_type,
    case = config$case
  )
}

generate_nongaussian_data <- function(config) {
  genNongaussianData(
    N = config$N,
    m = config$m,
    p = config$p,
    s = config$s,
    pc = config$pc,
    p_flip = config$p_flip,
    rho = config$rho,
    sigma2 = config$sigma2,
    type = config$noise_type,
    case = config$case
  )
}

single_sim <- function(n = 100,
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
  config <- list(
    n = n,
    m = m,
    N = N,
    s = s,
    p = p,
    p_flip = p_flip,
    sigma2 = sigma2,
    rho = rho,
    noise_type = noise_type,
    case = case,
    pc = pc,
    T_inner = T_inner,
    nlambda = nlambda,
    quiet = quiet
  )

  data <- generate_standard_data(config)
  run_estimators(data$X, data$y, data$betaT, data$graph, config)
}

single_sim_non <- function(n = 100,
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
  config <- list(
    n = n,
    m = m,
    N = N,
    s = s,
    p = p,
    p_flip = p_flip,
    sigma2 = sigma2,
    rho = rho,
    noise_type = noise_type,
    case = case,
    pc = pc,
    T_inner = T_inner,
    nlambda = nlambda,
    quiet = quiet
  )

  data <- generate_nongaussian_data(config)
  run_estimators(data$X, data$y, data$betaT, data$graph, config)
}

