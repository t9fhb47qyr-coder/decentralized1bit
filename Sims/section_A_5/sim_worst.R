library(decentralized1bit)
library(igraph)
library(glmnet)


#' Run the worst-case decentralized 1-bit simulation (Section A.5).
#'
#' @param n Samples per node.
#' @param m Number of nodes.
#' @param N Total sample size; defaults to `m * n`.
#' @param s Sparsity level of the true model.
#' @param p Ambient dimension.
#' @param p_flip Probability of label flip noise.
#' @param sigma2 Variance parameter for Gaussian design.
#' @param rho Correlation parameter for Gaussian design.
#' @param noise_type Distribution of the observation noise.
#' @param case Scenario indicator passed to `genData()`.
#' @param pc Edge probability for the Erdős–Rényi communication graph.
#' @param T_inner Number of inner iterations for the decentralized solvers.
#' @param nlambda Number of regularization parameters.
#' @param quiet Logical; suppress verbose output in solvers.
#' @param batch Optional batch identifier (kept for compatibility).
#' @param seed Optional seed for reproducibility.
#'
#' @return A data.frame with error/F1 metrics for each estimator.
sim_worst <- function(n = 100,
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
                      batch = 1,
                      seed = NULL) {
  if (!is.null(seed)) {
    seed_backup <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (!is.null(seed_backup)) {
        assign(".Random.seed", seed_backup, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  invisible(batch) # retained for backward compatibility

  if (!isTRUE(all.equal(N, m * n))) {
    stop("Total sample size N must equal m * n for this simulation setup.")
  }

  data <- decentralized1bit::genData(
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
  beta_true <- data$betaT
  beta_norm <- sqrt(sum(beta_true^2))
  beta_unit <- if (beta_norm > 0) beta_true / beta_norm else beta_true
  support_true <- which(abs(beta_unit) > 0)

  adjacency_matrix <- as.matrix(igraph::as_adjacency_matrix(data$graph))

  B_init <- vapply(
    seq_len(m),
    function(node) {
      idx <- decentralized1bit::calN_j(n, node)
      decentralized1bit::bic.glmnet(
        X = X[idx, ],
        y = y[idx],
        nlambda = nlambda
      )
    },
    FUN.VALUE = numeric(p)
  )

  error_local <- compute_max_error(B_init, beta_unit)
  f1_local <- compute_min_f1(B_init, support_true)

  B_lr <- decentralized1bit::deLR(
    X = X,
    y = y,
    adjacency_matrix = adjacency_matrix,
    B_init = B_init,
    betaT = beta_unit,
    tau_penalty_factor = 10 / (3 * m),
    T_inner = T_inner,
    lambda_max = 0.01,
    nlambda = nlambda,
    quiet = quiet
  )$B
  error_lr <- compute_max_error(B_lr, beta_unit)
  f1_lr <- compute_min_f1(B_lr, support_true)

  B_subgd <- decentralized1bit::deSubGD_LS(
    X = X,
    y = y,
    adjacency_matrix = adjacency_matrix,
    B_init = B_init,
    betaT = beta_unit,
    T = T_inner,
    nlambda = nlambda,
    lambda_max = 20,
    quiet = quiet
  )$B
  error_subgd <- compute_max_error(B_subgd, beta_unit)
  f1_subgd <- compute_min_f1(B_subgd, support_true)

  methods <- factor(c("Local", "D-subGD", "deLR"), levels = c("Local", "D-subGD", "deLR"))
  data.frame(
    errors = c(error_local, error_subgd, error_lr),
    f1s = c(f1_local, f1_subgd, f1_lr),
    methods = methods
  )
}


compute_max_error <- function(param_matrix, true_params) {
  coefs <- ensure_matrix(param_matrix)
  penalties <- apply(
    coefs,
    2,
    function(beta_est) {
      if (all(beta_est == 0)) {
        return(sqrt(sum(true_params^2)))
      }
      err <- decentralized1bit::computeError(beta_est, true_params, type = "F")
      if (is.nan(err) || !is.finite(err)) {
        sqrt(sum((beta_est - true_params)^2))
      } else {
        err
      }
    }
  )
  max(penalties)
}


compute_min_f1 <- function(param_matrix, true_support) {
  coefs <- ensure_matrix(param_matrix)
  scores <- apply(
    coefs,
    2,
    function(beta_est) {
      estimated_support <- which(abs(beta_est) > 0)
      result <- decentralized1bit::computeF1(estimated_support, true_support)
      if (!is.null(result$f1)) result$f1 else NA_real_
    }
  )
  if (all(is.na(scores))) {
    NA_real_
  } else {
    min(scores, na.rm = TRUE)
  }
}


ensure_matrix <- function(x) {
  if (is.null(dim(x))) {
    matrix(x, ncol = 1)
  } else {
    x
  }
}
