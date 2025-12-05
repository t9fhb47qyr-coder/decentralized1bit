library(decentralized1bit)
library(glmnet)
library(igraph)

if (!exists("SIM_METHOD_LEVELS", inherits = TRUE)) {
  SIM_METHOD_LEVELS <- c("Pooled", "Local", "Avg", "D-subGD", "deLR")
}

if (!exists("summarise_support_recovery_results", inherits = TRUE)) {
  source("Sims/sim_utils.R")
}

#' Run a single simulation focused on support recovery metrics.
#'
#' @inheritParams sim_single_simulation
#' @return Data frame with columns `exact_support`, `recall`, `precision`,
#'   and `methods`.
sim_single_support_recovery <- function(n = 100,
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
  compute_support_metrics <- function(coefs, support, p_dim) {
    estimated_support <- which(abs(coefs) > 0)
    stats <- decentralized1bit::computeF1(estimated_support, support)
    fp <- stats$fp
    fn <- stats$fn
    tp <- length(intersect(estimated_support, support))
    tn <- p_dim - length(support) - fp
    fpr <- if ((fp + tn) == 0) 0 else fp / (fp + tn)
    c(
      exact_support = as.numeric(setequal(estimated_support, support)),
      recall = stats$recall,
      precision = stats$precision,
      tpr = stats$recall,
      fpr = fpr
    )
  }

  mean_support_metrics <- function(coef_matrix, support, p_dim) {
    metrics <- vapply(
      seq_len(ncol(coef_matrix)),
      function(idx) compute_support_metrics(coef_matrix[, idx], support, p_dim),
      FUN.VALUE = c(exact_support = 0, recall = 0, precision = 0, tpr = 0, fpr = 0)
    )
    rowMeans(metrics)
  }

  RNGkind("L'Ecuyer-CMRG")
  set.seed(42)
  random_seed <- fixRNGStream(.Random.seed, batch)
  .GlobalEnv$.Random.seed <- random_seed

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

  beta_pooled <- bic.glmnet(y = y, X = X, nlambda = nlambda)
  pooled_metrics <- compute_support_metrics(beta_pooled, support_true, p)

  B_init <- matrix(0, nrow = p, ncol = m)
  for (node in seq_len(m)) {
    local_index <- calN_j(n, node)
    B_init[, node] <- bic.glmnet(X = X[local_index, ], y = y[local_index], nlambda = nlambda)
  }
  local_metrics <- mean_support_metrics(B_init, support_true, p)

  B_average <- rowMeans(B_init)
  average_metrics <- compute_support_metrics(B_average, support_true, p)

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
  lr_metrics <- mean_support_metrics(B_lr, support_true, p)

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
  subgd_metrics <- mean_support_metrics(B_subgd, support_true, p)

  methods <- factor(SIM_METHOD_LEVELS, levels = SIM_METHOD_LEVELS)
  metrics_matrix <- rbind(pooled_metrics, local_metrics, average_metrics, subgd_metrics, lr_metrics)

  metrics_df <- data.frame(
    exact_support = metrics_matrix[, "exact_support"],
    recall = metrics_matrix[, "recall"],
    precision = metrics_matrix[, "precision"],
    tpr = metrics_matrix[, "tpr"],
    fpr = metrics_matrix[, "fpr"],
    methods = methods
  )

  rownames(metrics_df) <- NULL
  metrics_df
}
