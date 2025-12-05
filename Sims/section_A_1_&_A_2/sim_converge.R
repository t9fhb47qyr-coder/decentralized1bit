library(decentralized1bit)
library(glmnet)
library(igraph)
library(ggplot2)
library(tidyr)

compute_mean_f1 <- function(coef_matrix, support_true) {
  column_f1 <- apply(
    coef_matrix,
    2,
    function(coefs) decentralized1bit::computeF1(which(abs(coefs) > 0), support_true)[[1]]
  )
  mean(unlist(column_f1))
}

build_local_estimates <- function(X, y, n, m, nlambda) {
  estimates <- matrix(0, nrow = ncol(X), ncol = m)
  for (node in seq_len(m)) {
    idx <- decentralized1bit::calN_j(n, node)
    estimates[, node] <- decentralized1bit::bic.glmnet(
      X = X[idx, ],
      y = y[idx],
      nlambda = nlambda
    )
  }
  estimates
}

extract_f1_history <- function(result_list, support_true) {
  vapply(
    result_list,
    compute_mean_f1,
    support_true = support_true,
    FUN.VALUE = numeric(1)
  )
}

plot_history <- function(data, y_label) {
  ggplot(data, aes(x = Iteration, y = Value, colour = Algorithm, linetype = Algorithm)) +
    geom_line(linewidth = 0.8, alpha = 0.85, show.legend = FALSE) +
    scale_color_manual(
      values = c("DeLR" = "#E69F00", "D-subGD" = "#56B4E9")
    ) +
    scale_linetype_manual(values = c("DeLR" = "solid", "D-subGD" = "dashed")) +
    scale_x_continuous(breaks = seq(0, max(data$Iteration), by = 100)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = y_label) +
    theme_minimal(base_size = 16)
}

config <- list(
  n = 100,
  m = 20,
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
  batch = 1
)

RNGkind("L'Ecuyer-CMRG")
set.seed(42)
.GlobalEnv$.Random.seed <- decentralized1bit::fixRNGStream(.Random.seed, config$batch)

data <- decentralized1bit::genData(
  N = config$m * config$n,
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

X <- data$X
y <- data$y
beta_true <- data$betaT / base::norm(data$betaT, "F")
support_true <- which(abs(beta_true) > 0)
adjacency_matrix <- as.matrix(igraph::as_adjacency_matrix(data$graph))

B_init <- build_local_estimates(X, y, config$n, config$m, config$nlambda)

delr <- decentralized1bit::deLR(
  X = X,
  y = y,
  adjacency_matrix = adjacency_matrix,
  B_init = B_init,
  betaT = beta_true,
  tau_penalty_factor = 10 / (3 * config$m),
  T_inner = config$T_inner,
  lambda_max = 0.01,
  nlambda = config$nlambda,
  quiet = config$quiet
)

dsubgd <- decentralized1bit::deSubGD_LS(
  X = X,
  y = y,
  adjacency_matrix = adjacency_matrix,
  B_init = B_init,
  betaT = beta_true,
  T = config$T_inner,
  nlambda = config$nlambda,
  lambda_max = 20,
  eps = 1e-7
)

iterations <- seq_len(config$T_inner)

error_history <- data.frame(
  Iteration = iterations,
  `DeLR` = delr$history$errors_inner,
  `D-subGD` = dsubgd$history$errors_inner
)
error_history <- tidyr::pivot_longer(
  error_history,
  cols = -Iteration,
  names_to = "Algorithm",
  values_to = "Value"
)

f1_history <- data.frame(
  Iteration = iterations,
  `DeLR` = extract_f1_history(delr$result, support_true),
  `D-subGD` = extract_f1_history(dsubgd$history$result, support_true)
)
f1_history <- tidyr::pivot_longer(
  f1_history,
  cols = -Iteration,
  names_to = "Algorithm",
  values_to = "Value"
)

error_plot <- plot_history(error_history, "L2 Error")
f1_plot <- plot_history(f1_history, "F1 Score")

print(error_plot)
print(f1_plot)
