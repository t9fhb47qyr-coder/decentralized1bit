library(rslurm)
library(decentralized1bit)
library(glmnet)
library(igraph)
library(MASS)
library(copula)

source("Sims/sim_utils.R")
source("Sims/sim_single_simulation.R")

compute_mean_f1 <- function(coef_matrix, support_true) {
  column_f1 <- apply(
    coef_matrix,
    2,
    function(coefs) decentralized1bit::computeF1(which(abs(coefs) > 0), support_true)[[1]]
  )
  mean(unlist(column_f1))
}

evaluate_estimators <- function(X,
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

  beta_pooled <- decentralized1bit::bic.glmnet(y = y, X = X, nlambda = config$nlambda)
  error_pooled <- decentralized1bit::computeError(beta_pooled, beta_unit)
  f1_pooled <- decentralized1bit::computeF1(which(abs(beta_pooled) > 0), support_true)[[1]]

  error_local <- decentralized1bit::computeError(B_init, beta_unit)
  f1_local <- compute_mean_f1(B_init, support_true)

  B_average <- rowMeans(B_init)
  error_average <- decentralized1bit::computeError(B_average, beta_unit)
  f1_average <- decentralized1bit::computeF1(which(abs(B_average) > 0), support_true)[[1]]

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
    errors = c(error_pooled, error_local, error_average, error_subgd, error_lr),
    f1s = c(f1_pooled, f1_local, f1_average, f1_subgd, f1_lr),
    methods = factor(SIM_METHOD_LEVELS, levels = SIM_METHOD_LEVELS)
  )
}

generate_noise <- function(n, type) {
  switch(
    type,
    Cauchy = stats::rcauchy(n),
    Normal = stats::rnorm(n),
    T2 = stats::rt(n, df = 2),
    stop("Unsupported noise type: ", type)
  )
}

generate_connected_graph <- function(m, pc) {
  graph <- igraph::sample_gnp(m, pc)
  while (!igraph::is_connected(graph)) {
    graph <- igraph::sample_gnp(m, pc)
  }
  graph
}

build_beta <- function(p, s) {
  beta <- rep(0, p)
  beta[sample.int(p, s)] <- 2 * stats::rbinom(s, 1, 0.5) - 1
  beta
}

generate_baseline_dataset <- function(config) {
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

generate_heterogeneous_dataset <- function(config) {
  proportion <- config$proportion %||% 0.4
  n <- config$N / config$m
  if (abs(n - round(n)) > .Machine$double.eps^0.5) {
    stop("N must be a multiple of m for heterogeneous generator.")
  }

  beta <- build_beta(config$p, config$s)
  graph <- generate_connected_graph(config$m, config$pc)

  Sigma_base <- config$sigma2 * stats::toeplitz(config$rho^seq(0, config$p - 1))
  Sigma_alt <- 3 * stats::toeplitz(0.3^seq(0, config$p - 1))

  case <- config$case
  X <- matrix(0, nrow = config$N, ncol = config$p)
  y <- numeric(config$N)

  if (case == 0) {
    n_alt <- round(config$N * proportion)
    n_base <- config$N - n_alt
    X[seq_len(n_base), ] <- MASS::mvrnorm(n_base, mu = rep(0, config$p), Sigma = Sigma_base)
    X[(n_base + 1):config$N, ] <- MASS::mvrnorm(n_alt, mu = rep(0, config$p), Sigma = Sigma_alt)
    noise <- generate_noise(config$N, config$noise_type)
    y <- sign(X %*% beta + noise)
    flips <- sample.int(config$N, size = floor(config$N * config$p_flip))
    y[flips] <- -y[flips]
  } else if (case == 1) {
    sigma2_set <- c(1, 3)
    rho_set <- c(0.1, 0.3)
    noise <- generate_noise(config$N, config$noise_type)
    for (node in seq_len(config$m)) {
      idx <- decentralized1bit::calN_j(n, node)
      Sigma_node <- sample(sigma2_set, 1) * stats::toeplitz(sample(rho_set, 1)^seq(0, config$p - 1))
      X[idx, ] <- MASS::mvrnorm(n, mu = rep(0, config$p), Sigma = Sigma_node)
      y[idx] <- sign(X[idx, ] %*% beta + noise[idx])
    }
    flips <- sample.int(config$N, size = floor(config$N * config$p_flip))
    y[flips] <- -y[flips]
  } else if (case == 2) {
    X <- MASS::mvrnorm(config$N, mu = rep(0, config$p), Sigma = Sigma_base)
    noise_type_ids <- sample(1:3, config$m, replace = TRUE)
    noise <- rep(0, config$N)
    for (node in seq_len(config$m)) {
      idx <- decentralized1bit::calN_j(n, node)
      noise[idx] <- switch(
        noise_type_ids[node],
        generate_noise(n, "Cauchy"),
        generate_noise(n, "Normal"),
        generate_noise(n, "T2")
      )
    }
    y <- sign(X %*% beta + noise)
    flips <- sample.int(config$N, size = floor(config$N * config$p_flip))
    y[flips] <- -y[flips]
  } else if (case == 3) {
    X <- MASS::mvrnorm(config$N, mu = rep(0, config$p), Sigma = Sigma_base)
    noise <- generate_noise(config$N, config$noise_type)
    y <- sign(X %*% beta + noise)
    p_flip_set <- c(config$p_flip, config$p_flip * 2)
    for (node in seq_len(config$m)) {
      node_idx <- decentralized1bit::calN_j(n, node)
      flips <- sample(node_idx, size = floor(length(node_idx) * sample(p_flip_set, 1)))
      y[flips] <- -y[flips]
    }
  } else {
    stop("Unsupported heterogeneity case: ", case)
  }

  list(X = X, y = y, betaT = matrix(beta, ncol = 1), graph = graph)
}

run_single_iteration <- function(generator, config, extra_args = list(), batch_id = 1L) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(42)
  .GlobalEnv$.Random.seed <- decentralized1bit::fixRNGStream(.Random.seed, batch_id)

  generator_config <- modifyList(config, extra_args)
  data <- generator(generator_config)
  evaluate_estimators(data$X, data$y, data$betaT, data$graph, generator_config)
}

simulate_heterogeneity <- function(label, generator, proportion, batch) {
  generator_fun <- get(generator, mode = "function")
  extra_args <- list()
  if (!is.na(proportion)) {
    extra_args$proportion <- proportion
  }

  result <- run_single_iteration(
    generator = generator_fun,
    config = base_config,
    extra_args = extra_args,
    batch_id = batch
  )

  result$label <- label
  result$batch <- batch
  result$generator <- generator
  result$proportion <- proportion
  existing_cols <- names(result)
  metadata_cols <- c("label", "batch", "generator", "proportion")
  other_cols <- existing_cols[!(existing_cols %in% metadata_cols)]
  result[, c(metadata_cols, other_cols)]
}

base_config <- list(
  n = 100,
  m = 20,
  N = 2000,
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
  replications = 100L
)

simulation_name <- "heterogeneity_supplement"
run_context <- initialize_simulation_run(simulation_name)
submission_plan <- default_submission_plan(local_reps = 5L, cluster_reps = 100L)
base_config$replications <- submission_plan$reps

scenario_grid <- data.frame(
  label = c("baseline", "heter_prop_0_2", "heter_prop_0_4", "heter_prop_0_6", "heter_prop_0_8"),
  generator = c(
    "generate_baseline_dataset",
    rep("generate_heterogeneous_dataset", 4)
  ),
  proportion = c(NA, 0.2, 0.4, 0.6, 0.8),
  stringsAsFactors = FALSE
)

params <- scenario_grid[rep(seq_len(nrow(scenario_grid)), each = submission_plan$reps), ]
params$batch <- rep(seq_len(submission_plan$reps), times = nrow(scenario_grid))
params$label <- factor(params$label, levels = scenario_grid$label)

sjob <- run_slurm_simulation(
  f = simulate_heterogeneity,
  params = params,
  jobname = run_context$jobname,
  submit = submission_plan$submit,
  nodes = 50,
  cpus_per_node = 1,
  global_objects = ls()
)

if (!submission_plan$submit) {
  rslurm::local_slurm_array(sjob)
}

save_slurm_handles(
  sjob = sjob,
  output_dir = run_context$output_dir,
  jobname = run_context$jobname,
  simulation_name = simulation_name
)

status <- rslurm::get_job_status(sjob)$queue
message("Current Slurm job status: ", status)

summarise_slurm_results(
  sjob = sjob,
  params = params,
  varying_param = "label",
  output_dir = run_context$output_dir,
  jobname = run_context$jobname,
  simulation_name = simulation_name
)

slurm_output <- rslurm::get_slurm_out(sjob, type = "table", wait = TRUE)
scenario_results <- setNames(
  lapply(scenario_grid$label, function(lbl) subset(slurm_output, label == lbl)),
  scenario_grid$label
)

save(
  list = c("scenario_results", "simulation_name"),
  file = file.path(run_context$output_dir, paste0(simulation_name, ".RData"))
)
