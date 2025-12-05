library(rslurm)
library(decentralized1bit)
library(glmnet)
library(igraph)

source("Sims/sim_utils.R")
source("Sims/sim_single_simulation.R")

simulate_highdimension <- function(label, n, m, N, p, s, batch, quiet = TRUE) {
  result <- sim_single_simulation(
    n = n,
    m = m,
    N = N,
    p = p,
    s = s,
    batch = batch,
    quiet = quiet
  )

  result$label <- label
  result$batch <- batch
  result$n <- n
  result$m <- m
  result$N <- N
  result$p <- p
  result$s <- s
  existing_cols <- names(result)
  metadata_cols <- c("label", "batch", "n", "m", "N", "p", "s")
  other_cols <- existing_cols[!(existing_cols %in% metadata_cols)]
  result[, c(metadata_cols, other_cols)]
}

simulation_name <- "highdimension"
run_context <- initialize_simulation_run(simulation_name)
submission_plan <- default_submission_plan(local_reps = 5L, cluster_reps = 100L)

scenario_grid <- data.frame(
  label = c("N200_p600", "N300_p600", "N400_p600", "N500_p600", "N600_p600"),
  n = rep(100L, 5),
  m = 5:9,
  N = seq(200L, 600L, by = 100L),
  p = rep(600L, 5),
  s = rep(10L, 5),
  stringsAsFactors = FALSE
)

params <- scenario_grid[rep(seq_len(nrow(scenario_grid)), each = submission_plan$reps), ]
params$batch <- rep(seq_len(submission_plan$reps), times = nrow(scenario_grid))
params$label <- factor(params$label, levels = scenario_grid$label)

sjob <- run_slurm_simulation(
  f = simulate_highdimension,
  params = params,
  jobname = run_context$jobname,
  submit = submission_plan$submit,
  nodes = 50,
  cpus_per_node = 1,
  global_objects = ls(),
  constant_args = list(quiet = TRUE)
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
