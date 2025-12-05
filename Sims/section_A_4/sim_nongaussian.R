library(rslurm)
library(decentralized1bit)
library(glmnet)
library(igraph)
library(copula)

source("Sims/sim_utils.R")
source("Sims/section_A_4/single_Nongaussian.R")

simulate_nongaussian <- function(label, n, m, N, p, s, batch, quiet = TRUE) {
  result <- single_sim_non(
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

simulation_name <- "nongaussian"
run_context <- initialize_simulation_run(simulation_name)
submission_plan <- default_submission_plan(local_reps = 5L, cluster_reps = 100L)

scenario_grid <- data.frame(
  label = c("N1000", "N2000", "N3000", "N4000", "N5000", "N6000"),
  n = seq(50L, 300L, by = 50L),
  m = rep(20L, 6),
  N = seq(1000L, 6000L, by = 1000L),
  p = rep(100L, 6),
  s = rep(10L, 6),
  stringsAsFactors = FALSE
)

params <- scenario_grid[rep(seq_len(nrow(scenario_grid)), each = submission_plan$reps), ]
params$batch <- rep(seq_len(submission_plan$reps), times = nrow(scenario_grid))
params$label <- factor(params$label, levels = scenario_grid$label)

sjob <- run_slurm_simulation(
  f = simulate_nongaussian,
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

slurm_output <- rslurm::get_slurm_out(sjob, type = "table", wait = TRUE)
scenario_results <- setNames(
  lapply(scenario_grid$label, function(lbl) subset(slurm_output, label == lbl)),
  scenario_grid$label
)

metric_columns <- setdiff(names(slurm_output), c("label", "batch", "n", "m", "N", "p", "s"))
summary_table <- aggregate(
  slurm_output[, metric_columns],
  by = list(label = slurm_output$label),
  FUN = mean
)

utils::write.csv(
  summary_table,
  file = file.path(run_context$output_dir, paste0(run_context$jobname, "_summary.csv")),
  row.names = FALSE
)

save(
  list = c("scenario_results", "simulation_name"),
  file = file.path(run_context$output_dir, paste0(simulation_name, ".RData"))
)
