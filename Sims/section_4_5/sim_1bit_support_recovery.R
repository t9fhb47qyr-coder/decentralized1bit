library(rslurm)

source("Sims/sim_utils.R")
source("Sims/sim_support_recovery_single.R")

#' Experiment 4.5 (support recovery): vary sparsity and track support metrics.
simulation_name <- "sparsity_support_recovery"
run_context <- initialize_simulation_run(simulation_name)
submission_plan <- default_submission_plan()
#submission_plan$reps <- 100

sparsity_levels <- seq(5, 30, by = 5)
params <- expand.grid(
  s = sparsity_levels,
  batch = seq_len(submission_plan$reps)
)

sjob <- run_slurm_simulation(
  f = sim_single_support_recovery,
  params = params,
  jobname = run_context$jobname,
  submit = submission_plan$submit,
  nodes = 100,
  cpus_per_node = 2,
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

summarise_support_recovery_results(
  sjob = sjob,
  params = params,
  varying_param = "s",
  output_dir = run_context$output_dir,
  jobname = run_context$jobname,
  simulation_name = simulation_name
)


combined_csv <- file.path(run_context$output_dir, paste0(run_context$jobname, "_support_recovery.csv"))
if (file.exists(combined_csv)) {
  plot_support_recovery(combined_csv = combined_csv, varying_param = "s")
}

save.image(file = file.path(run_context$output_dir, paste0(run_context$jobname, ".RData")))

# Plot
