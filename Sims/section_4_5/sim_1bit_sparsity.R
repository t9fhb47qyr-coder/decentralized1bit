library(rslurm)

source("Sims/sim_utils.R")
source("Sims/sim_single_simulation.R")

#' Experiment 4.5: vary the sparsity level of the ground-truth signal.
simulation_name <- "sparsity"
run_context <- initialize_simulation_run(simulation_name)
submission_plan <- default_submission_plan()
submission_plan$reps <- 100

sparsity_levels <- c(60) #c(seq(5, 30, by = 5), 40, 50)
params <- expand.grid(
  s = sparsity_levels,
  batch = seq_len(submission_plan$reps)
)

sjob <- run_slurm_simulation(
  f = sim_single_simulation,
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

summarise_slurm_results(
  sjob = sjob,
  params = params,
  varying_param = "s",
  output_dir = run_context$output_dir,
  jobname = run_context$jobname,
  simulation_name = simulation_name
)

save.image(file = file.path(run_context$output_dir, paste0(run_context$jobname, ".RData")))
