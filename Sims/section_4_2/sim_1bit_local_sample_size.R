library(rslurm)

source("Sims/sim_utils.R")
source("Sims/sim_single_simulation.R")

#' Experiment 4.2: vary the local sample size while fixing other parameters.
simulation_name <- "local_sample_size"
run_context <- initialize_simulation_run(simulation_name)
submission_plan <- default_submission_plan()

local_sample_sizes <- seq(50, 300, by = 50)
params <- expand.grid(
  n = local_sample_sizes,
  batch = seq_len(submission_plan$reps)
)

sjob <- run_slurm_simulation(
  f = sim_single_simulation,
  params = params,
  jobname = run_context$jobname,
  submit = submission_plan$submit,
  nodes = 20,
  cpus_per_node = 1,
  global_objects = ls(),
  constant_args = list(m = 20, s = 10)
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
  varying_param = "n",
  output_dir = run_context$output_dir,
  jobname = run_context$jobname,
  simulation_name = simulation_name
)

save.image(file = file.path(run_context$output_dir, paste0(run_context$jobname, ".RData")))
