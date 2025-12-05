library(rslurm)

source("Sims/sim_utils.R")
source("Sims/sim_single_simulation.R")

#' App. A.1–A.2: vary the number of inner iterations for the decentralized solvers.
simulation_name <- "T_inner"
run_context <- initialize_simulation_run(simulation_name)
submission_plan <- default_submission_plan()

T_inner_values <- c(500, 1000, 1500)
params <- expand.grid(
  T_inner = T_inner_values,
  batch = seq_len(submission_plan$reps)
)

sjob <- run_slurm_simulation(
  f = sim_single_simulation,
  params = params,
  jobname = run_context$jobname,
  submit = submission_plan$submit,
  nodes = 100,
  cpus_per_node = 2,
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
  varying_param = "T_inner",
  output_dir = run_context$output_dir,
  jobname = run_context$jobname,
  simulation_name = simulation_name
)

save.image(file = file.path(run_context$output_dir, paste0(run_context$jobname, ".RData")))
