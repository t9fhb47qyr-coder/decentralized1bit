#' Helper utilities for running and summarizing decentralized 1-bit simulations.
#'
#' These helpers centralize repeated patterns that appear across the scripts in
#' `Sims/`. They provide small building blocks for configuring Slurm jobs,
#' persisting outputs, and producing tidy summary tables.

SIM_METHOD_LEVELS <- c("Pooled", "Local", "Avg", "D-subGD", "deLR")

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}


#' Check whether the current session is running on macOS.
#'
#' @return Logical scalar.
is_running_on_mac <- function() {
  identical(Sys.info()[["sysname"]], "Darwin")
}

#' Compute default submission settings based on the host machine.
#'
#' @param local_reps Number of Monte Carlo repetitions when not submitting jobs.
#' @param cluster_reps Number of repetitions when submitting to Slurm.
#'
#' @return A list with elements `submit` and `reps`.
default_submission_plan <- function(local_reps = 2L, cluster_reps = 200L) {
  if (is_running_on_mac()) {
    list(submit = FALSE, reps = local_reps)
  } else {
    list(submit = TRUE, reps = cluster_reps)
  }
}

#' Create a timestamped job name and ensure the output directory exists.
#'
#' @param simulation_name Short identifier (e.g., "local_sample_size").
#' @param output_dir Directory where artifacts should be stored.
#'
#' @return A list containing `jobname`, `timestamp`, and `output_dir`.
initialize_simulation_run <- function(simulation_name, output_dir = "Output") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  jobname <- paste(simulation_name, timestamp, sep = "_")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  list(jobname = jobname, timestamp = timestamp, output_dir = output_dir)
}

#' Wrapper around `rslurm::slurm_apply()` with clearer argument structure.
#'
#' @param f Simulation function to execute.
#' @param params Data frame describing the parameter grid.
#' @param jobname Identifier produced by `initialize_simulation_run()`.
#' @param submit Logical flag passed to `slurm_apply()`.
#' @param nodes Number of nodes to request.
#' @param cpus_per_node Number of CPUs per node to request.
#' @param global_objects Character vector of objects to export to Slurm workers.
#' @param slurm_options Named list forwarded to `slurm_options` argument.
#' @param constant_args Named list forwarded to the simulation function through
#'   `...`.
#'
#' @return The Slurm job object returned by `slurm_apply()`.
run_slurm_simulation <- function(f,
                                 params,
                                 jobname,
                                 submit,
                                 nodes = 20L,
                                 cpus_per_node = 1L,
                                 global_objects = NULL,
                                 slurm_options = list(),
                                 constant_args = list()) {
  apply_args <- list(
    f = f,
    params = params,
    jobname = jobname,
    submit = submit,
    nodes = nodes,
    cpus_per_node = cpus_per_node,
    slurm_options = slurm_options
  )

  if (!is.null(global_objects)) {
    apply_args$global_objects <- global_objects
  }

  apply_args <- c(apply_args, constant_args)
  do.call(rslurm::slurm_apply, apply_args)
}

#' Persist a Slurm job handle using both timestamped and canonical filenames.
#'
#' @param sjob Slurm job object.
#' @param output_dir Destination directory.
#' @param jobname Timestamped job name.
#' @param simulation_name Canonical simulation identifier.
save_slurm_handles <- function(sjob, output_dir, jobname, simulation_name) {
  saveRDS(sjob, file = file.path(output_dir, paste0(jobname, ".RDS")))
  saveRDS(sjob, file = file.path(output_dir, paste0(simulation_name, ".RDS")))
}

#' Summarise errors and F1 scores across Monte Carlo repetitions.
#'
#' @param sjob Slurm job object.
#' @param params Parameter grid used to launch the job.
#' @param varying_param Name of the varying parameter column in `params`.
#' @param output_dir Directory where tables should be written.
#' @param jobname Timestamped job identifier.
#' @param simulation_name Canonical simulation identifier.
#' @param method_levels Factor levels for the `methods` column.
#' @param digits List with optional entries `errors`, `f1s`, and `combined`
#'   defining the digits vector supplied to `xtable::xtable()`.
#' @param wait Forwarded to `rslurm::get_slurm_out()`.
#' @param print_tables Logical; when `TRUE`, formatted tables are printed.
#'
#' @return Invisibly returns a list with wide-format data frames for errors,
#'   F1 scores, and the combined table.
summarise_slurm_results <- function(sjob,
                                    params,
                                    varying_param,
                                    output_dir,
                                    jobname,
                                    simulation_name,
                                    method_levels = SIM_METHOD_LEVELS,
                                    digits = list(),
                                    wait = TRUE,
                                    print_tables = TRUE) {
  result <- rslurm::get_slurm_out(sjob, outtype = "table", wait = wait)
  if (!nrow(result)) {
    stop("No results returned by get_slurm_out(). Did the Slurm job finish?")
  }

  result$methods <- factor(result$methods, levels = method_levels)
  n_methods <- length(method_levels)
  expanded_params <- params[rep(seq_len(nrow(params)), each = n_methods), , drop = FALSE]
  combined <- cbind(expanded_params, result)

  build_summary <- function(metric) {
    agg_formula <- stats::as.formula(sprintf("%s ~ methods + %s", metric, varying_param))
    summary_df <- stats::aggregate(agg_formula, data = combined, FUN = mean)
    reshape(summary_df,
            idvar = varying_param,
            timevar = "methods",
            direction = "wide")
  }

  rename_metric_columns <- function(df, metric) {
    new_names <- gsub(sprintf("%s\\.", metric), paste0(metric, "_"), names(df), fixed = FALSE)
    stats::setNames(df, new_names)
  }

  errors_wide <- rename_metric_columns(build_summary("errors"), "errors")
  f1s_wide <- rename_metric_columns(build_summary("f1s"), "f1s")

  combined_table <- merge(errors_wide, f1s_wide, by = varying_param, sort = TRUE)

  errors_path <- file.path(output_dir, paste0(jobname, "_errors.csv"))
  f1s_path <- file.path(output_dir, paste0(jobname, "_f1s.csv"))
  combined_path <- file.path(output_dir, paste0(jobname, "_errors_f1s.csv"))

  utils::write.csv(errors_wide, errors_path, row.names = FALSE)
  utils::write.csv(f1s_wide, f1s_path, row.names = FALSE)
  utils::write.csv(combined_table, combined_path, row.names = FALSE)

  if (print_tables) {
    if (!requireNamespace("xtable", quietly = TRUE)) {
      message("Package 'xtable' not installed; skipping LaTeX table output.")
    } else {
      errors_digits <- digits$errors %||% c(0, 0, rep(4, ncol(errors_wide) - 1))
      f1_digits <- digits$f1s %||% c(0, 0, rep(4, ncol(f1s_wide) - 1))
      combined_digits <- digits$combined %||% c(0, 0, rep(4, ncol(combined_table) - 1))

      print(xtable::xtable(errors_wide, digits = errors_digits), include.rownames = FALSE)
      print(xtable::xtable(f1s_wide, digits = f1_digits), include.rownames = FALSE)
      print(xtable::xtable(combined_table, digits = combined_digits), include.rownames = FALSE)
    }
  }

  invisible(list(
    errors = errors_wide,
    f1s = f1s_wide,
    combined = combined_table,
    output_paths = list(errors = errors_path, f1s = f1s_path, combined = combined_path)
  ))
}

#' Summarise support recovery metrics (exact recovery, recall, precision).
#'
#' @param sjob Slurm job object.
#' @param params Parameter grid used to launch the job.
#' @param varying_param Name of the varying parameter column in `params`.
#' @param output_dir Directory where tables should be written.
#' @param jobname Timestamped job identifier.
#' @param simulation_name Canonical simulation identifier.
#' @param method_levels Factor levels for the `methods` column.
#' @param digits Optional list with entries `exact_support`, `recall`,
#'   `precision`, and `combined` controlling xtable rounding.
#' @param wait Forwarded to `rslurm::get_slurm_out()`.
#' @param print_tables Logical; when `TRUE`, formatted tables are printed.
#'
#' @return Invisibly returns a list with wide-format data frames for each metric
#'   and the combined table.
summarise_support_recovery_results <- function(sjob,
                                               params,
                                               varying_param,
                                               output_dir,
                                               jobname,
                                               simulation_name,
                                               method_levels = SIM_METHOD_LEVELS,
                                               digits = list(),
                                               wait = TRUE,
                                               print_tables = TRUE) {
  result <- rslurm::get_slurm_out(sjob, outtype = "table", wait = wait)
  if (!nrow(result)) {
    stop("No results returned by get_slurm_out(). Did the Slurm job finish?")
  }

  metrics <- c("exact_support", "recall", "precision", "tpr", "fpr")

  result$methods <- factor(result$methods, levels = method_levels)
  n_methods <- length(method_levels)
  expanded_params <- params[rep(seq_len(nrow(params)), each = n_methods), , drop = FALSE]
  combined <- cbind(expanded_params, result)

  build_summary <- function(metric) {
    agg_formula <- stats::as.formula(sprintf("%s ~ methods + %s", metric, varying_param))
    summary_df <- stats::aggregate(agg_formula, data = combined, FUN = mean)
    reshape(summary_df,
            idvar = varying_param,
            timevar = "methods",
            direction = "wide")
  }

  rename_metric_columns <- function(df, metric) {
    new_names <- gsub(sprintf("%s\\.", metric), paste0(metric, "_"), names(df), fixed = FALSE)
    stats::setNames(df, new_names)
  }

  summary_tables <- lapply(metrics, function(metric) rename_metric_columns(build_summary(metric), metric))
  names(summary_tables) <- metrics

  combined_table <- Reduce(function(x, y) merge(x, y, by = varying_param, sort = TRUE), summary_tables)

  output_paths <- vapply(metrics, function(metric) {
    path <- file.path(output_dir, paste0(jobname, "_", metric, ".csv"))
    utils::write.csv(summary_tables[[metric]], path, row.names = FALSE)
    path
  }, FUN.VALUE = character(1))

  combined_path <- file.path(output_dir, paste0(jobname, "_support_recovery.csv"))
  utils::write.csv(combined_table, combined_path, row.names = FALSE)

  if (print_tables) {
    if (!requireNamespace("xtable", quietly = TRUE)) {
      message("Package 'xtable' not installed; skipping LaTeX table output.")
    } else {
      digits_defaults <- c(0, 0, rep(4, ncol(combined_table) - 1))
      print_digits <- function(metric) digits[[metric]] %||% digits_defaults
      lapply(metrics, function(metric) {
        print(
          xtable::xtable(summary_tables[[metric]], digits = print_digits(metric)),
          include.rownames = FALSE
        )
      })
      print(xtable::xtable(combined_table, digits = digits$combined %||% digits_defaults),
            include.rownames = FALSE)
    }
  }

  invisible(list(
    summary_tables = summary_tables,
    combined = combined_table,
    output_paths = c(as.list(output_paths), combined = combined_path)
  ))
}

#' Plot support recovery metrics from the combined CSV produced by
#' `summarise_support_recovery_results()`.
#'
#' @param combined_csv Path to the combined support recovery CSV.
#' @param varying_param Name of the varying parameter column (default: "s").
#' @param metric_prefixes Metric prefixes to plot.
#' @param output_path Optional path to save the PNG plot. Defaults to
#'   `<combined_csv>_plot.png`.
#'
#' @return The path to the saved plot (invisibly).
plot_support_recovery <- function(combined_csv,
                                  varying_param = "s",
                                  metric_prefixes = c("exact_support", "recall", "precision", "tpr", "fpr"),
                                  output_path = NULL) {
  if (!file.exists(combined_csv)) {
    stop("Combined support recovery CSV not found: ", combined_csv)
  }

  suppressPackageStartupMessages({
    library(ggplot2)
    library(tidyr)
    library(dplyr)
  })

  df <- read.csv(combined_csv, stringsAsFactors = FALSE)

  build_long <- function(metric) {
    cols <- grep(paste0("^", metric, "_"), names(df), value = TRUE)
    if (!length(cols)) return(NULL)
    df_metric <- df[, c(varying_param, cols), drop = FALSE]
    long <- tidyr::pivot_longer(
      df_metric,
      cols = cols,
      names_to = "method",
      values_to = "value"
    )
    long$method <- sub(paste0("^", metric, "_"), "", long$method)
    long$method <- factor(long$method, levels = SIM_METHOD_LEVELS)
    long$metric <- metric
    long[, c(varying_param, "method", "metric", "value")]
  }

  long_df <- bind_rows(lapply(metric_prefixes, build_long))

  if (!nrow(long_df)) {
    stop("No metrics found in the provided CSV. Expected prefixes: ",
         paste(metric_prefixes, collapse = ", "))
  }

  plot_obj <- ggplot(long_df, aes_string(x = varying_param, y = "value", color = "method")) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    labs(
      title = "Support Recovery Metrics",
      x = varying_param,
      y = "Metric value",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )

  if (is.null(output_path)) {
    output_path <- sub("\\.csv$", "_plot.png", combined_csv)
  }
  ggsave(output_path, plot_obj, width = 9, height = 6, dpi = 300)
  message("Saved plot to: ", output_path)
  invisible(output_path)
}
