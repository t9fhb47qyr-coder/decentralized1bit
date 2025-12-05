#!/usr/bin/env Rscript

# Plot support recovery results (exact recovery, recall/precision, TPR/FPR)
# produced by `sim_1bit_support_recovery.R`.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript Sims/section_4_5/plot_support_recovery.R <combined_support_recovery_csv> [output_png]")
}

combined_csv <- args[[1]]
output_png <- if (length(args) >= 2) args[[2]] else NULL

source("Sims/sim_utils.R")

plot_support_recovery(combined_csv = combined_csv, output_path = output_png)
