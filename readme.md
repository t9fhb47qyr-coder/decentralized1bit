# Efficient Decoding from Heterogeneous 1-Bit Compressive Measurements over Networks

This repository accompanies the paper *Efficient Decoding from Heterogeneous 1-Bit Compressive Measurements over Networks*. It includes the source code of the `decentralized1bit` R package together with simulation and plotting scripts that reproduce all empirical results in the main manuscript and supplementary material.

## Getting Started

- Install [R](https://www.r-project.org/) (tested with R 4.2+) and a working C/C++ toolchain for building `Rcpp` code.
- Install the package dependencies. From an interactive R session:

```r
install.packages(c(
  "hqreg", "Rcpp", "RcppArmadillo",
  "glmnet", "igraph", "rslurm",
  "xtable", "reshape2", "data.table", "tidyr",
  "ggplot2", "ggthemes"
))
```

- Install the local package:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_local(".")
```

The installation step builds the `decentralized1bit` package and makes helper functions such as `deLR`, `deSubGD_LS`, and `genData` available to the simulation scripts.

## Repository Layout

- `R/`, `src/`, `man/`: package sources and documentation.
- `Sims/`: scripts used to generate simulation results and figures.
- `Output/`: default destination for intermediate results, tables, and figures created by the scripts.
- `realdata/`: resources for the real-world experiment (Section B of the paper).
- `decentralized1bit.Rproj`: optional RStudio project file.

## Reproducing the Manuscript & Supplement

Each script can be run from the project root with `Rscript path/to/script.R`. On macOS the scripts default to `submit = FALSE` and use `rslurm::local_slurm_array()` with a small number of repetitions (`Nreps = 2`) so that the code runs locally. To reproduce the full-scale experiments used in the paper, edit the script(s) to set `submit = TRUE` and increase `Nreps` before running on a Slurm-enabled cluster.

| Manuscript section | Description | Script(s) | Notes |
| --- | --- | --- | --- |
| §4.2 | Effect of local sample size | `Sims/section_4_2/sim_1bit_local_sample_size.R` | Writes CSV tables and `*.RDS` snapshots under `Output/`. |
| §4.3 | Proportion of common support | `Sims/section_4_3/sim_1bit_pc.R` | Adjust the `pc` vector to explore other levels. |
| §4.4 | Additional heterogeneity study | `Sims/section_A_6/sim_heter_supp.R` |  | 
| §4.5 | Dimensionality, sparsity, sign flips | `Sims/section_4_5/sim_1bit_dimension.R`, `sim_1bit_sparsity.R`, `sim_1bit_p_flip.R` | Batch grids are defined at the top of each script. |
| App. §A.1–A.2 | Convergence diagnostics; inner-iteration sensitivity | `Sims/section_A_1_&_A_2/sim_converge.R`, `Sims/section_A_1_&_A_2/sim_T_inner.R` | Produces diagnostics saved in `Output/`. |
| App. §A.3 | High-dimensional regime | `Sims/section_A_3/sim_highdimension.R` | |
| App. §A.4 | Non-Gaussian noise | `Sims/section_A_4/sim_nongaussian.R`, `Sims/section_A_4/single_Nongaussian.R` | Use the single-run script to inspect an individual setting. |
| App. §A.5 | Single simulation trace | `Sims/sim_single_simulation.R` | Standalone execution returns a tibble with error metrics. |
| App. §A.6 | Heterogeneous networks | `Sims/section_4_4/sim_1bit_heter1.R`, `sim_1bit_heter2.R`, `sim_1bit_heter3.R` | Each script targets a different heterogeneity scenario. |
| Section B | Real EEG data application | `realdata/realdata.R` | Loads `realdata/realdata.RData` and saves figures as PDFs. |
| Figs. (main) | Publication-ready plots | `Sims/plots/ss_plot.R` | Assumes the relevant outputs already exist in `Output/`. |
| Figs. (supp.) | Supplementary plots | `Sims/plots/ss_plot_supp.R` | Writes figures to `Output/figs/`. |

All scripts assume the working directory is the repository root so relative paths resolve correctly. The generated artifacts inherit time-stamped file names; you can tweak `outputdir` at the top of each file if you prefer a different destination.

## Real Data Workflow

1. Generate simulation outputs required by the real-data plots (see table above).
2. Run `Rscript realdata/realdata.R`. The script uses the pre-packaged `realdata.RData` file and produces the boxplots reported in the paper (e.g., `realdata/positive.pdf`, `realdata/neutral.pdf`).

The plotting script depends on `ggplot2`, `ggthemes`, and `reshape2`.

## Troubleshooting

- Ensure `rslurm` is configured to point to your Slurm scheduler. When running locally, keep `submit = FALSE`.
- When executing on a cluster, verify that `Output/` is writable and shared across nodes if needed.
- Re-run `remotes::install_local(".")` after modifying the C++ sources under `src/`.

## Citation

If you use this code, please cite the accompanying paper:

> C. Chen, L. Zhu, and Z. Zhu, “Efficient Decoding from Heterogeneous 1-Bit Compressive Measurements over Networks,” 2025.

## License

For commercial or redistribution questions, please contact the authors.

