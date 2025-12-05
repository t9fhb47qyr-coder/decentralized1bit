library(ggplot2)
library(latex2exp)

output_dir <- file.path("Output", "plots", "supplement")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

build_long_form <- function(matrix_values, x_values, method_ids) {
  stopifnot(ncol(matrix_values) == length(method_ids))
  data.frame(
    x = rep(x_values, times = length(method_ids)),
    method = rep(method_ids, each = length(x_values)),
    value = as.vector(matrix_values)
  )
}

plot_metric <- function(data,
                        x_label,
                        y_label,
                        x_breaks,
                        x_limits,
                        y_breaks,
                        y_limits,
                        shapes,
                        linetypes,
                        colors,
                        file_name) {
  plot_obj <- ggplot(
    data = data,
    mapping = aes(x = x, y = value, colour = method, shape = method, linetype = method)
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = linetypes) +
    scale_colour_manual(values = colors) +
    scale_x_continuous(name = x_label, breaks = x_breaks, limits = x_limits) +
    scale_y_continuous(name = y_label, breaks = y_breaks, limits = y_limits) +
    theme_minimal(base_size = 18) +
    theme(legend.position = "none")

  ggsave(
    filename = file_name,
    plot = plot_obj,
    path = output_dir,
    width = 8,
    height = 6
  )
}

plot_pair <- function(error_matrix,
                      score_matrix,
                      x_values,
                      method_ids,
                      shapes,
                      linetypes,
                      colors,
                      x_label,
                      x_breaks,
                      x_limits,
                      error_breaks,
                      error_limits,
                      score_breaks,
                      score_limits,
                      error_file,
                      score_file) {
  error_df <- build_long_form(error_matrix, x_values, method_ids)
  score_df <- build_long_form(score_matrix, x_values, method_ids)

  plot_metric(
    data = error_df,
    x_label = x_label,
    y_label = TeX("$L_2$-error"),
    x_breaks = x_breaks,
    x_limits = x_limits,
    y_breaks = error_breaks,
    y_limits = error_limits,
    shapes = shapes,
    linetypes = linetypes,
    colors = colors,
    file_name = error_file
  )

  plot_metric(
    data = score_df,
    x_label = x_label,
    y_label = TeX("$F_1$-score"),
    x_breaks = x_breaks,
    x_limits = x_limits,
    y_breaks = score_breaks,
    y_limits = score_limits,
    shapes = shapes,
    linetypes = linetypes,
    colors = colors,
    file_name = score_file
  )
}

standard_methods <- c("A", "B", "C", "D", "E")
standard_shapes <- c(A = 15, B = 16, C = 17, D = 18, E = 6)
standard_linetypes <- c(A = "solid", B = "solid", C = "dotted", D = "dashed", E = "dashed")
standard_colors <- c(
  A = "#1f77b4",
  B = "#ff7f0e",
  C = "#2ca02c",
  D = "#d62728",
  E = "#9467bd"
)

# High-dimension (setting 1)
error_high_1 <- matrix(
  c(
    0.2995, 1.0410, 0.8526, 0.9071, 0.4972,
    0.2537, 1.0468, 0.9660, 0.8582, 0.3684,
    0.2256, 1.0446, 0.9697, 0.8075, 0.2761,
    0.1994, 1.0404, 0.9689, 0.7839, 0.2276,
    0.1756, 1.0405, 0.9393, 0.7701, 0.2106
  ),
  byrow = TRUE,
  ncol = length(standard_methods)
)

score_high_1 <- matrix(
  c(
    0.8916, 0.0590, 0.2487, 0.0198, 0.2557,
    0.9214, 0.0506, 0.2931, 0.0198, 0.3732,
    0.9359, 0.0506, 0.2951, 0.0198, 0.4504,
    0.9277, 0.0627, 0.3640, 0.0198, 0.5915,
    0.9313, 0.0538, 0.3545, 0.0198, 0.7027
  ),
  byrow = TRUE,
  ncol = length(standard_methods)
)

plot_pair(
  error_matrix = error_high_1,
  score_matrix = score_high_1,
  x_values = 2:6,
  method_ids = standard_methods,
  shapes = standard_shapes,
  linetypes = standard_linetypes,
  colors = standard_colors,
  x_label = "Number of nodes",
  x_breaks = 2:6,
  x_limits = c(2, 6),
  error_breaks = seq(0, 1.25, by = 0.25),
  error_limits = c(0, 1.25),
  score_breaks = seq(0, 1, by = 0.2),
  score_limits = c(0, 1),
  error_file = "high_dimension_1_error.pdf",
  score_file = "high_dimension_1_score.pdf"
)

# High-dimension (setting 2)
error_high_2 <- matrix(
  c(
    0.3918, 0.9693, 0.7643, 1.1124, 0.7637,
    0.2617, 0.9650, 0.7618, 0.8310, 0.5327,
    0.1935, 0.9824, 0.7421, 0.7308, 0.3342,
    0.1643, 0.9669, 0.6415, 0.6836, 0.2125
  ),
  byrow = TRUE,
  ncol = length(standard_methods)
)

score_high_2 <- matrix(
  c(
    0.8787, 0.3278, 0.5398, 0.0198, 0.1480,
    0.9243, 0.3368, 0.6456, 0.0198, 0.2910,
    0.9295, 0.3060, 0.7111, 0.0198, 0.4639,
    0.9335, 0.3227, 0.7686, 0.0198, 0.5903
  ),
  byrow = TRUE,
  ncol = length(standard_methods)
)

plot_pair(
  error_matrix = error_high_2,
  score_matrix = score_high_2,
  x_values = 2:5,
  method_ids = standard_methods,
  shapes = standard_shapes,
  linetypes = standard_linetypes,
  colors = standard_colors,
  x_label = "Number of nodes",
  x_breaks = 2:5,
  x_limits = c(2, 5),
  error_breaks = seq(0, 1.25, by = 0.25),
  error_limits = c(0, 1.25),
  score_breaks = seq(0, 1, by = 0.2),
  score_limits = c(0, 1),
  error_file = "high_dimension_2_error.pdf",
  score_file = "high_dimension_2_score.pdf"
)

# Non-Gaussian design
error_nongaussian <- matrix(
  c(
    0.3721, 1.0648, 1.1101, 0.8116, 1.0074,
    0.2396, 1.0454, 0.9504, 0.7079, 0.5605,
    0.1818, 1.0379, 0.8165, 0.6403, 0.2898,
    0.1463, 1.0338, 0.7138, 0.5986, 0.2447,
    0.1238, 1.0147, 0.6044, 0.5456, 0.1984,
    0.1289, 0.9862, 0.4936, 0.5011, 0.2350
  ),
  byrow = TRUE,
  ncol = length(standard_methods)
)

score_nongaussian <- matrix(
  c(
    0.8371, 0.0290, 0.3375, 0.1818, 0.0105,
    0.8718, 0.0508, 0.4994, 0.1818, 0.5856,
    0.8797, 0.0926, 0.6624, 0.1818, 0.8093,
    0.8766, 0.1318, 0.7303, 0.1818, 0.7809,
    0.8867, 0.1959, 0.7540, 0.1818, 0.7776,
    0.8780, 0.2714, 0.6923, 0.1818, 0.7868
  ),
  byrow = TRUE,
  ncol = length(standard_methods)
)

plot_pair(
  error_matrix = error_nongaussian,
  score_matrix = score_nongaussian,
  x_values = seq(50, 300, by = 50),
  method_ids = standard_methods,
  shapes = standard_shapes,
  linetypes = standard_linetypes,
  colors = standard_colors,
  x_label = "Local sample size",
  x_breaks = seq(50, 300, by = 50),
  x_limits = c(50, 300),
  error_breaks = seq(0, 1.25, by = 0.25),
  error_limits = c(0, 1.25),
  score_breaks = seq(0, 1, by = 0.2),
  score_limits = c(0, 1),
  error_file = "nongaussian_error.pdf",
  score_file = "nongaussian_score.pdf"
)

# Worst-case comparison (six methods for error, five for F1)
worst_methods_error <- c("A", "B", "C", "D", "E", "F")
worst_shapes_error <- c(A = 15, B = 15, C = 16, D = 16, E = 17, F = 17)
worst_linetypes_error <- c(A = "solid", B = "dotted", C = "solid", D = "dotted", E = "solid", F = "dotted")
worst_colors_error <- c(
  A = "#1f77b4",
  B = "#1f77b4",
  C = "#ff7f0e",
  D = "#ff7f0e",
  E = "#2ca02c",
  F = "#2ca02c"
)

error_worst <- matrix(
  c(
    1.0533, 1.0732, 0.4295, 0.4543, 0.2829, 0.2931,
    1.0201, 1.0402, 0.2813, 0.2912, 0.1443, 0.1532,
    0.9246, 0.9323, 0.2108, 0.2203, 0.0993, 0.1012,
    0.6964, 0.7034, 0.1655, 0.1783, 0.0744, 0.0854,
    0.4869, 0.4932, 0.1398, 0.1442, 0.0635, 0.0722,
    0.3868, 0.3932, 0.1191, 0.1222, 0.0553, 0.0633
  ),
  byrow = TRUE,
  ncol = length(worst_methods_error)
)

error_worst_df <- build_long_form(error_worst, seq(50, 300, by = 50), worst_methods_error)

plot_metric(
  data = error_worst_df,
  x_label = "Local sample size",
  y_label = TeX("$L_2$-error"),
  x_breaks = seq(50, 300, by = 50),
  x_limits = c(50, 300),
  y_breaks = seq(0, 1.25, by = 0.25),
  y_limits = c(0, 1.25),
  shapes = worst_shapes_error,
  linetypes = worst_linetypes_error,
  colors = worst_colors_error,
  file_name = "worst_case_error.pdf"
)

worst_methods_score <- c("A", "B", "C", "D", "E")
worst_shapes_score <- c(A = 15, B = 15, C = 16, D = 17, E = 17)
worst_linetypes_score <- c(A = "solid", B = "dotted", C = "solid", D = "solid", E = "dotted")
worst_colors_score <- c(
  A = "#1f77b4",
  B = "#1f77b4",
  C = "#ff7f0e",
  D = "#2ca02c",
  E = "#2ca02c"
)

score_worst <- matrix(
  c(
    0.0628, 0.0523, 0.1818, 0.8433, 0.8231,
    0.1636, 0.1533, 0.1818, 0.8290, 0.8122,
    0.3846, 0.3667, 0.1818, 0.8267, 0.8099,
    0.7115, 0.7022, 0.1818, 0.8208, 0.8001,
    0.8662, 0.8553, 0.1818, 0.8494, 0.8339,
    0.8957, 0.8823, 0.1818, 0.8652, 0.8553
  ),
  byrow = TRUE,
  ncol = length(worst_methods_score)
)

score_worst_df <- build_long_form(score_worst, seq(50, 300, by = 50), worst_methods_score)

plot_metric(
  data = score_worst_df,
  x_label = "Local sample size",
  y_label = TeX("$F_1$-score"),
  x_breaks = seq(50, 300, by = 50),
  x_limits = c(50, 300),
  y_breaks = seq(0, 1, by = 0.2),
  y_limits = c(0, 1),
  shapes = worst_shapes_score,
  linetypes = worst_linetypes_score,
  colors = worst_colors_score,
  file_name = "worst_case_score.pdf"
)

# Heterogeneity proportion study
error_heter <- matrix(
  c(
    0.0871, 1.0108, 0.5667, 0.2626, 0.1398,
    0.0889, 1.0117, 0.5705, 0.2730, 0.1438,
    0.0895, 1.0206, 0.5820, 0.2788, 0.1530,
    0.0906, 1.0235, 0.5943, 0.2809, 0.1603,
    0.0911, 1.0298, 0.5994, 0.2885, 0.1786
  ),
  byrow = TRUE,
  ncol = length(standard_methods)
)

score_heter <- matrix(
  c(
    0.9531, 0.1882, 0.8432, 0.1818, 0.9210,
    0.9444, 0.1831, 0.8322, 0.1818, 0.9100,
    0.9421, 0.1716, 0.8263, 0.1818, 0.8810,
    0.9338, 0.1656, 0.8203, 0.1818, 0.8699,
    0.9223, 0.1578, 0.8182, 0.1818, 0.8553
  ),
  byrow = TRUE,
  ncol = length(standard_methods)
)

plot_pair(
  error_matrix = error_heter,
  score_matrix = score_heter,
  x_values = seq(0, 0.8, by = 0.2),
  method_ids = standard_methods,
  shapes = standard_shapes,
  linetypes = standard_linetypes,
  colors = standard_colors,
  x_label = "Proportion",
  x_breaks = seq(0, 0.8, by = 0.2),
  x_limits = c(0, 0.8),
  error_breaks = seq(0, 1.25, by = 0.25),
  error_limits = c(0, 1.25),
  score_breaks = seq(0, 1, by = 0.2),
  score_limits = c(0, 1),
  error_file = "heterogeneity_error.pdf",
  score_file = "heterogeneity_score.pdf"
)

