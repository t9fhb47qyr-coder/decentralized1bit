library(ggplot2)
library(latex2exp)

method_labels <- c("A", "B", "C", "D", "E")
method_shapes <- c(A = 15, B = 16, C = 17, D = 18, E = 6)
method_linetypes <- c(A = "solid", B = "solid", C = "dotted", D = "dashed", E = "dashed")
method_colors <- c(
  A = "#1f77b4",
  B = "#ff7f0e",
  C = "#2ca02c",
  D = "#d62728",
  E = "#9467bd"
)

output_dir <- file.path("Output", "plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

build_long_form <- function(matrix_values, x_values, metric_name) {
  stopifnot(ncol(matrix_values) == length(method_labels))
  metric_values <- as.vector(matrix_values)
  data.frame(
    x = rep(x_values, times = length(method_labels)),
    method = rep(method_labels, each = length(x_values)),
    value = metric_values,
    metric = metric_name
  )
}

plot_metric <- function(data, x_label, y_label, y_breaks, y_limits, file_name) {
  plot_obj <- ggplot(
    data = data,
    mapping = aes(x = x, y = value, colour = method, shape = method, linetype = method)
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    scale_shape_manual(values = method_shapes) +
    scale_linetype_manual(values = method_linetypes) +
    scale_colour_manual(values = method_colors) +
    scale_x_continuous(name = x_label, breaks = unique(data$x)) +
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
                      x_label,
                      error_breaks,
                      error_limits,
                      score_breaks,
                      score_limits,
                      error_file,
                      score_file) {
  error_df <- build_long_form(error_matrix, x_values, "error")
  score_df <- build_long_form(score_matrix, x_values, "score")

  plot_metric(
    data = error_df,
    x_label = x_label,
    y_label = TeX("$L_2$-error"),
    y_breaks = error_breaks,
    y_limits = error_limits,
    file_name = error_file
  )

  plot_metric(
    data = score_df,
    x_label = x_label,
    y_label = TeX("$F_1$-score"),
    y_breaks = score_breaks,
    y_limits = score_limits,
    file_name = score_file
  )
}

# Figure set 1: local sample size
error_local_sample <- matrix(
  c(
    0.1390, 1.0533, 0.8716, 0.4295, 0.2829,
    0.0909, 1.0201, 0.5644, 0.2813, 0.1443,
    0.0692, 0.9246, 0.2965, 0.2108, 0.0993,
    0.0572, 0.6964, 0.1586, 0.1655, 0.0744,
    0.0511, 0.4869, 0.1143, 0.1398, 0.0635,
    0.0465, 0.3868, 0.0916, 0.1191, 0.0553
  ),
  byrow = TRUE,
  ncol = length(method_labels)
)

score_local_sample <- matrix(
  c(
    0.9362, 0.0628, 0.5946, 0.1818, 0.8433,
    0.9512, 0.1636, 0.8310, 0.1818, 0.8290,
    0.9578, 0.3846, 0.7350, 0.1818, 0.8267,
    0.9612, 0.7115, 0.4783, 0.1818, 0.8208,
    0.9642, 0.8662, 0.3867, 0.1818, 0.8494,
    0.9628, 0.8957, 0.3703, 0.1818, 0.8652
  ),
  byrow = TRUE,
  ncol = length(method_labels)
)

plot_pair(
  error_matrix = error_local_sample,
  score_matrix = score_local_sample,
  x_values = seq(50, 300, by = 50),
  x_label = "Local sample size",
  error_breaks = seq(0, 1.25, by = 0.25),
  error_limits = c(0, 1.25),
  score_breaks = seq(0, 1, by = 0.2),
  score_limits = c(0, 1),
  error_file = "sample_error.pdf",
  score_file = "sample_score.pdf"
)

# Figure set 2: number of nodes
error_nodes <- matrix(
  c(
    0.2321, 1.0178, 0.8793, 0.5675, 0.2297,
    0.1390, 1.0188, 0.7168, 0.4150, 0.1823,
    0.1045, 1.0203, 0.6392, 0.3161, 0.1682,
    0.0909, 1.0201, 0.5664, 0.2813, 0.1443,
    0.0761, 1.0188, 0.5250, 0.2583, 0.1187,
    0.0692, 1.0214, 0.4879, 0.2585, 0.1096,
    0.0628, 1.0198, 0.4580, 0.2782, 0.1003,
    0.0572, 1.0212, 0.4242, 0.2886, 0.0944,
    0.0562, 1.0208, 0.4069, 0.2839, 0.0912,
    0.0511, 1.0216, 0.3950, 0.2739, 0.0879
  ),
  byrow = TRUE,
  ncol = length(method_labels)
)

score_nodes <- matrix(
  c(
    0.9188, 0.1574, 0.5867, 0.1818, 0.3839,
    0.9362, 0.1647, 0.7415, 0.1818, 0.7848,
    0.9524, 0.1576, 0.8211, 0.1818, 0.8595,
    0.9512, 0.1636, 0.8310, 0.1818, 0.8290,
    0.9556, 0.1529, 0.8467, 0.1818, 0.7906,
    0.9578, 0.1571, 0.8385, 0.1818, 0.7785,
    0.9627, 0.1600, 0.8163, 0.1818, 0.7621,
    0.9612, 0.1590, 0.7996, 0.1818, 0.7527,
    0.9591, 0.1566, 0.8050, 0.1818, 0.7347,
    0.9642, 0.1613, 0.7628, 0.1818, 0.7318
  ),
  byrow = TRUE,
  ncol = length(method_labels)
)

plot_pair(
  error_matrix = error_nodes,
  score_matrix = score_nodes,
  x_values = seq(5, 50, by = 5),
  x_label = "Number of nodes",
  error_breaks = seq(0, 1.25, by = 0.25),
  error_limits = c(0, 1.25),
  score_breaks = seq(0, 1, by = 0.2),
  score_limits = c(0, 1),
  error_file = "nodes_error.pdf",
  score_file = "nodes_score.pdf"
)

# Figure set 3: sparsity levels
error_sparsity <- matrix(
  c(
    0.0909, 1.0201, 0.5664, 0.2813, 0.1443,
    0.0909, 1.0201, 0.5644, 0.2745, 0.1463,
    0.0909, 1.0201, 0.5644, 0.2738, 0.1521
  ),
  byrow = TRUE,
  ncol = length(method_labels)
)

score_sparsity <- matrix(
  c(
    0.9512, 0.1636, 0.8310, 0.1818, 0.8290,
    0.9536, 0.1705, 0.8442, 0.1818, 0.8138,
    0.9515, 0.1687, 0.8520, 0.1818, 0.7992
  ),
  byrow = TRUE,
  ncol = length(method_labels)
)

plot_pair(
  error_matrix = error_sparsity,
  score_matrix = score_sparsity,
  x_values = c(5, 10, 15),
  x_label = "Sparsity level",
  error_breaks = seq(0, 1, by = 0.2),
  error_limits = c(0, 1),
  score_breaks = seq(0, 1, by = 0.2),
  score_limits = c(0, 1),
  error_file = "sparsity_error.pdf",
  score_file = "sparsity_score.pdf"
)
