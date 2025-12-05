#' @title Generate synthetic 1-bit data
#' @description Generate synthetic 1-bit data for decentralized 1-bit
#' compressed sensing.
#' @name genData
#' @param N total sample size
#' @param m number of nodes
#' @param p number of features
#' @param s sparsity level of the true signal
#' @param pc connection probability for the network graph
#' @param p_flip flipping probability for the 1-bit measurements
#' @param rho correlation parameter for the design matrix
#' @param sigma2 variance parameter for the design matrix
#' @param type error distribution type, can be "Cauchy", "Normal" or "T2"
#' @param sigma2_set set of variance parameters for heterogeneous data
#' @param rho_set set of correlation parameters for heterogeneous data
#' @param p_flip_set set of flipping probabilities for heterogeneous data
#' @param case data generation case, can be 0, 1, 2 or 3
#' @return A list consists of
#' \item{X}{design matrix}
#' \item{y}{response vector}
#' \item{betaT}{true signal vector}
#' \item{graph}{network graph}
#' @export
genData <- function(N = 2000, m = 10, p = 100, s = 5,
                    pc = 0.3, p_flip = 0.2, rho = 0.1, sigma2 = 1,
                    type = c("Cauchy", "Normal", "T2"),
                    sigma2_set = c(1, 3),
                    rho_set = c(0.1, 0.3),
                    p_flip_set = c(0.05, 0.1),
                    case = 1) {
  type <- match.arg(type)
  if (!case %in% 0:3) {
    stop("`case` must be one of 0, 1, 2 or 3.")
  }
  if (pc < 0 || pc > 1) {
    stop("`pc` must belong to [0, 1].")
  }
  if (p_flip < 0 || p_flip > 1) {
    stop("`p_flip` must belong to [0, 1].")
  }
  if (any(p_flip_set < 0 | p_flip_set > 1)) {
    stop("All `p_flip_set` elements must belong to [0, 1].")
  }
  if (N %% m != 0) {
    stop("`N` must be divisible by `m` so that each node has the same sample size.")
  }
  if (s > p) {
    stop("`s` cannot exceed the ambient dimension `p`.")
  }

  n <- as.integer(N / m)
  sigma_builder <- function(var, corr) {
    var * toeplitz(corr^seq(0, p - 1, by = 1))
  }
  draw_noise <- function(dist, n_obs) {
    switch(dist,
           Cauchy = rcauchy(n_obs),
           Normal = rnorm(n_obs),
           T2 = rt(n_obs, df = 2),
           stop("Unsupported noise type: ", dist))
  }
  generate_beta <- function() {
    beta <- matrix(0, p, 1)
    idx <- sample.int(p, s)
    beta[idx] <- 2 * rbinom(s, 1, 0.5) - 1
    beta
  }
  apply_global_flips <- function(y_vec, prob) {
    if (prob < 0 || prob > 1) {
      stop("Flip probabilities must belong to [0, 1].")
    }
    n_flip <- floor(length(y_vec) * prob)
    if (n_flip > 0) {
      flip_id <- sample.int(length(y_vec), n_flip)
      y_vec[flip_id] <- -y_vec[flip_id]
    }
    y_vec
  }

  case_result <- switch(case + 1,
    {
      betaT <- generate_beta()
      Sigma <- sigma_builder(sigma2, rho)
      X <- MASS::mvrnorm(N, rep(0, p), Sigma)
      noise <- draw_noise(type, N)
      y <- sign(drop(X %*% betaT) + noise)
      y <- apply_global_flips(y, p_flip)
      list(X = X, y = y, betaT = betaT)
    },
    {
      betaT <- generate_beta()
      X <- matrix(NA_real_, nrow = N, ncol = p)
      noise <- draw_noise(type, N)
      y <- numeric(N)
      for (j in seq_len(m)) {
        idx <- ((j - 1) * n + 1):(j * n)
        Sigma <- sigma_builder(sample(sigma2_set, 1), sample(rho_set, 1))
        X[idx, ] <- MASS::mvrnorm(n, rep(0, p), Sigma)
        y[idx] <- sign(drop(X[idx, , drop = FALSE] %*% betaT) + noise[idx])
      }
      y <- apply_global_flips(y, p_flip)
      list(X = X, y = y, betaT = betaT)
    },
    { # case == 2
      betaT <- generate_beta()
      Sigma <- sigma_builder(sigma2, rho)
      X <- MASS::mvrnorm(N, rep(0, p), Sigma)
      block_ids <- split(seq_len(N), rep(seq_len(m), each = n))
      dist_map <- c("Cauchy", "Normal", "T2")
      node_types <- sample.int(length(dist_map), m, replace = TRUE)
      noise <- numeric(N)
      for (j in seq_len(m)) {
        block <- block_ids[[j]]
        noise[block] <- draw_noise(dist_map[node_types[j]], length(block))
      }
      y <- sign(drop(X %*% betaT) + noise)
      y <- apply_global_flips(y, p_flip)
      list(X = X, y = y, betaT = betaT)
    },
    {
      betaT <- generate_beta()
      Sigma <- sigma_builder(sigma2, rho)
      X <- MASS::mvrnorm(N, rep(0, p), Sigma)
      noise <- draw_noise(type, N)
      y <- sign(drop(X %*% betaT) + noise)
      block_ids <- split(seq_len(N), rep(seq_len(m), each = n))
      for (block in block_ids) {
        p_flip_node <- sample(p_flip_set, 1)
        n_flip <- floor(length(block) * p_flip_node)
        if (n_flip > 0) {
          flip_id <- sample(block, n_flip)
          y[flip_id] <- -y[flip_id]
        }
      }
      list(X = X, y = y, betaT = betaT)
    }
  )

  graph <- igraph::sample_gnp(m, pc)
  if (m > 1) {
    tries <- 0L
    while (!igraph::is_connected(graph)) {
      graph <- igraph::sample_gnp(m, pc)
      tries <- tries + 1L
      if (tries > 1000L) {
        stop("Failed to sample a connected graph after 1000 attempts. Consider increasing `pc`.")
      }
    }
  }

  c(case_result, list(graph = graph))
}

#' @title Create directory if not exist
#' @description Create directory if not exist
#' @name createdir
#' @param path directory path
#' @param recursive logical, whether to create directories recursively
#' @return TRUE if directory is created, FALSE otherwise
#' @export
createdir <- function(path, recursive = TRUE) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = recursive)
  } else {
    FALSE
  }
}

#' @title Combine multiple lists into one list
#' @description Combine multiple lists into one list by concatenating the
#' elements with the same index.
#' @name comb
#' @export
comb <- function(x, ...) {
  if (!is.list(x)) {
    stop("`x` must be a list.")
  }
  others <- list(...)
  if (length(others) == 0) {
    return(x)
  }
  inputs <- c(list(x), others)
  lengths <- vapply(inputs, length, integer(1), USE.NAMES = FALSE)
  if (!all(lengths == lengths[1])) {
    stop("All lists must have the same length.")
  }
  lapply(seq_len(lengths[1]), function(i) {
    elements <- lapply(inputs, `[[`, i)
    do.call(c, elements)
  })
}

#' @title Compute F1 score of two set
#' @description Compute \eqn{F_1} score of two set. Note that \eqn{F_1} score is
#' symmetric.
#' @name computeF1
#' @param esupp the estimated support
#' @param supp the true support
#' @return A list consists of
#' \item{f1}{the computed F1 score}
#' \item{precision}{the precision}
#' \item{recall}{the recall}
#' \item{fp}{false positve}
#' \item{fn}{false negative}
#' @references
#' \url{https://en.wikipedia.org/wiki/F-score}
#' @export
computeF1 <- function(esupp, supp) {
  esupp <- unique(esupp[!is.na(esupp)])
  supp <- unique(supp[!is.na(supp)])
  tp <- length(intersect(esupp, supp))
  fp <- length(setdiff(esupp, supp))
  fn <- length(setdiff(supp, esupp))
  precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
  f1 <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  list(f1 = f1,
       precision = precision,
       recall = recall,
       fp = fp,
       fn = fn)
}


#' @title Compute error between two vectors or two matrices
#' @description Compute the error between two vectors or two matrices using
#' specified norm.
#' @name computeError
#' @param x first vector or matrix
#' @param y second vector or matrix
#' @param type norm type, default is "F" (Frobenius norm)
#' @return The computed error
#' @export
computeError <- function(x, y, type = "F") {
  stopifnot(is.numeric(x) && is.numeric(y))
  if ((length(x) != length(y)) && is.matrix(x)) {
    m <- ncol(x)
    for (j in 1:m) {
      if (sum(abs(x[, j])) > 1e-12) {
        x[, j] <- x[, j] / norm(x[, j, drop = F], "F")
      }
    }
    sqrt(norm(x - matrix(rep(y, ncol(x)), ncol = ncol(x)), type = type)^2 / m)
  } else {
    stopifnot(length(x) == length(y))
    x <- x / norm(matrix(x), "F")
    return(drop(norm(matrix(x - y), type = type)))
  }
}


#' @title BIC for glmnet
#' @description Compute BIC for glmnet and return the estimated coefficient
#' @name bic.glmnet
#' @param X design matrix
#' @param y response vector
#' @param nlambda number of lambda values
#' @return The estimated coefficient vector selected by BIC
#' @export
bic.glmnet <- function(X, y, nlambda = 100) {
  X <- as.matrix(X)
  stopifnot(is.numeric(y))
  N <- nrow(X)
  if (N != length(y)) {
    stop("`y` must have length equal to the number of rows of `X`.")
  }
  fit <- glmnet::glmnet(
    x = X,
    y = y,
    intercept = FALSE,
    nlambda = nlambda
  )
  beta_mat <- as.matrix(fit$beta)
  if (ncol(beta_mat) == 0) {
    stop("glmnet returned no coefficient estimates.")
  }
  fitted <- X %*% beta_mat
  residuals <- matrix(y, nrow = N, ncol = ncol(beta_mat)) - fitted
  rss <- colSums(residuals^2)
  df <- colSums(beta_mat != 0)
  bic_values <- rss / N + log(N) / N * df
  best <- which.min(bic_values)
  matrix(beta_mat[, best], ncol = 1)
}

#' @title Fix RNG stream for parallel computing
#' @description Given an initial RNG stream, fix the RNG stream after certain
#' number of streams.
#' @name fixRNGStream
#' @param s initial RNG stream
#' @param batch number of streams to advance
#' @return The fixed RNG stream after certain number of streams.
#'
#' @export
fixRNGStream <- function(s, batch) {
  if (batch < 0) {
    stop("`batch` must be non-negative.")
  }
  if (batch == 0) {
    return(s)
  }
  for (i in seq_len(batch)) {
    s <- parallel::nextRNGStream(s)
  }
  s
}

#' @title Calculate indices for node j
#' @description Calculate the indices of samples assigned to node j
#' @name calN_j
#' @param n number of samples per node
#' @param j node index
#' @return A vector of indices assigned to node j
#' @export
calN_j <- function(n, j) {
  return((n * (j - 1) + 1):(n * j))
}