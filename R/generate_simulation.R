#' Generate simulated datasets for classifier exploration
#'
#' @description Generates reproducible synthetic data for testing and visualizing classification boundaries.
#'
#' @param n An integer, the number of observations to generate.
#' @param type A character string indicating the simulation type (e.g., 'linear', 'spiral', 'clusters').
#' @param ... Additional arguments for simulation parameters.
#'
#' @return A data frame containing simulated features and a 'class' label column.
#' @keywords internal
generate_simulation <- function(n, type, ...) {
  # TODO: Implement data simulation functions
  stop("generate_simulation() is not yet implemented.")
}

#' Simulate data using Gaussian mixture models from MixSim
#'
#' @param n Sample size.
#' @param K Number of classes.
#' @param p Number of dimensions.
#' @param MaxOmega Maximum overlap between components.
#' @param class_names Optional character vector of class labels.
#' @param seed Optional integer for reproducibility.
#' @param noise_ratio Numeric between 0 and 1. Proportion of `n` to add as uniform background noise.
#' @param test_ratio Numeric between 0 and 1. If greater than 0, independently generates an additional test dataset using the exact same parameters.
#'
#' @return If `test_ratio == 0`, a data frame containing simulated features and a 'Sim' label column.
#'   If `test_ratio > 0`, a list with `$train` and `$test` data frames.
#' @export
simulate_mixsim <- function(n, K, p, MaxOmega, class_names = NULL, seed = NULL, noise_ratio = 0, test_ratio = 0) {
  if (!requireNamespace("MixSim", quietly = TRUE)) {
    stop("Package 'MixSim' must be installed to use simulate_mixsim().", call. = FALSE)
  }

  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- globalenv()$.Random.seed
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }

  if (is.null(class_names)) {
    class_names <- paste0("Class ", 1:K)
  }

  sim_params <- MixSim::MixSim(MaxOmega = MaxOmega, K = K, p = p)
  
  # Helper to generate and format a single dataset
  generate_df <- function(num_obs) {
    raw <- MixSim::simdataset(n = num_obs, Pi = sim_params$Pi, Mu = sim_params$Mu, S = sim_params$S)
    df <- data.frame(raw$X)
    colnames(df) <- paste0("X", seq_len(p))
    df$Sim <- factor(class_names[raw$id], levels = class_names)
    df <- df[, c("Sim", colnames(df)[-ncol(df)])]
    
    if (noise_ratio > 0) {
      n_noise <- ceiling(num_obs * noise_ratio)
      noise_df <- data.frame(matrix(
        stats::runif(n_noise * p, min = apply(df[, -1, drop = FALSE], 2, min), max = apply(df[, -1, drop = FALSE], 2, max)), 
        ncol = p, byrow = TRUE
      ))
      colnames(noise_df) <- colnames(df)[-1]
      noise_df$Sim <- sample(levels(df$Sim), n_noise, replace = TRUE)
      df <- rbind(df, noise_df[, colnames(df)])
    }
    return(df)
  }

  train_df <- generate_df(n)
  
  if (test_ratio > 0) {
    test_n <- max(1, round(n * test_ratio))
    test_df <- generate_df(test_n)
    return(list(train = train_df, test = test_df))
  }

  return(train_df)
}

#' Simulate multivariate normal distribution data for N classes
#'
#' @param means A list of numeric vectors containing the means for each class.
#' @param covs A list of covariance matrices for each class.
#' @param ns A numeric vector of sample sizes for each class.
#' @param class_names A character vector of class labels.
#' @param seed Optional integer for reproducibility.
#' @param noise_ratio Numeric between 0 and 1. Proportion of `sum(ns)` to add as uniform background noise.
#' @param test_ratio Numeric between 0 and 1. If greater than 0, independently generates an additional test dataset using the exact same parameters.
#'
#' @return If `test_ratio == 0`, a data frame containing simulated features and a 'Sim' label column.
#'   If `test_ratio > 0`, a list with `$train` and `$test` data frames.
#' @export
simu_n <- function(means, covs, ns, class_names = NULL, seed = NULL, noise_ratio = 0, test_ratio = 0) {
  num_classes <- length(ns)

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' must be installed to use simu_n().", call. = FALSE)
  }
  
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- globalenv()$.Random.seed
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }

  if (is.null(class_names)) {
    class_names <- paste0("Class ", 1:num_classes)
  }

  generate_df <- function(sizes) {
    sim_data <- lapply(1:num_classes, function(i) {
      if (sizes[i] == 0) return(NULL)
      sim <- MASS::mvrnorm(n = sizes[i], mu = means[[i]], Sigma = covs[[i]])
      if (sizes[i] == 1) sim <- matrix(sim, nrow = 1)
      df <- data.frame(Sim = class_names[i], sim)
      colnames(df)[-1] <- paste0("X", seq_len(ncol(sim)))
      df
    })
  
    df <- do.call(rbind, sim_data)
    if (is.null(df)) return(NULL)
    
    df$Sim <- factor(df$Sim, levels = class_names)
    
    if (noise_ratio > 0) {
      n_noise <- ceiling(sum(sizes) * noise_ratio)
      p <- ncol(df) - 1
      noise_df <- data.frame(matrix(
        stats::runif(n_noise * p, min = apply(df[, -1, drop = FALSE], 2, min), max = apply(df[, -1, drop = FALSE], 2, max)), 
        ncol = p, byrow = TRUE
      ))
      colnames(noise_df) <- colnames(df)[-1]
      noise_df$Sim <- sample(levels(df$Sim), n_noise, replace = TRUE)
      df <- rbind(df, noise_df[, colnames(df)])
    }
    
    return(df)
  }

  train_df <- generate_df(ns)
  
  if (test_ratio > 0) {
    test_ns <- pmax(1, round(ns * test_ratio))
    test_df <- generate_df(test_ns)
    return(list(train = train_df, test = test_df))
  }
  
  return(train_df)
}
