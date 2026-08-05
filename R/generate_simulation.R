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
#'
#' @return A data frame containing simulated features and a 'Sim' label column.
#' @export
simulate_mixsim <- function(n, K, p, MaxOmega, class_names = NULL) {
  if (!requireNamespace("MixSim", quietly = TRUE)) {
    stop("Package 'MixSim' must be installed to use simulate_mixsim().", call. = FALSE)
  }

  if (is.null(class_names)) {
    class_names <- paste0("Class ", 1:K)
  }

  sim_params <- MixSim::MixSim(MaxOmega = MaxOmega, K = K, p = p)
  sim_data_raw <- MixSim::simdataset(n = n, Pi = sim_params$Pi, Mu = sim_params$Mu, S = sim_params$S)

  # Format data frame
  df <- data.frame(sim_data_raw$X)
  colnames(df) <- paste0("X", seq_len(p))

  # Convert class IDs to class names
  class_ids <- sim_data_raw$id
  df$Sim <- factor(class_names[class_ids], levels = class_names)

  # Reorder columns to put 'Sim' first
  df <- df[, c("Sim", colnames(df)[-ncol(df)])]

  return(df)
}

#' Simulate multivariate normal distribution data for N classes
#'
#' @param means A list of numeric vectors containing the means for each class.
#' @param covs A list of covariance matrices for each class.
#' @param ns A numeric vector of sample sizes for each class.
#' @param class_names A character vector of class labels.
#'
#' @return A data frame containing simulated features and a 'Sim' label column.
#' @export
simu_n <- function(means, covs, ns, class_names = NULL) {
  num_classes <- length(ns)

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' must be installed to use simu_n().", call. = FALSE)
  }

  if (is.null(class_names)) {
    class_names <- paste0("Class ", 1:num_classes)
  }

  sim_data <- lapply(1:num_classes, function(i) {
    sim <- MASS::mvrnorm(n = ns[i], mu = means[[i]], Sigma = covs[[i]])
    df <- data.frame(Sim = class_names[i], sim)
    colnames(df)[-1] <- paste0("X", seq_len(ncol(sim)))
    df
  })

  do.call(rbind, sim_data)
}
