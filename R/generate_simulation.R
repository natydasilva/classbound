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

  if (is.null(class_names)) {
    if (num_classes <= 26) {
      class_names <- paste0("Class ", LETTERS[1:num_classes])
    } else {
      class_names <- paste0("Class ", 1:num_classes)
    }
  }

  sim_data <- lapply(1:num_classes, function(i) {
    sim <- MASS::mvrnorm(n = ns[i], mu = means[[i]], Sigma = covs[[i]])
    df <- data.frame(Sim = class_names[i], sim)
    colnames(df)[-1] <- paste0("X", seq_len(ncol(sim)))
    df
  })

  do.call(rbind, sim_data)
}
