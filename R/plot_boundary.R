#' Visualize the classification decision boundary
#'
#' @description Renders a 2D plot for boundary exploration.
#'
#' @param boundary The boundary data frame returned by `boundary_compute()`.
#'   Must contain columns named `x`, `y`, and `prediction`.
#' @param obs_data Optional data frame of training observations to overlay.
#' @param x_col Column name in `obs_data` for the x-axis feature. Also used as the x-axis label.
#'   Required if `obs_data` is provided.
#' @param y_col Column name in `obs_data` for the y-axis feature. Also used as the y-axis label.
#'   Required if `obs_data` is provided.
#' @param true_label Column name in `obs_data` representing the true class labels. Required if `obs_data` is provided.
#' @param type The type of visualization to generate. Only '2D' is supported.
#' @param ... Additional visualization parameters.
#'
#' @return A `ggplot2` object.
#' @examples
#' \donttest{
#' data(data69_1)
#' data69_1$Y <- as.factor(data69_1$Y)
#' model <- fit_model(data69_1, Y ~ V1 + V2, rpart::rpart)
#' grid <- boundary_compute(model, list(V1 = c(-1, 1), V2 = c(-1, 1)))
#' plot_boundary(grid, data69_1, "V1", "V2", "Y")
#' }
#' @export
plot_boundary <- function(boundary, obs_data = NULL, x_col = NULL, y_col = NULL, true_label = NULL, type = "2D", ...) {
  if (type != "2D") {
    stop("Only type='2D' is currently supported.", call. = FALSE)
  }

  if (!inherits(boundary, "data.frame")) {
    stop("boundary must be a data.frame returned by boundary_compute().", call. = FALSE)
  }

  req_cols <- c("x", "y", "prediction")
  if (!all(req_cols %in% colnames(boundary))) {
    stop("boundary must contain 'x', 'y', and 'prediction' columns.", call. = FALSE)
  }

  # Basic plot
  p <- ggplot2::ggplot(boundary, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = prediction), alpha = 0.3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = if (!is.null(x_col)) x_col else "Feature 1",
      y = if (!is.null(y_col)) y_col else "Feature 2",
      fill = "Prediction"
    )

  # Overlay training observations if provided
  if (!is.null(obs_data)) {
    if (is.null(x_col) || is.null(y_col) || is.null(true_label)) {
      stop("When providing obs_data, you must also specify x_col, y_col, and true_label.", call. = FALSE)
    }

    if (!inherits(obs_data, "data.frame")) {
      stop("obs_data must be a data.frame.", call. = FALSE)
    }

    if (!all(c(x_col, y_col, true_label) %in% colnames(obs_data))) {
      stop(sprintf("obs_data must contain columns '%s', '%s', and '%s'.", x_col, y_col, true_label), call. = FALSE)
    }

    # Map the dynamic columns to standard names for ggplot evaluation
    obs_df <- obs_data[, c(x_col, y_col, true_label)]
    colnames(obs_df) <- c("x_val", "y_val", "true_class")

    p <- p + ggplot2::geom_point(
      data = obs_df,
      ggplot2::aes(x = x_val, y = y_val, color = true_class),
      size = 2, shape = 16
    ) +
      ggplot2::labs(color = "True Class")
  }

  p
}
