#' Compute the classification decision boundary
#'
#' @description Computes the decision boundary for a fitted classifier over a specified feature space.
#'
#' @param model A fitted model object of class `classbound_model`.
#' @param range A named list defining the feature ranges for a 2D space. Each element should be a numeric vector of length 2 (e.g., `list(Sepal.Length = c(4, 8), Sepal.Width = c(2, 5))`).
#' @param resolution An integer specifying the grid resolution.
#' @param ... Additional computation parameters.
#'
#' @return A data frame containing the grid points (`x`, `y`), predicted `prediction`, and probabilities.
#' @export
boundary_compute <- function(model, range, resolution = 100, ...) {
  if (!inherits(model, "classbound_model")) {
    stop("model must be a 'classbound_model' object.", call. = FALSE)
  }
  
  if (!is.list(range) || length(range) != 2 || is.null(names(range))) {
    stop("range must be a named list of length 2 specifying feature ranges.", call. = FALSE)
  }
  
  var_names <- names(range)
  
  # Generate grid points for each dimension
  seq_x <- seq(from = min(range[[1]]), to = max(range[[1]]), length.out = resolution)
  seq_y <- seq(from = min(range[[2]]), to = max(range[[2]]), length.out = resolution)
  
  # Create grid
  grid_df <- expand.grid(seq_x, seq_y)
  colnames(grid_df) <- var_names
  
  # Predict using the standardized predict_model function
  preds <- predict_model(model, newdata = grid_df, ...)
  
  # Combine grid and predictions
  res <- data.frame(
    x = grid_df[[1]],
    y = grid_df[[2]],
    prediction = preds$class
  )
  
  if (!is.null(preds$probs)) {
    res <- cbind(res, as.data.frame(preds$probs))
  }
  
  res
}
