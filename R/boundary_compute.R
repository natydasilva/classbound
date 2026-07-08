#' Compute the classification decision boundary
#'
#' @description Computes the decision boundary for a fitted classifier over a specified feature space.
#'
#' @param model A fitted model object of class `classbound`.
#' @param range A named list defining the feature ranges for a 2D space. Each element should be a numeric vector of length 2 (e.g., `list(Sepal.Length = c(4, 8), Sepal.Width = c(2, 5))`).
#' @param resolution An integer specifying the grid resolution.
#' @param predfun A custom function to generate predictions for non-standard models.
#' @param ... Additional computation parameters.
#'
#' @return A data frame containing the grid points (`x`, `y`), predicted `prediction`, and probabilities.
#' @examples
#' \donttest{
#' data(data69_1)
#' data69_1$Y <- as.factor(data69_1$Y)
#' model <- fit_model(data69_1, Y ~ V1 + V2, rpart::rpart)
#' grid <- boundary_compute(model, list(V1 = c(-1, 1), V2 = c(-1, 1)))
#' }
#' @export
boundary_compute <- function(model, range, resolution = 100, predfun = NULL, ...) {
  if (!inherits(model, "classbound")) {
    stop("model must be a 'classbound' object.", call. = FALSE)
  }

  if (!is.list(range) || length(range) != 2 || is.null(names(range))) {
    stop("range must be a named list of length 2 specifying feature ranges.", call. = FALSE)
  }

  var_names <- names(range)

  if (any(duplicated(var_names))) {
    stop("Duplicate feature names found in `range`.", call. = FALSE)
  }

  if (is.null(model$features)) {
    stop("Model metadata is missing. Please retrain the model.", call. = FALSE)
  }

  expected_names <- model$features$names
  
  if (length(var_names) != length(expected_names)) {
    stop(sprintf("Expected %d boundary variables, but got %d.", length(expected_names), length(var_names)), call. = FALSE)
  }
  
  missing_names <- setdiff(expected_names, var_names)
  invalid_names <- setdiff(var_names, expected_names)
  
  if (length(missing_names) > 0 || length(invalid_names) > 0) {
    msg <- "Names in `range` do not match the training features."
    if (length(missing_names) > 0) {
      msg <- paste0(msg, "\nMissing features: ", paste(missing_names, collapse = ", "))
    }
    if (length(invalid_names) > 0) {
      msg <- paste0(msg, "\nInvalid features: ", paste(invalid_names, collapse = ", "))
    }
    stop(msg, call. = FALSE)
  }

  # Ensure requested boundary features are numeric
  feature_types <- model$features$types[var_names]
  is_numeric <- vapply(feature_types, function(x) any(c("numeric", "integer") %in% x), logical(1))
  
  if (any(!is_numeric)) {
    stop("Boundary generation requires numeric features. Categorical ranges are not supported.", call. = FALSE)
  }
  # Generate grid points for each dimension
  seq_x <- seq(from = min(range[[1]]), to = max(range[[1]]), length.out = resolution)
  seq_y <- seq(from = min(range[[2]]), to = max(range[[2]]), length.out = resolution)

  # Create grid
  grid_df <- expand.grid(seq_x, seq_y)
  colnames(grid_df) <- var_names

  # Predict using the standardized predict_model function
  preds <- predict_model(model, newdata = grid_df, predfun = predfun, ...)

  # Class predictions are now strictly guaranteed to be correctly leveled factors 
  # by the predict_model() contract, so we directly construct the frame.

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
