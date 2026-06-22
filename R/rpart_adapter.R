#' Fit an rpart model
#'
#' @description Adapter function to fit a recursive partitioning tree.
#'
#' @param object An S3 method spec object.
#' @param data A data frame containing the training features.
#' @param labels A vector (usually a factor) of target labels.
#' @param ... Additional arguments passed to `rpart::rpart()`.
#'
#' @return A fitted `rpart` object.
#' @export
fit_adapter.classbound_rpart <- function(object, data, labels, ...) {
  if (!requireNamespace("rpart", quietly = TRUE)) {
    stop("Package 'rpart' must be installed to use method = 'rpart'.", call. = FALSE)
  }

  # Combine data and labels into a single formula for rpart
  train_data <- data
  train_data$.label. <- labels

  # Fit the model
  rpart::rpart(.label. ~ ., data = train_data, ...)
}

#' Predict using a fitted rpart model
#'
#' @description Adapter function to generate standardized predictions from an rpart model.
#'
#' @param model A fitted `rpart` object.
#' @param newdata A data frame of new observations to predict on.
#' @param ... Additional arguments passed to `predict()`.
#'
#' @return A list containing `class` (predicted labels) and `probs` (probabilities).
#' @importFrom stats predict
#' @export
predict_model.classbound_rpart <- function(model, newdata, ...) {
  raw_model <- model$model
  if (!inherits(raw_model, "rpart")) {
    stop("Model must be an 'rpart' object.", call. = FALSE)
  }

  # Predict classes (discrete labels)
  preds <- predict(raw_model, newdata, type = "class", ...)

  # Predict probabilities (continuous matrix)
  probs <- predict(raw_model, newdata, type = "prob", ...)

  list(
    class = preds,
    probs = probs
  )
}
