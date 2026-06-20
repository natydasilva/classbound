#' Fit a randomForest model
#'
#' @description Adapter function to fit a random forest model.
#'
#' @param object An S3 method spec object.
#' @param data A data frame containing the training features.
#' @param labels A vector (usually a factor) of target labels.
#' @param ... Additional arguments passed to `randomForest::randomForest()`.
#'
#' @return A fitted `randomForest` object.
#' @export
fit_adapter.classbound_randomForest <- function(object, data, labels, ...) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' must be installed to use method = 'randomForest'.", call. = FALSE)
  }

  train_data <- data
  train_data$.label. <- as.factor(labels)

  randomForest::randomForest(.label. ~ ., data = train_data, ...)
}

#' Predict using a fitted randomForest model
#'
#' @description Adapter function to generate standardized predictions from a randomForest model.
#'
#' @param model A fitted `randomForest` object.
#' @param newdata A data frame of new observations to predict on.
#' @param ... Additional arguments passed to `predict()`.
#'
#' @return A list containing `class` (predicted labels) and `probs` (probabilities).
#' @importFrom stats predict
#' @export
predict_model.classbound_randomForest <- function(model, newdata, ...) {
  raw_model <- model$model
  if (!inherits(raw_model, "randomForest")) {
    stop("Model must be a 'randomForest' object.", call. = FALSE)
  }

  preds <- predict(raw_model, newdata, type = "response", ...)
  probs <- predict(raw_model, newdata, type = "prob", ...)

  list(
    class = preds,
    probs = probs
  )
}
