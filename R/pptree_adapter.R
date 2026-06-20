#' Fit a PPtreeExt model
#'
#' @description Adapter function to fit a Projection Pursuit classification tree.
#'
#' @param object An S3 method spec object.
#' @param data A data frame containing the training features.
#' @param labels A vector (usually a factor) of target labels.
#' @param ... Additional arguments passed to `PPtreeExt::PPtreeExtclass()`.
#'
#' @return A fitted `PPtreeExtclass` object.
#' @export
fit_adapter.classbound_pptree <- function(object, data, labels, ...) {
  if (!requireNamespace("PPtreeExt", quietly = TRUE)) {
    stop("Package 'PPtreeExt' must be installed to use method = 'pptree'.", call. = FALSE)
  }

  # Combine data and labels into a single formula for PPtreeExtclass
  train_data <- data
  train_data$.label. <- as.factor(labels)

  # Fit the model
  PPtreeExt::PPtreeExtclass(.label. ~ ., data = train_data, ...)
}

#' Predict using a fitted PPtreeExt model
#'
#' @description Adapter function to generate standardized predictions from a PPtreeExt model.
#'
#' @param model A fitted `PPtreeExtclass` object.
#' @param newdata A data frame of new observations to predict on.
#' @param ... Additional arguments passed to `predict()`.
#'
#' @return A list containing `class` (predicted labels) and `probs` (probabilities, NULL for pptree).
#' @importFrom stats predict
#' @export
predict_model.classbound_pptree <- function(model, newdata, ...) {
  raw_model <- model$model
  if (!inherits(raw_model, "PPtreeExtclass")) {
    stop("Model must be a 'PPtreeExtclass' object.", call. = FALSE)
  }

  # Predict classes
  preds_raw <- predict(raw_model, newdata, ...)
  preds <- as.factor(preds_raw$predict.class)

  # Probability output is typically not provided natively by PPtreeExtclass
  probs <- NULL

  list(
    class = preds,
    probs = probs
  )
}
