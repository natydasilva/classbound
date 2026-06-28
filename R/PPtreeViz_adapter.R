#' Fit a PPtreeViz model
#'
#' @description Adapter function to fit a Projection Pursuit classification tree (PPtreeViz).
#'
#' @param object An S3 method spec object.
#' @param data A data frame containing the training features.
#' @param labels A vector (usually a factor) of target labels.
#' @param ... Additional arguments passed to \code{PPtreeViz::PPTreeclass()}.
#'
#' @return A fitted \code{PPtreeclass} object.
#' @export
fit_adapter.classbound_PPtreeViz <- function(object, data, labels, ...) {
  if (!requireNamespace("PPtreeViz", quietly = TRUE)) {
    stop("Package 'PPtreeViz' must be installed to use method = 'PPtreeViz'.", call. = FALSE)
  }

  # Combine data and labels into a single formula for PPtreeViz
  train_data <- data
  train_data$.label. <- labels

  # Fit the model
  PPtreeViz::PPTreeclass(.label. ~ ., data = train_data, ...)
}

#' Predict using a fitted PPtreeViz model
#'
#' @description Adapter function to generate standardized predictions from a PPtreeViz model.
#'
#' @param model A fitted \code{PPtreeclass} object.
#' @param newdata A data frame of new observations to predict on.
#' @param ... Additional arguments passed to \code{predict()}.
#'
#' @return A list containing \code{class} (predicted labels) and \code{probs} (probabilities, NULL for PPtreeViz).
#' @importFrom stats predict
#' @export
predict_adapter.classbound_PPtreeViz <- function(model, newdata, ...) {
  raw_model <- model$model
  if (!inherits(raw_model, "PPtreeclass")) {
    stop("Model must be a 'PPtreeclass' object.", call. = FALSE)
  }

  newdata_df <- as.data.frame(newdata)

  # Predict classes
  preds <- as.factor(predict(raw_model, newdata = newdata_df, ...))

  # PPtreeViz does not support probability predictions natively
  probs <- NULL

  list(
    class = preds,
    probs = probs
  )
}
