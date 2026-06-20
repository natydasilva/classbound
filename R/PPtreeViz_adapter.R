#' Fit a PPtreeViz model
#'
#' @description Adapter function to fit a Projection Pursuit classification tree (PPtreeViz).
#'
#' @param data A data frame containing the training features.
#' @param labels A vector (usually a factor) of target labels.
#' @param ... Additional arguments passed to \code{PPtreeViz::PPTreeclass()}.
#'
#' @return A fitted \code{PPtreeclass} object.
#' @export
fit_PPtreeViz <- function(data, labels, ...) {
  if (!requireNamespace("PPtreeViz", quietly = TRUE)) {
    stop("Package 'PPtreeViz' must be installed to use method = 'PPtreeViz'.", call. = FALSE)
  }

  train_data <- as.data.frame(data)
  train_data$.label. <- as.factor(labels)

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
predict_PPtreeViz <- function(model, newdata, ...) {
  if (!inherits(model, "PPtreeclass")) {
    stop("Model must be a 'PPtreeclass' object.", call. = FALSE)
  }

  newdata_df <- as.data.frame(newdata)
  
  # PPtreeViz predict requires only the feature columns. We drop non-numeric columns.
  # Assuming the prediction grid features are all numeric.
  numeric_cols <- vapply(newdata_df, is.numeric, logical(1))
  newdata_numeric <- newdata_df[, numeric_cols, drop = FALSE]

  # Predict classes
  preds <- predict(model, newdata = newdata_numeric, ...)

  # PPtreeViz does not support probability predictions natively
  probs <- NULL

  list(
    class = preds,
    probs = probs
  )
}
