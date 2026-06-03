#' Predict using a fitted model
#'
#' @description Generates predictions using a unified interface across all classifiers.
#'
#' @param model A fitted model object typically returned by `fit_model`.
#' @param newdata A data frame of new observations to predict on.
#' @param ... Additional arguments passed to the specific model adapter.
#'
#' @return A list containing `class` (predicted labels) and `probs` (probabilities).
#' @export
predict_model <- function(model, newdata, ...) {
  if (!inherits(model, "classbound_model")) {
    stop("Model must be a 'classbound_model' object returned from fit_model().", call. = FALSE)
  }

  # Route prediction to specific adapter based on method
  preds <- switch(
    model$method,
    "rpart" = predict_rpart(model$model, newdata, ...),
    stop(sprintf("Prediction for method '%s' is not supported.", model$method), call. = FALSE)
  )

  preds
}
