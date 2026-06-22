#' Predict using a fitted model
#'
#' @description Generates predictions using a unified interface across all classifiers.
#'
#' @param model A fitted model object typically returned by `fit_model`.
#' @param newdata A data frame of new observations to predict on.
#' @param ... Additional arguments passed to the specific model adapter.
#'
#' @return A list containing `class` (a factor of predicted labels) and `probs` 
#' (a probability matrix, or strictly `NULL` if the classifier lacks probability support).
#' Downstream functions like `boundary_compute()` are designed to handle `probs = NULL` gracefully.
#' @export
predict_model <- function(model, newdata, ...) {
  UseMethod("predict_model")
}

#' @export
predict_model.default <- function(model, newdata, ...) {
  if (!inherits(model, "classbound_model")) {
    stop("Model must be a 'classbound_model' object returned from fit_model().", call. = FALSE)
  }
  stop(sprintf("Prediction for method '%s' is not supported.", model$method), call. = FALSE)
}
