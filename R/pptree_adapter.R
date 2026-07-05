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
predict_adapter.PPtreeExtclass <- function(model, newdata, ...) {
  # Predict classes
  preds_raw <- predict(model, newdata, ...)
  if (is.null(preds_raw$predict.class)) {
    stop("Unexpected prediction output from PPtreeExt model: 'predict.class' not found.", call. = FALSE)
  }
  preds <- as.factor(preds_raw$predict.class)

  # Probability output is typically not provided natively by PPtreeExtclass
  probs <- NULL

  list(
    class = preds,
    probs = probs
  )
}
