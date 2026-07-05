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
predict_adapter.randomForest <- function(model, newdata, ...) {
  preds <- predict(model, newdata, type = "response", ...)
  probs <- predict(model, newdata, type = "prob", ...)

  list(
    class = preds,
    probs = probs
  )
}
