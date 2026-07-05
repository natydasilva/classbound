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
predict_adapter.rpart <- function(model, newdata, ...) {
  # Predict classes (discrete labels)
  preds <- predict(model, newdata, type = "class", ...)

  # Predict probabilities (continuous matrix)
  probs <- predict(model, newdata, type = "prob", ...)

  list(
    class = preds,
    probs = probs
  )
}
