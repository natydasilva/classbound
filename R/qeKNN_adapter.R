#' @export
predict_adapter.qeKNN <- function(model, newdata, ...) {
  # qeKNN from qeML returns an object with predClasses for standard predictions
  preds <- predict(model, newdata, ...)
  list(class = as.factor(preds$predClasses), probs = NULL)
}
