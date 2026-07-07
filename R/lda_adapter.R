#' @export
predict_adapter.lda <- function(model, newdata, ...) {
  # lda from MASS returns a list containing $class, $posterior, and $x
  preds <- predict(model, newdata, ...)
  list(class = as.factor(preds$class), probs = NULL)
}
