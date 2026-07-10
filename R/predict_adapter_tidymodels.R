#' @export
predict_adapter.workflow <- function(model, newdata, ...) {
  if (!requireNamespace("workflows", quietly = TRUE)) {
    stop("Package 'workflows' is required to process workflow models.")
  }

  preds_class <- predict(model, newdata, type = "class")[[".pred_class"]]

  preds_prob <- tryCatch(
    {
      as.matrix(predict(model, newdata, type = "prob"))
    },
    error = function(e) {
      NULL
    }
  )

  list(class = preds_class, probs = preds_prob)
}

#' @export
predict_adapter.model_fit <- function(model, newdata, ...) {
  if (!requireNamespace("parsnip", quietly = TRUE)) {
    stop("Package 'parsnip' is required to process model_fit models.")
  }

  preds_class <- predict(model, newdata, type = "class")[[".pred_class"]]

  preds_prob <- tryCatch(
    {
      as.matrix(predict(model, newdata, type = "prob"))
    },
    error = function(e) {
      NULL
    }
  )

  list(class = preds_class, probs = preds_prob)
}
