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
  if (!inherits(model, "classbound_model")) {
    stop("Model must be a 'classbound_model' object returned from fit_model().", call. = FALSE)
  }
  
  res <- predict_adapter(model, newdata, ...)
  
  if (!is.null(model$class_levels)) {
    res$class <- factor(res$class, levels = model$class_levels)
  } else if (!is.factor(res$class)) {
    res$class <- as.factor(res$class)
  }
  
  res
}

#' Internal generic for predicting classifier adapters
#'
#' @description This generic is exported for S3 dispatch but is not intended for direct use.
#' Use \code{\link{predict_model}} instead.
#'
#' @param model A fitted model object.
#' @param newdata A data frame of new observations.
#' @param ... Additional arguments.
#'
#' @keywords internal
#' @export
predict_adapter <- function(model, newdata, ...) {
  UseMethod("predict_adapter")
}

#' @export
predict_adapter.default <- function(model, newdata, ...) {
  stop(sprintf("Prediction for method '%s' is not supported.", model$method), call. = FALSE)
}
