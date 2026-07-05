#' Predict using a classbound model
#'
#' @description Generates predictions using a unified interface across all classifiers.
#' Dispatches natively on objects of class \code{"classbound"}.
#'
#' @param object A fitted model object returned by \code{fit_model} (of class \code{"classbound"}).
#' @param newdata A data frame of new observations to predict on.
#' @param predict_args A named list of additional arguments passed to \code{predict_adapter}.
#' @param ... Additional arguments passed to the specific model adapter.
#'
#' @return A list containing \code{class} (a factor of predicted labels) and \code{probs} 
#' (a probability matrix, or strictly \code{NULL} if the classifier lacks probability support).
#' Downstream functions like \code{boundary_compute()} are designed to handle \code{probs = NULL} gracefully.
#' @export
predict.classbound <- function(object, newdata, predict_args = list(), ...) {
  # Call predict_adapter on the native model, dispatching on its native class
  args <- c(list(model = object$fit, newdata = newdata), predict_args, list(...))
  res <- do.call(predict_adapter, args)
  
  if (!is.null(object$class_levels)) {
    res$class <- factor(res$class, levels = object$class_levels)
  } else if (!is.factor(res$class)) {
    res$class <- as.factor(res$class)
  }
  
  res
}

#' Predict using a fitted model (Compatibility Wrapper)
#'
#' @description Generates predictions using a unified interface across all classifiers.
#' This is a compatibility wrapper around \code{\link{predict.classbound}}.
#'
#' @param model A fitted model object typically returned by \code{fit_model} (of class \code{"classbound"}).
#' @param newdata A data frame of new observations to predict on.
#' @param predict_args A named list of additional arguments passed to \code{predict_adapter}.
#' @param ... Additional arguments passed to the specific model adapter.
#'
#' @return A list containing \code{class} (a factor of predicted labels) and \code{probs} 
#' (a probability matrix, or strictly \code{NULL} if the classifier lacks probability support).
#' @export
predict_model <- function(model, newdata, predict_args = list(), ...) {
  if (!inherits(model, "classbound")) {
    stop("Model must be a 'classbound' object returned from fit_model().", call. = FALSE)
  }
  predict(model, newdata = newdata, predict_args = predict_args, ...)
}

#' Internal generic for predicting classifier adapters
#'
#' @description This generic is exported for S3 dispatch but is not intended for direct use.
#' Use \code{\link{predict.classbound}} or \code{\link{predict_model}} instead.
#'
#' @param model A fitted native model object.
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
  stop(sprintf("Prediction adapter for native class '%s' is not implemented.", class(model)[1]), call. = FALSE)
}
