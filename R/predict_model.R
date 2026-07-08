#' Predict using a classbound model
#'
#' @description Generates predictions using a unified interface across all classifiers.
#' Dispatches natively on objects of class \code{"classbound"}.
#'
#' @param object A fitted classbound model. This is the object returned by \code{fit_model()} or \code{classbound()}.
#' @param newdata A data frame of new observations to predict on.
#' @param predict_args A named list of additional arguments passed to \code{predict_adapter}.
#' @param predfun A custom function to generate predictions for non-standard models.
#' @param ... Additional arguments passed to the specific model adapter.
#'
#' @return A list containing \code{class} (a factor of predicted labels) and \code{probs}
#' (a probability matrix, or strictly \code{NULL} if the classifier lacks probability support).
#' Downstream functions like \code{boundary_compute()} are designed to handle \code{probs = NULL} gracefully.
#' @export
predict.classbound <- function(object, newdata, predict_args = list(), predfun = NULL, ...) {
  if (!is.null(predfun)) {
    # Use the user-supplied custom prediction function
    args <- c(list(object$fit, newdata), predict_args, list(...))
    preds_raw <- do.call(predfun, args)

    if (is.list(preds_raw) && !is.null(preds_raw$class)) {
      res <- preds_raw
    } else {
      res <- list(class = as.factor(preds_raw), probs = NULL)
    }
  } else {
    # Call predict_adapter on the native model, dispatching on its native class
    args <- c(list(model = object$fit, newdata = newdata), predict_args, list(...))
    res <- do.call(predict_adapter, args)
  }

  if (!is.null(object$class_levels)) {
    res$class <- factor(res$class, levels = object$class_levels)
  } else if (!is.factor(res$class)) {
    res$class <- as.factor(res$class)
  }

  res
}

#' Predict using a fitted classbound model
#'
#' @description Generates predictions using a unified interface across all classifiers.
#' This function is a compatibility wrapper around the standard \code{predict()} method for "classbound" objects.
#'
#' @param model A fitted classbound model. This corresponds to the \code{object} argument used by the standard R \code{predict()} generic. This wrapper simply calls \code{predict(model, ...)}.
#' @param newdata A data frame of new observations to predict on.
#' @param predict_args A named list of additional arguments passed to \code{predict_adapter}.
#' @param predfun A custom function to generate predictions for non-standard models.
#' @param ... Additional arguments passed to the specific model adapter.
#'
#' @return A list containing \code{class} (a factor of predicted labels) and \code{probs}
#' (a probability matrix, or strictly \code{NULL} if the classifier lacks probability support).
#' @export
predict_model <- function(model, newdata, predict_args = list(), predfun = NULL, ...) {
  if (!inherits(model, "classbound")) {
    stop("Model must be a 'classbound' object returned from fit_model().", call. = FALSE)
  }
  predict(model, newdata = newdata, predict_args = predict_args, predfun = predfun, ...)
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
  # Standard R predict convention: predict(model) usually returns a vector or factor of classes
  # We wrap this in a tryCatch to provide an informative error if it fails
  preds <- tryCatch(
    predict(model, newdata, ...),
    error = function(e) {
      stop(sprintf(
        "Default prediction failed for native class '%s'. Original error: %s",
        class(model)[1], e$message
      ), call. = FALSE)
    }
  )

  # Explicitly guard against non-standard outputs (e.g. lists, data frames, matrices)
  # that some models return. These require specific predict_adapter implementations.
  if (is.list(preds) || is.matrix(preds) || is.data.frame(preds)) {
    stop(sprintf(
      "The model returned a non-standard prediction object (class '%s' returned %s). Please supply `predfun` to extract the desired class predictions.",
      class(model)[1], class(preds)[1]
    ), call. = FALSE)
  }

  list(class = as.factor(preds), probs = NULL)
}
