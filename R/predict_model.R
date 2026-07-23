#' Predict using a classbound model
#'
#' @description Generates predictions using a unified interface across all classifiers.
#' Dispatches natively on objects of class \code{"classbound"}. Attempting to call
#' \code{predict()} directly on a \code{"classbound_multi"} object will result in an error,
#' as multi-model boundaries are evaluated internally by \code{boundary_compute()}.
#'
#' @param object A fitted classbound model. This is the object returned by \code{fit_model()} or \code{classbound()}.
#' @param newdata A data frame of new observations to predict on.
#' @param predict_args A named list of additional arguments passed to \code{predict_adapter}.
#' @param predfun A custom function to generate predictions for non-standard models.
#'   The function must accept at least two arguments: \code{model} (the fitted native model) 
#'   and \code{newdata} (a data frame of new observations). It should return either a vector/factor 
#'   of predicted classes, or a list containing \code{class} (predicted labels) and \code{probs} (a probability matrix).
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
    # Isolate strictly the features the model was trained on
    # This prevents fragile native predict methods (like PPtreeViz) from crashing
    # when non-numeric target variables or extra metadata columns are passed in newdata.
    if (!is.null(object$metadata$features$names)) {
      expected_feats <- object$metadata$features$names
      present_feats <- intersect(expected_feats, colnames(newdata))
      if (length(present_feats) > 0) {
        newdata <- newdata[, present_feats, drop = FALSE]
      }
    }

    # Call predict_adapter on the native model, dispatching on its native class
    args <- c(list(model = object$fit, newdata = newdata), predict_args, list(...))
    res <- do.call(predict_adapter, args)
  }

  if (!is.null(object$metadata$class_levels)) {
    res$class <- factor(res$class, levels = object$metadata$class_levels)
  } else if (!is.factor(res$class)) {
    res$class <- as.factor(res$class)
  }

  res
}

#' @rdname predict.classbound
#' @export
predict.classbound_multi <- function(object, newdata, predict_args = list(), predfun = NULL, ...) {
  stop("Cannot call predict() directly on a multi-model boundary comparison object.", call. = FALSE)
}

#' Predict using a fitted classbound model
#'
#' @description Generates predictions using a unified interface across all classifiers.
#' This function is a compatibility wrapper around the standard \code{predict()} method for "classbound" objects.
#'
#' @param model A fitted classbound model. This corresponds to the \code{object} argument used by the
#'   standard R \code{predict()} generic. This wrapper simply calls \code{predict(model, ...)}.
#' @param newdata A data frame of new observations to predict on.
#' @param predict_args A named list of additional arguments passed to \code{predict_adapter}.
#' @param predfun A custom function to generate predictions for non-standard models.
#'   The function must accept at least two arguments: \code{model} (the fitted native model) 
#'   and \code{newdata} (a data frame of new observations). It should return either a vector/factor 
#'   of predicted classes, or a list containing \code{class} (predicted labels) and \code{probs} (a probability matrix).
#' @param ... Additional arguments passed to the specific model adapter.
#'
#' @return A list containing \code{class} (a factor of predicted labels) and \code{probs}
#' (a probability matrix, or strictly \code{NULL} if the classifier lacks probability support).
#' @examples
#' \donttest{
#' library(palmerpenguins)
#' data(penguins)
#' peng_data <- na.omit(penguins[, c("species", "bill_length_mm", "bill_depth_mm")])
#' 
#' m_rpart <- fit_model(peng_data, species ~ ., rpart::rpart)
#' preds <- predict_model(m_rpart, newdata = peng_data[1:5, ])
#' }
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
    stop(
      sprintf(
        paste0(
          "The model returned a non-standard prediction object (class '%s' returned %s). ",
          "Please supply `predfun` to extract the desired class predictions."
        ),
        class(model)[1], class(preds)[1]
      ),
      call. = FALSE
    )
  }

  list(class = as.factor(preds), probs = NULL)
}
