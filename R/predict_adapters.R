#' Predict using a fitted rpart model
#'
#' @description Adapter function to generate standardized predictions from an rpart model.
#'
#' @param model A fitted `rpart` object.
#' @param newdata A data frame of new observations to predict on.
#' @param ... Additional arguments passed to `predict()`.
#'
#' @return A list containing `class` (predicted labels) and `probs` (probabilities).
#' @keywords internal
#' @importFrom stats predict
#' @export
predict_adapter.rpart <- function(model, newdata, ...) {
  preds <- predict(model, newdata, type = "class", ...)

  probs <- predict(model, newdata, type = "prob", ...)

  list(
    class = preds,
    probs = probs
  )
}

#' Predict using a fitted randomForest model
#'
#' @description Adapter function to generate standardized predictions from a randomForest model.
#'
#' @param model A fitted `randomForest` object.
#' @param newdata A data frame of new observations to predict on.
#' @param ... Additional arguments passed to `predict()`.
#'
#' @return A list containing `class` (predicted labels) and `probs` (probabilities).
#' @keywords internal
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

#' Predict using a fitted PPtreeExt model
#'
#' @description Adapter function to generate standardized predictions from a PPtreeExt model.
#'
#' @param model A fitted `PPtreeExtclass` object.
#' @param newdata A data frame of new observations to predict on.
#' @param ... Additional arguments passed to `predict()`.
#'
#' @return A list containing `class` (predicted labels) and `probs` (probabilities, NULL for pptree).
#' @keywords internal
#' @importFrom stats predict
#' @export
predict_adapter.PPtreeExtclass <- function(model, newdata, ...) {
  preds_raw <- predict(model, newdata, ...)
  if (is.null(preds_raw$predict.class)) {
    stop("Unexpected prediction output from PPtreeExt model: 'predict.class' not found.", call. = FALSE)
  }
  preds <- as.factor(preds_raw$predict.class)

  probs <- NULL

  list(
    class = preds,
    probs = probs
  )
}

#' Predict using a fitted PPtreeViz model
#'
#' @description Adapter function to generate standardized predictions from a PPtreeViz model.
#'
#' @param model A fitted \code{PPtreeclass} object.
#' @param newdata A data frame of new observations to predict on.
#' @param ... Additional arguments passed to \code{predict()}.
#'
#' @return A list containing \code{class} (predicted labels) and \code{probs} (probabilities, NULL for PPtreeViz).
#' @keywords internal
#' @importFrom stats predict
#' @export
predict_adapter.PPtreeclass <- function(model, newdata, ...) {
  newdata_df <- as.data.frame(newdata)

  preds <- as.factor(predict(model, newdata = newdata_df, ...))

  probs <- NULL

  list(
    class = preds,
    probs = probs
  )
}

#' Predict using a fitted ppforest2 model
#'
#' @description Adapter function to generate standardized predictions from a ppforest2 model.
#'
#' \strong{Developer Note:} During development, an issue was observed with
#' \code{ppforest2} v0.1.3 when training data is contiguous by class but
#' ordered so that the lowest class does not occur first (e.g.,
#' \code{Class 3, Class 3, Class 1, Class 1}). In this case, the underlying
#' C++ code can produce \code{Grouping::init: partition must be rooted at row 0}.
#' \code{classbound} does not reorder or otherwise modify the training data
#' specifically to avoid this behavior. If a future version of \code{ppforest2}
#' changes this behavior, no corresponding change to \code{classbound} is
#' expected to be necessary.
#' @param model A fitted \code{pprf_classification} object.
#' @param newdata A data frame of new observations to predict on.
#' @param ... Additional arguments passed to \code{predict()}.
#'
#' @return A list containing \code{class} (predicted labels) and \code{probs} (probabilities).
#' @keywords internal
#' @importFrom stats predict
#' @export
predict_adapter.pprf_classification <- function(model, newdata, ...) {
  newdata_df <- as.data.frame(newdata)

  expected_preds <- colnames(model$x)
  if (!is.null(expected_preds) && all(expected_preds %in% colnames(newdata_df))) {
    new_data_mat <- as.matrix(newdata_df[, expected_preds, drop = FALSE])
  } else {
    # Fallback if expected_preds cannot be resolved
    target <- all.vars(model$formula)[1]
    if (!is.null(target) && !(target %in% colnames(newdata_df))) {
      dummy_val <- if (!is.null(model$groups)) model$groups[1] else 1
      newdata_df[[target]] <- dummy_val
    }
    new_data_mat <- newdata_df
  }

  preds <- predict(model, new_data_mat, type = "class", ...)
  probs <- as.matrix(predict(model, new_data_mat, type = "prob", ...))

  list(
    class = preds,
    probs = probs
  )
}

#' @export
predict_adapter.workflow <- function(model, newdata, ...) {
  if (!requireNamespace("workflows", quietly = TRUE)) {
    stop("Package 'workflows' is required to process workflow models.")
  }

  preds_class <- predict(model, newdata, type = "class")[[".pred_class"]]

  preds_prob <- tryCatch(
    {
      p_mat <- as.matrix(predict(model, newdata, type = "prob"))
      colnames(p_mat) <- gsub("^\\.pred_", "", colnames(p_mat))
      p_mat
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
      p_mat <- as.matrix(predict(model, newdata, type = "prob"))
      colnames(p_mat) <- gsub("^\\.pred_", "", colnames(p_mat))
      p_mat
    },
    error = function(e) {
      NULL
    }
  )

  list(class = preds_class, probs = preds_prob)
}
