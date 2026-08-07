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

#' Predict using a fitted PPforest model
#'
#' @description Adapter function to generate standardized predictions from a PPforest model.
#'
#' @param model A fitted \code{PPforest} object.
#' @param newdata A data frame of new observations to predict on.
#' @param ... Additional arguments passed to \code{predict()}.
#'
#' @return A list containing \code{class} (predicted labels) and \code{probs} (probabilities, NULL for PPforest).
#' @importFrom stats predict
#' @export
predict_adapter.PPforest <- function(model, newdata, ...) {
  # Subset newdata to model features.

  newdata_df <- as.data.frame(newdata)

  # Suppress PPforest codetools false-positive warning.
  preds_raw <- withCallingHandlers(
    predict(model, newdata = newdata_df, ...),
    warning = function(w) {
      if (grepl("may be used in an incorrect context", w$message)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  # Validate prediction output structure
  if (!is.list(preds_raw) || length(preds_raw) < 3) {
    stop("Unexpected prediction output from PPforest model.", call. = FALSE)
  }

  preds <- as.factor(preds_raw[[3]])

  list(
    class = preds,
    probs = NULL
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
      as.matrix(predict(model, newdata, type = "prob"))
    },
    error = function(e) {
      NULL
    }
  )

  list(class = preds_class, probs = preds_prob)
}
