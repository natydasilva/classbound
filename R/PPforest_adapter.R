#' Fit a PPforest model
#'
#' @description Adapter function to fit a Projection Pursuit Random Forest.
#'
#' @param object An S3 method spec object.
#' @param data A data frame containing the training features.
#' @param labels A vector (usually a factor) of target labels.
#' @param ... Additional arguments passed to \code{PPforest::PPforest()}.
#'
#' @return A fitted \code{PPforest} object.
#' @export
fit_adapter.classbound_PPforest <- function(object, data, labels, ...) {
  if (!requireNamespace("PPforest", quietly = TRUE)) {
    stop("Package 'PPforest' must be installed to use method = 'PPforest'.", call. = FALSE)
  }

  train_data <- as.data.frame(data)
  train_data$.label. <- as.factor(labels)


  # Extract dot args
  args <- list(...)
  if (is.null(args$size.tr)) args$size.tr <- 1
  if (is.null(args$m)) args$m <- 50
  if (is.null(args$size.p)) args$size.p <- 1
  if (is.null(args$PPmethod)) args$PPmethod <- "LDA"

  args$data <- train_data
  args$y <- ".label."

  do.call(PPforest::PPforest, args)
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
predict_model.classbound_PPforest <- function(model, newdata, ...) {
  raw_model <- model$model
  if (!inherits(raw_model, "PPforest")) {
    stop("Model must be a 'PPforest' object.", call. = FALSE)
  }

  # Ensure newdata only contains features used in model
  # The training data was saved in the model but without .label.
  # But PPforest predict requires just the features.

  # remove any non-numeric columns if PPforest expects it,
  # or just pass newdata directly since our pipeline guarantees matching features.
  newdata_df <- as.data.frame(newdata)

  # PPforest has a known internal bug where foreach/codetools emits a false-positive
  # warning about `...` being used in an incorrect context. We muffle only this specific warning.
  preds_raw <- withCallingHandlers(
    predict(raw_model, newdata = newdata_df, ...),
    warning = function(w) {
      if (grepl("may be used in an incorrect context", w$message)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  # The predicted classes are in the 3rd element of the list
  preds <- as.factor(preds_raw[[3]])

  list(
    class = preds,
    probs = NULL
  )
}
