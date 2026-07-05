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
  # Ensure newdata only contains features used in model
  # The training data was saved in the model but without .label.
  # But PPforest predict requires just the features.

  # remove any non-numeric columns if PPforest expects it,
  # or just pass newdata directly since our pipeline guarantees matching features.
  newdata_df <- as.data.frame(newdata)

  # PPforest has a known internal bug where foreach/codetools emits a false-positive
  # warning about `...` being used in an incorrect context. We muffle only this specific warning.
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

  # The predicted classes are in the 3rd element of the list
  preds <- as.factor(preds_raw[[3]])

  list(
    class = preds,
    probs = NULL
  )
}
