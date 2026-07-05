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

  # Predict classes
  preds <- as.factor(predict(model, newdata = newdata_df, ...))

  # PPtreeViz does not support probability predictions natively
  probs <- NULL

  list(
    class = preds,
    probs = probs
  )
}
