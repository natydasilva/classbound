#' Fit a PPforest model
#'
#' @description Adapter function to fit a Projection Pursuit Random Forest.
#'
#' @param data A data frame containing the training features.
#' @param labels A vector (usually a factor) of target labels.
#' @param ... Additional arguments passed to \code{PPforest::PPforest()}.
#'
#' @return A fitted \code{PPforest} object.
#' @export
fit_PPforest <- function(data, labels, ...) {
  if (!requireNamespace("PPforest", quietly = TRUE)) {
    stop("Package 'PPforest' must be installed to use method = 'PPforest'.", call. = FALSE)
  }

  train_data <- as.data.frame(data)
  train_data$.label. <- as.factor(labels)

  # Determine number of features for default size.p
  n_features <- ncol(data)

  # Extract dot args
  args <- list(...)
  if (is.null(args$size.tr)) args$size.tr <- 1
  if (is.null(args$m)) args$m <- 50
  if (is.null(args$size.p)) args$size.p <- 1
  if (is.null(args$PPmethod)) args$PPmethod <- 'LDA'

  args$data <- train_data
  args$y <- '.label.'

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
predict_PPforest <- function(model, newdata, ...) {
  if (!inherits(model, "PPforest")) {
    stop("Model must be a 'PPforest' object.", call. = FALSE)
  }

  # Ensure newdata only contains features used in model
  # The training data was saved in the model but without .label.
  # But PPforest predict requires just the features.
  
  # remove any non-numeric columns if PPforest expects it, 
  # or just pass newdata directly since our pipeline guarantees matching features.
  newdata_df <- as.data.frame(newdata)
  preds_raw <- predict(model, newdata = newdata_df, ...)
  
  # The predicted classes are in the 3rd element of the list
  preds <- as.factor(preds_raw[[3]])

  list(
    class = preds,
    probs = NULL
  )
}
