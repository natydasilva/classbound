#' Fit a machine learning model
#'
#' @description Fits a machine learning classifier using a standard interface.
#'
#' @param data A data frame containing the training features and response.
#' @param formula A formula specifying the response and predictors.
#' @param classifier The classification function to use (e.g., \code{rpart::rpart}).
#' @param interface A string specifying how to invoke the classifier.
#' \itemize{
#'   \item \code{"formula"} (default): For classifiers that accept a formula and data frame natively (verified for \code{rpart::rpart}, \code{e1071::svm}, \code{stats::glm}).
#'   \item \code{"matrix"}: For classifiers that expect a predictor matrix \code{x} and response vector \code{y} (verified for \code{randomForest::randomForest} default method).
#'   \item \code{"custom"}: For non-standard APIs where you must pass all arguments explicitly via \code{fit_args} (verified for \code{qeML::qeKNN}).
#' }
#' @param fit_args A named list of additional arguments passed to the classifier during fitting.
#'
#' @return A fitted model object with a normalized structure of class "classbound", containing the raw model and extracted feature metadata.
#' @importFrom stats model.frame model.matrix model.response
#' @export
fit_model <- function(data, formula, classifier, interface = c("formula", "matrix", "custom"), fit_args = list()) {
  if (missing(classifier)) {
    stop("Please specify a classifier function.", call. = FALSE)
  }
  if (missing(formula)) {
    stop("Please specify a formula.", call. = FALSE)
  }

  interface <- match.arg(interface)

  # Extract labels based on formula
  # model.frame handles NAs implicitly, but our preprocess_data strictly rejects them.
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)

  if (is.null(y)) {
    stop("Could not extract response variable from formula.", call. = FALSE)
  }

  # Capture original class levels before preprocessing drops unused levels
  orig_class_levels <- if (is.factor(y)) levels(y) else sort(unique(as.character(y)))

  # Centralized validation and preprocessing
  processed <- preprocess_data(data, y)
  data <- processed$data
  y <- processed$labels

  # Update mf with preprocessed data
  mf <- model.frame(formula, data = data)
  # Route fitting based on selected interface
  if (interface == "formula") {
    model_fit <- do.call(classifier, c(list(formula = formula, data = mf), fit_args))
  } else if (interface == "matrix") {
    x_matrix <- model.matrix(formula, data = mf)[, -1, drop = FALSE] # remove intercept
    model_fit <- do.call(classifier, c(list(x = x_matrix, y = y), fit_args))
  } else if (interface == "custom") {
    model_fit <- do.call(classifier, fit_args)
  }

  # Extract feature metadata using the processed model frame (excluding response)
  response_var <- all.vars(formula[[2]])
  predictors_df <- mf[, setdiff(colnames(mf), response_var), drop = FALSE]
  feature_meta <- list(
    names = colnames(predictors_df),
    types = lapply(predictors_df, class),
    levels = lapply(predictors_df, levels)
  )

  # Wrap the fitted model and metadata in the single public 'classbound' class
  structure(
    list(
      fit = model_fit,
      features = feature_meta,
      class_levels = orig_class_levels
    ),
    class = "classbound"
  )
}
