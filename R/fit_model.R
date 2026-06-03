#' Fit a machine learning model
#'
#' @description Fits a machine learning classifier using a unified interface.
#'
#' @param data A data frame containing the training features.
#' @param labels A vector of target labels corresponding to `data`.
#' @param method A string specifying the classifier method (e.g., "rpart").
#' @param ... Additional arguments passed to the specific model adapter.
#'
#' @return A fitted model object with a normalized structure of class "classbound_model".
#' @export
fit_model <- function(data, labels, method, ...) {
  if (missing(method)) {
    stop("Please specify a classification method.", call. = FALSE)
  }

  # Route fitting to specific adapter based on method
  model_fit <- switch(
    method,
    "rpart" = fit_rpart(data, labels, ...),
    stop(sprintf("Classifier method '%s' is not supported.", method), call. = FALSE)
  )

  # Wrap the fitted model and metadata to maintain context for prediction
  structure(
    list(
      model = model_fit,
      method = method
    ),
    class = "classbound_model"
  )
}
