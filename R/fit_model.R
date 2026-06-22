#' Fit a machine learning model
#'
#' @description Fits a machine learning classifier using a unified interface.
#'
#' @param data A data frame containing the training features.
#' @param labels A vector of target labels corresponding to `data`.
#' @param method A string specifying the classifier method (e.g., "rpart").
#' @param ... Additional arguments passed to the specific model adapter.
#'
#' @return A fitted model object with a normalized structure of class "classbound_model", containing the model, method, and extracted feature metadata.
#' @export
fit_model <- function(data, labels, method, ...) {
  if (missing(method)) {
    stop("Please specify a classification method.", call. = FALSE)
  }

  if (is.null(labels)) {
    stop("Please provide classification labels.", call. = FALSE)
  }

  if (is.data.frame(data) && ".label." %in% colnames(data)) {
    stop("Column name '.label.' is reserved for internal use. Please rename this column.", call. = FALSE)
  }

  # Convert method string into an S3 method spec for dispatch
  method_spec <- structure(method, class = c(paste0("classbound_", method), "classbound_method"))

  # Centralized validation and preprocessing
  processed <- preprocess_data(data, labels)
  data <- processed$data
  labels <- processed$labels

  # Route fitting to specific adapter using S3 dispatch
  model_fit <- fit_adapter(method_spec, data, labels, ...)

  # Extract feature metadata
  feature_meta <- list(
    names = colnames(data),
    types = lapply(data, class),
    levels = lapply(data, levels)
  )

  # Wrap the fitted model and metadata to maintain context for prediction
  structure(
    list(
      model = model_fit,
      method = method,
      features = feature_meta
    ),
    class = c(paste0("classbound_", method), "classbound_model")
  )
}

#' Internal generic for fitting classifier adapters
#'
#' @description This generic is exported for S3 dispatch but is not intended for direct use.
#' Use \code{\link{fit_model}} instead. Adapter authors should implement S3 methods
#' for this generic (e.g., \code{fit_adapter.classbound_mymodel}).
#'
#' @param object An S3 object defining the method spec.
#' @param data A data frame containing the training features.
#' @param labels A vector of target labels.
#' @param ... Additional arguments.
#'
#' @keywords internal
#' @export
fit_adapter <- function(object, data, labels, ...) {
  UseMethod("fit_adapter")
}

#' @export
fit_adapter.default <- function(object, data, labels, ...) {
  stop(sprintf("Classifier method '%s' is not supported.", as.character(object)), call. = FALSE)
}
