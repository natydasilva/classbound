#' Preprocess data before model fitting
#'
#' @description Handles common preprocessing steps like validating inputs, converting factors, 
#' handling missing values, and optional scaling. This centralizes validation so classifier 
#' adapters can assume validated inputs.
#'
#' @param data A data frame of raw training or prediction data.
#' @param labels Optional. A vector of target labels corresponding to `data`. If provided, it is converted to a factor with unused levels dropped.
#' @param scale Logical. If `TRUE`, scales numeric columns in `data` to mean 0 and standard deviation 1. Defaults to `FALSE`. Note that this scaling is not persisted for predictions; it is an isolated transformation.
#' @param ... Additional processing arguments (currently unused).
#'
#' @return A list containing the validated and processed `data` and `labels`.
#' @export
preprocess_data <- function(data, labels = NULL, scale = FALSE, ...) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  
  # Coerce to a standard data.frame to remove tibble or other sub-classes
  # This prevents indexing bugs in older packages that expect df[, col] to return a vector
  data <- as.data.frame(data)

  if (!is.null(labels)) {
    if (length(labels) != nrow(data)) {
      stop("Length of `labels` must match the number of rows in `data`.", call. = FALSE)
    }
  }

  if (any(is.na(data)) || (!is.null(labels) && any(is.na(labels)))) {
    stop("Missing values are not supported. Please handle them before modeling.", call. = FALSE)
  }

  # Factor handling for labels
  if (!is.null(labels)) {
    labels <- droplevels(as.factor(labels))
  }

  # Factor handling for character columns in data
  # Character columns are automatically converted to factors to ensure predictable modeling behavior.
  char_cols <- sapply(data, is.character)
  if (any(char_cols)) {
    data[char_cols] <- lapply(data[char_cols], as.factor)
  }

  # Drop unused levels for all factor columns in data
  factor_cols <- sapply(data, is.factor)
  if (any(factor_cols)) {
    data[factor_cols] <- lapply(data[factor_cols], droplevels)
  }

  # Optional scaling for numeric columns
  if (scale) {
    num_cols <- sapply(data, is.numeric)
    if (any(num_cols)) {
      data[num_cols] <- lapply(data[num_cols], function(x) as.numeric(scale(x)))
    }
  }

  list(data = data, labels = labels)
}
