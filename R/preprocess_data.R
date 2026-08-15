#' Preprocess data before model fitting
#'
#' @description
#' Applies standard preprocessing steps to a training dataset: validates structure,
#' coerces character columns to factors, drops unused factor levels, rejects missing
#' and infinite values, and converts response labels to a factor.
#'
#' @details
#' This function is called automatically by `fit_model()`. **Do not call it manually
#' before `fit_model()`** — doing so will result in feature metadata being extracted
#' from the pre-processed data rather than the original data, which corrupts the
#' imputation values used by `boundary_compute()` for 2D slicing.
#'
#' It is exported for use in custom workflows and testing, but it is not typically
#' needed in the standard `fit_model()` → `boundary_compute()` pipeline.
#'
#' @param data A data frame of raw training data.
#' @param labels Optional vector of target labels. Converted to a factor with
#'   unused levels dropped.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with two elements:
#'   - `$data`: the preprocessed data frame
#'   - `$labels`: the processed factor labels (or `NULL` if not supplied)
#' @export
preprocess_data <- function(data, labels = NULL, ...) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  # Coerce to a standard data.frame to remove tibble or other sub-classes
  # This prevents indexing bugs in older packages that expect df[, col] to return a vector
  data <- as.data.frame(data)

  if (nrow(data) < 2) {
    stop("`data` must have at least two rows to train a classification boundary.", call. = FALSE)
  }

  if (any(duplicated(colnames(data)))) {
    stop("Duplicate column names found in `data`.", call. = FALSE)
  }

  if (!is.null(labels)) {
    if (length(labels) != nrow(data)) {
      stop("Length of `labels` must match the number of rows in `data`.", call. = FALSE)
    }
  }

  if (any(is.na(data)) || (!is.null(labels) && any(is.na(labels)))) {
    stop("Missing values are not supported. Please handle them before modeling.", call. = FALSE)
  }

  # Reject Infinite values which crash downstream C/C++ model backends
  is_inf <- function(x) is.numeric(x) && any(is.infinite(x))
  if (any(sapply(data, is_inf))) {
    stop("Infinite values are not supported. Please handle them before modeling.", call. = FALSE)
  }

  # Factor handling for labels
  if (!is.null(labels)) {
    labels <- droplevels(as.factor(labels))
  }

  # Coerce character columns to factors.
  char_cols <- sapply(data, is.character)
  if (any(char_cols)) {
    data[char_cols] <- lapply(data[char_cols], as.factor)
  }

  # Drop unused factor levels.
  factor_cols <- sapply(data, is.factor)
  if (any(factor_cols)) {
    data[factor_cols] <- lapply(data[factor_cols], droplevels)
  }


  list(data = data, labels = labels)
}
