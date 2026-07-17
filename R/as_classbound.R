#' Convert an existing fitted model into a classbound object
#'
#' @description Wraps a pre-fitted model with the metadata required by the `classbound`
#'   visualization pipeline. This enables "bring your own model" workflows and supports
#'   native integration with frameworks like `tidymodels`.
#'
#' @param model A fitted model object.
#' @param data A data frame containing the training data used to fit the model.
#'   This is strictly used to extract feature metadata (names, types, and levels).
#'   It is not passed through `preprocess_data()`.
#' @param response A string specifying the name of the response column in `data`.
#'   If provided, this column is excluded from feature metadata, and its levels
#'   are recorded as the class levels.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `classbound` object.
#' @export
as_classbound <- function(model, data, response = NULL, ...) {
  UseMethod("as_classbound")
}

#' @export
as_classbound.default <- function(model, data, response = NULL, ...) {
  if (missing(data)) {
    stop("`data` must be provided to extract feature metadata.", call. = FALSE)
  }

  if (!is.null(response) && response %in% colnames(data)) {
    y <- data[[response]]
    orig_class_levels <- if (is.factor(y)) levels(y) else sort(unique(as.character(y)))
    predictors_df <- data[, setdiff(colnames(data), response), drop = FALSE]
  } else {
    orig_class_levels <- NULL
    predictors_df <- data
  }

  feature_meta <- extract_feature_metadata(predictors_df)

  structure(
    list(
      fit = model,
      metadata = list(
        features = feature_meta,
        class_levels = orig_class_levels
      ),
      boundary_data = NULL
    ),
    class = "classbound"
  )
}

#' @export
as_classbound.workflow <- function(model, data, response = NULL, ...) {
  if (!requireNamespace("workflows", quietly = TRUE)) {
    stop("Package 'workflows' is required to process workflow objects.", call. = FALSE)
  }
  
  if (!workflows::is_trained_workflow(model)) {
    stop(
      "The provided workflow is not trained. `as_classbound()` requires a pre-fitted model.\n",
      "Please run `fit()` on your workflow first, or use `boundary_workflow_set()` to have classbound train it automatically.",
      call. = FALSE
    )
  }

  # Natively wraps a fitted tidymodels workflow
  NextMethod()
}

#' @export
as_classbound.model_fit <- function(model, data, response = NULL, ...) {
  # Natively wraps a parsnip model_fit
  NextMethod()
}

# Internal helper used by both fit_model() and as_classbound()
extract_feature_metadata <- function(data) {
  imputation_values <- lapply(data, function(col) {
    if (is.numeric(col)) {
      stats::median(col, na.rm = TRUE)
    } else {
      tbl <- table(col)
      if (length(tbl) > 0) {
        mode_val <- names(tbl)[which.max(tbl)]
        if (is.factor(col)) factor(mode_val, levels = levels(col)) else mode_val
      } else {
        if (is.factor(col)) factor(levels(col)[1], levels = levels(col)) else NA
      }
    }
  })

  list(
    names = colnames(data),
    types = lapply(data, class),
    levels = lapply(data, levels),
    imputation = imputation_values
  )
}
