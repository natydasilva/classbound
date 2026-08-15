#' Convert a fitted model into a classbound object
#'
#' @description
#' Wraps a pre-fitted model in a `classbound` object, adding the feature metadata
#' (names, types, ranges, imputation values) required by `boundary_compute()` and
#' `plot_boundary()`. This is the entry point for "bring your own model" (BYO) workflows,
#' including native `tidymodels` workflows and `parsnip` model fits.
#'
#' @details
#' Unlike `fit_model()`, `as_classbound()` does **not** refit the model or call
#' `preprocess_data()`. Feature metadata is extracted solely from the `data` argument,
#' which should be the same data used to fit the model. If the model was fitted on
#' scaled or transformed data, ensure `data` reflects that transformation.
#'
#' `as_classbound()` dispatches on the class of `model`:
#' - **Default**: wraps any fitted model object (e.g., `rpart`, `lda`).
#' - **`workflow`**: requires the workflow to already be trained (via `fit()`). Use
#'   `boundary_workflow_set()` to train and wrap an entire workflow set at once.
#' - **`model_fit`**: wraps a fitted `parsnip` model.
#'
#' @param model A fitted model object. For `tidymodels` objects, must be a trained
#'   `workflow` or `model_fit` — not a model specification.
#' @param data A data frame of the training data. Used only to extract feature metadata;
#'   the data is not passed through `preprocess_data()`. Must contain all features that
#'   the model was trained on.
#' @param response Optional string naming the response column in `data`. If provided,
#'   this column is excluded from feature metadata and its levels are stored as class
#'   levels. If `NULL`, class levels will be `NULL` in the returned object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `classbound` object ready for use with `boundary_compute()` and
#'   `plot_boundary()`.
#'
#' @examples
#' \donttest{
#' library(palmerpenguins)
#' data(penguins)
#' peng_data <- na.omit(penguins[, c("species", "bill_length_mm", "bill_depth_mm")])
#'
#' # Wrap any pre-fitted model
#' raw_model <- rpart::rpart(species ~ ., data = peng_data)
#' cb_model <- as_classbound(raw_model, data = peng_data, response = "species")
#'
#' # Tidymodels: wrap a fitted workflow
#' # (boundary_workflow_set() does this automatically for entire workflow sets)
#' }
#' @seealso [fit_model()], [boundary_workflow_set()], [boundary_compute()]
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

  range_values <- lapply(data, function(col) {
    if (is.numeric(col)) {
      val <- suppressWarnings(range(col, na.rm = TRUE))
      if (any(is.infinite(val))) NULL else val
    } else {
      NULL
    }
  })

  list(
    names = colnames(data),
    types = lapply(data, class),
    levels = lapply(data, levels),
    imputation = imputation_values,
    range = range_values
  )
}
