#' Fit a machine learning model
#'
#' @description
#' Fits a classification model and wraps it in a `classbound` object, which carries
#' the feature metadata needed by `boundary_compute()` and `plot_boundary()`.
#'
#' @details
#' ## Interface modes
#'
#' `fit_model()` supports three calling conventions via the `interface` argument:
#'
#' - **`"formula"`** (default): passes `formula` and `data` directly to the classifier.
#'   Works for the vast majority of R classifiers (e.g., `rpart::rpart`, `e1071::svm`,
#'   `stats::lda`).
#' - **`"matrix"`**: constructs a predictor matrix `x` and response vector `y` from
#'   the formula, then calls `classifier(x, y, ...)`. Required for classifiers whose
#'   primary interface is matrix-based (e.g., `randomForest::randomForest`).
#' - **`"custom"`**: passes only `fit_args` to the classifier, giving you full control
#'   over the call. Use this for non-standard APIs.
#'
#' ## Tidymodels support
#'
#' If `classifier` is a `parsnip` model specification (`model_spec`), a fitted
#' `model_fit`, or a `tidymodels` `workflow`, `fit_model()` dispatches to the
#' appropriate method automatically. The `interface` argument is not needed for these
#' objects.
#'
#' ## Preprocessing
#'
#' `fit_model()` calls `preprocess_data()` internally to coerce response labels to
#' a factor, handle missing values, and extract feature metadata. Do not call
#' `preprocess_data()` manually before calling `fit_model()`; this will corrupt the
#' stored metadata.
#'
#' @param data A data frame containing the training features and response variable.
#'   All columns referenced in `formula` must be present.
#' @param formula A formula specifying the response and predictors, e.g.,
#'   `species ~ bill_length_mm + bill_depth_mm` or `species ~ .` to use all columns.
#' @param classifier The classification function or model specification to use.
#'   Pass a function (e.g., `rpart::rpart`), a string name (e.g., `"rpart::rpart"`),
#'   or a `parsnip` model spec / fitted `workflow`.
#' @param interface A string specifying how to invoke the classifier: `"formula"`
#'   (default), `"matrix"`, or `"custom"`. See Details.
#' @param fit_args A named list of additional arguments forwarded to the classifier
#'   during fitting (e.g., `list(cp = 0.01)` for `rpart`).
#' @param ... Additional arguments passed to methods.
#'
#' @return A `classbound` object (a list of class `"classbound"`) containing:
#'   - `$fit`: the raw fitted model object returned by the classifier
#'   - `$metadata`: a list with `$features` (names, types, ranges, imputation values)
#'     and `$class_levels` (sorted character vector of class labels)
#'   - `$boundary_data`: `NULL` until `boundary_compute()` is called
#'
#' @examples
#' \donttest{
#' library(palmerpenguins)
#' data(penguins)
#' peng_data <- na.omit(penguins[, c("species", "bill_length_mm", "bill_depth_mm")])
#'
#' # Formula interface (most classifiers)
#' m_rpart <- fit_model(peng_data, species ~ ., rpart::rpart)
#'
#' # Matrix interface (randomForest)
#' m_rf <- fit_model(peng_data, species ~ ., randomForest::randomForest,
#'   interface = "matrix"
#' )
#'
#' # Additional fitting arguments via fit_args
#' m_rpart_cp <- fit_model(peng_data, species ~ ., rpart::rpart,
#'   fit_args = list(control = rpart::rpart.control(cp = 0.001))
#' )
#' }
#' @seealso [boundary_compute()], [predict_model()], [classbound()]
#' @importFrom stats model.frame model.matrix model.response na.pass
#' @export
fit_model <- function(data, formula, classifier, ...) {
  if (missing(classifier)) {
    stop("Please specify a classifier function.", call. = FALSE)
  }
  UseMethod("fit_model", classifier)
}

#' @export
#' @rdname fit_model
fit_model.default <- function(data, formula, classifier, interface = c("formula", "matrix", "custom"), fit_args = list(), ...) {
  if (missing(classifier)) {
    stop("Please specify a classifier function.", call. = FALSE)
  }
  if (missing(formula)) {
    stop("Please specify a formula.", call. = FALSE)
  }

  interface <- match.arg(interface)

  # Extract labels based on formula
  mf <- model.frame(formula, data = data, na.action = na.pass)
  y <- model.response(mf)

  if (is.null(y)) {
    stop("Could not extract response variable from formula.", call. = FALSE)
  }

  # Capture original class levels.
  orig_class_levels <- if (is.factor(y)) levels(y) else sort(unique(as.character(y)))

  # Preprocess data.
  processed <- preprocess_data(data, y)
  data <- processed$data
  y <- processed$labels

  mf <- model.frame(formula, data = data)
  # Dispatch to model interface.
  if (interface == "formula") {
    model_fit <- do.call(classifier, c(list(formula = formula, data = mf), fit_args))
  } else if (interface == "matrix") {
    x_matrix <- model.matrix(formula, data = mf)[, -1, drop = FALSE] # remove intercept
    model_fit <- do.call(classifier, c(list(x = x_matrix, y = y), fit_args))
  } else if (interface == "custom") {
    model_fit <- do.call(classifier, fit_args)
  }

  # Extract feature metadata.
  response_var <- all.vars(formula[[2]])
  predictors_df <- mf[, setdiff(colnames(mf), response_var), drop = FALSE]
  feature_meta <- extract_feature_metadata(predictors_df)

  # Wrap as classbound object.
  structure(
    list(
      fit = model_fit,
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
#' @rdname fit_model
fit_model.function <- function(data, formula, classifier, interface = c("formula", "matrix", "custom"), fit_args = list(), ...) {
  fit_model.default(data, formula, classifier, interface, fit_args, ...)
}

#' @export
#' @rdname fit_model
fit_model.character <- function(data, formula, classifier, interface = c("formula", "matrix", "custom"), fit_args = list(), ...) {
  # Resolve classifier to function.
  fn <- tryCatch(
    eval(str2lang(classifier)),
    error = function(e) {
      match.fun(classifier)
    }
  )

  if (!is.function(fn)) {
    stop(sprintf("Could not resolve '%s' to a valid function.", classifier), call. = FALSE)
  }

  fit_model.default(data, formula, fn, interface, fit_args, ...)
}
NA
