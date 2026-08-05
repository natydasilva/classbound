#' Visualize a classification decision boundary
#'
#' @description A high-level unified wrapper to fit a model, compute its 2D decision boundary,
#' and plot the results in a single step.
#'
#' @param data A data frame containing the full training dataset. The specific variables
#'   used for modeling and plotting are strictly determined by the `formula`.
#' @param formula A formula specifying the response and exactly two predictors for the 2D visualization.
#' @param classifier The classification function to use (e.g., \code{rpart::rpart}, \code{e1071::svm}).
#'   This works with any R package classification algorithm. If the classifier uses a non-standard 
#'   API, you can adapt it via the `predfun` argument.
#' @param interface A string specifying how to invoke the classifier: \code{"formula"},
#'   \code{"matrix"}, or \code{"custom"}.
#' @param fit_args A named list of additional arguments passed to the classifier during fitting.
#' @param predict_args A named list of additional arguments passed to \code{predict()} during boundary computation.
#' @param predfun A custom function to generate predictions for non-standard models.
#' @param resolution An integer specifying the grid resolution for the decision boundary.
#' @param ... Additional arguments passed to \code{\link{plot_boundary}}.
#'
#' @return A \code{ggplot} object visualizing the 2D decision boundary and original observations.
#' @examples
#' \donttest{
#' library(palmerpenguins)
#' data(penguins)
#' peng_data <- na.omit(penguins[, c("species", "bill_length_mm", "bill_depth_mm")])
#' 
#' # Quick 2D boundary visualization for an SVM
#' classbound(peng_data, species ~ bill_length_mm + bill_depth_mm, e1071::svm)
#' }
#' @export
classbound <- function(data,
                       formula,
                       classifier,
                       interface = c("formula", "matrix", "custom"),
                       fit_args = list(),
                       predict_args = list(),
                       predfun = NULL,
                       resolution = 100,
                       ...) {
  interface <- match.arg(interface)

  # 1. Fit the model
  model <- fit_model(
    data = data,
    formula = formula,
    classifier = classifier,
    interface = interface,
    fit_args = fit_args
  )

  # Compute feature ranges automatically from the preprocessed data
  # The model$features object tracked during fit_model has the column names.
  # We extract the predictors by parsing the formula.
  mf <- model.frame(formula, data = data)
  response_var <- all.vars(formula[[2]])
  predictor_vars <- setdiff(colnames(mf), response_var)

  if (length(predictor_vars) != 2) {
    stop("The classbound() wrapper currently supports exactly 2 predictors for simple 2D visualization. For higher-dimensional data (which requires projections or reference values), please use the separated pipeline: fit_model() -> boundary_compute() -> plot_boundary().", call. = FALSE)
  }

  feature_range <- list()
  for (v in predictor_vars) {
    feature_range[[v]] <- range(mf[[v]], na.rm = TRUE)
  }

  # 2. Compute the decision boundary
  model <- boundary_compute(
    model = model,
    range = feature_range,
    resolution = resolution,
    predfun = predfun,
    predict_args = predict_args
  )

  # 3. Plot the result
  plot_boundary(
    model = model,
    obs_data = data,
    x_col = predictor_vars[1],
    y_col = predictor_vars[2],
    true_label = response_var,
    ...
  )
}
