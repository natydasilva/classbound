#' Visualize a classification decision boundary
#'
#' @description A high-level unified wrapper to fit a model, compute its 2D decision boundary,
#' and plot the results in a single step.
#'
#' @param data A data frame containing the training features and response.
#' @param formula A formula specifying the response and predictors.
#' @param classifier The classification function to use (e.g., \code{rpart::rpart}).
#' @param interface A string specifying how to invoke the classifier: \code{"formula"},
#'   \code{"matrix"}, or \code{"custom"}.
#' @param fit_args A named list of additional arguments passed to the classifier during fitting.
#' @param predict_args A named list of additional arguments passed to \code{predict()} during boundary computation.
#' @param predfun A custom function to generate predictions for non-standard models.
#' @param resolution An integer specifying the grid resolution for the decision boundary.
#' @param ... Additional arguments passed to \code{\link{plot_boundary}}.
#'
#' @return A \code{ggplot} object visualizing the decision boundary and original observations.
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
    stop("classbound() currently supports exactly 2 predictors for 2D visualization.", call. = FALSE)
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
