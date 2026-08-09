#' Visualize a classification decision boundary
#'
#' @description A high-level unified wrapper to fit a model, compute its 2D decision boundary,
#' and plot the results in a single step.
#'
#' @param data A data frame containing the full training dataset. The specific variables
#'   used for modeling and plotting are strictly determined by the `formula`.
#' @param formula A formula specifying the response and predictors.
#' @param classifier The classification function to use (e.g., \code{rpart::rpart}, \code{e1071::svm}).
#'   This works with any R package classification algorithm. If the classifier uses a non-standard 
#'   API, you can adapt it via the `predfun` argument.
#' @param interface A string specifying how to invoke the classifier: \code{"formula"},
#'   \code{"matrix"}, or \code{"custom"}.
#' @param projection An optional list (e.g., \code{list(basis=..., center=..., scale=...)}) specifying a 2D projection for high-dimensional data.
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
                       projection = NULL,
                       fit_args = list(),
                       predict_args = list(),
                       predfun = NULL,
                       resolution = 100,
                       ...) {
  interface <- match.arg(interface)

  model <- fit_model(
    data = data,
    formula = formula,
    classifier = classifier,
    interface = interface,
    fit_args = fit_args
  )

  # Extract feature ranges.
  mf <- model.frame(formula, data = data)
  response_var <- all.vars(formula[[2]])
  predictor_vars <- setdiff(colnames(mf), response_var)

  if (is.null(projection)) {
    if (length(predictor_vars) != 2) {
      stop("To visualize boundaries for high-dimensional data (>2 predictors), please provide a 2D projection matrix (e.g., from PCA or tourr) via the `projection` argument.", call. = FALSE)
    }
    feature_range <- list()
    for (v in predictor_vars) {
      feature_range[[v]] <- range(mf[[v]], na.rm = TRUE)
    }
    x_col_plot <- predictor_vars[1]
    y_col_plot <- predictor_vars[2]
  } else {
    pred_vars <- rownames(projection$basis)
    if (is.null(pred_vars)) {
      stop("projection$basis must have rownames corresponding to the feature names.", call. = FALSE)
    }
    x_mat <- as.matrix(data[, pred_vars, drop = FALSE])
    if (!is.null(projection$center)) x_mat <- scale(x_mat, center = projection$center, scale = FALSE)
    if (!is.null(projection$scale)) x_mat <- scale(x_mat, center = FALSE, scale = projection$scale)
    
    proj_data <- x_mat %*% projection$basis
    pnames <- colnames(projection$basis)
    if (is.null(pnames)) pnames <- c("PC1", "PC2")
    
    feature_range <- list()
    feature_range[[pnames[1]]] <- range(proj_data[, 1], na.rm = TRUE)
    feature_range[[pnames[2]]] <- range(proj_data[, 2], na.rm = TRUE)
    
    x_col_plot <- pnames[1]
    y_col_plot <- pnames[2]
  }

  model <- boundary_compute(
    model = model,
    feature_range = feature_range,
    resolution = resolution,
    predfun = predfun,
    projection = projection,
    predict_args = predict_args
  )

  plot_boundary(
    model = model,
    obs_data = data,
    x_col = x_col_plot,
    y_col = y_col_plot,
    true_label = response_var,
    ...
  )
}
