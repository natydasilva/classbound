#' classbound: Visualization for Classification Decision Boundaries
#'
#' @description
#' The `classbound` package provides tools for exploring, visualizing, and comparing
#' classification decision boundaries in R. It supports both two-dimensional data and
#' high-dimensional data (via 2D slicing or linear projections), and works with native
#' R classifiers, `tidymodels` workflows, and user-supplied models.
#'
#' @details
#' ## Two main workflows
#'
#' **Interactive workflow**: Launch `explorapp()` to start the built-in Shiny application.
#' From there, you can import your own data, simulate datasets, draw data by hand, choose
#' classifiers, adjust parameters, compare decision boundaries side-by-side, inspect
#' probability surfaces, inject outliers, and export the results.
#'
#' **Programmatic workflow**: Use the modular API directly:
#'
#' ```r
#' model <- fit_model(data, formula, classifier)
#' model <- boundary_compute(model, feature_range, resolution = 100)
#' plot_boundary(model, obs_data = data, x_col = "x", y_col = "y", true_label = "class")
#' ```
#'
#' For a one-step wrapper, use `classbound()`.
#'
#' ## High-dimensional data
#'
#' When a model is trained on more than two features, `boundary_compute()` supports
#' two visualization strategies:
#'
#' - **2D Slice**: two features are selected for the axes; all other numeric features
#'   are fixed at their median and categorical features at their mode.
#' - **Projection**: a projection matrix maps the high-dimensional feature space to two
#'   dimensions (e.g., PCA or a tour basis from the `tourr` package). The boundary grid
#'   is generated in projection space and inverse-projected back for prediction.
#'
#' ## Supported classifiers
#'
#' Any classifier whose `predict()` method returns a vector or factor of class labels
#' works automatically. Built-in adapters are provided for `rpart`, `randomForest`,
#' `PPtreeViz`, `PPtreeExt`, and `ppforest2`. For classifiers that return complex
#' objects (such as lists), use the `predfun` argument to extract class labels.
#' Native `tidymodels` integration is available via `boundary_workflow_set()`.
#'
#' @seealso
#' - `classbound()` for the all-in-one wrapper
#' - `fit_model()` to fit a model
#' - `boundary_compute()` to compute a decision boundary
#' - `plot_boundary()` to visualize a decision boundary
#' - `explorapp()` for the interactive Shiny application
#' - `boundary_workflow_set()` for tidymodels multi-model comparison
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
