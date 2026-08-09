#' Compute classification boundaries for a workflow_set
#'
#' @description A dedicated helper to automatically fit and extract classification boundaries
#'   for an entire `workflow_set`. This avoids the need for manual iteration over models.
#'
#' @param wf_set A `workflow_set` object from the `workflowsets` package.
#' @param data A data frame containing the training data. This is required to extract feature
#'   metadata and to fit any workflows that are not yet trained.
#' @param feature_range An optional named list specifying the minimum and maximum values for each feature,
#'   or a character vector of feature names. If `NULL`, the ranges are automatically computed from the training data
#'   (if there are exactly 2 numeric features).
#' @param response A string specifying the name of the response column in `data`.
#' @param resolution An integer specifying the number of points along each axis (default = 100).
#' @param ... Additional arguments passed to `boundary_compute()`.
#'
#' @return A data frame containing the combined boundary grid for all models, with a `model` column
#'   indicating the `wflow_id`.
#' @examples
#' \donttest{
#' library(palmerpenguins)
#' library(workflowsets)
#' library(parsnip)
#' 
#' data(penguins)
#' peng_data <- na.omit(penguins[, c("species", "bill_length_mm", "bill_depth_mm")])
#' 
#' # Define multiple engines
#' spec_rpart <- decision_tree() |> set_engine("rpart") |> set_mode("classification")
#' spec_glm   <- multinom_reg() |> set_engine("nnet") |> set_mode("classification")
#' 
#' # Create a workflow set
#' wf_set <- workflow_set(
#'   preproc = list(base = species ~ bill_length_mm + bill_depth_mm),
#'   models = list(tree = spec_rpart, log_reg = spec_glm)
#' )
#' 
#' # Compute 2D boundaries for all models simultaneously (auto-range)
#' bounds <- boundary_workflow_set(wf_set, peng_data, response = "species", resolution = 30)
#' }
#' @export
boundary_workflow_set <- function(wf_set, data, feature_range = NULL, response, resolution = 100, ...) {
  if (!requireNamespace("workflowsets", quietly = TRUE)) {
    stop("Package 'workflowsets' is required to process workflow_set objects.", call. = FALSE)
  }
  if (!requireNamespace("workflows", quietly = TRUE)) {
    stop("Package 'workflows' is required to check workflow trained status.", call. = FALSE)
  }
  if (!requireNamespace("parsnip", quietly = TRUE)) {
    stop("Package 'parsnip' is required to fit workflows.", call. = FALSE)
  }

  if (missing(data)) {
    stop("`data` must be provided to extract metadata and fit workflows.", call. = FALSE)
  }

  if (missing(response)) {
    stop("`response` must be provided to guarantee that class levels are perfectly synchronized across all models in the set.", call. = FALSE)
  }

  if (!response %in% colnames(data)) {
    stop(sprintf("Response column '%s' not found in `data`.", response), call. = FALSE)
  }

  model_list <- lapply(wf_set$wflow_id, function(id) {
    wf <- workflowsets::extract_workflow(wf_set, id)

    if (!workflows::is_trained_workflow(wf)) {
      wf <- parsnip::fit(wf, data = data)
    }

    as_classbound(wf, data = data, response = response)
  })

  names(model_list) <- wf_set$wflow_id

  boundary_compute(model_list, feature_range = feature_range, resolution = resolution, ...)
}
