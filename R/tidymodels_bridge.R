#' Internal Registry for Tidymodels Shiny Bridge
#'
#' This registry maps engine strings (from the Shiny UI) to their canonical
#' `parsnip` model specifications.
#'
#' @details
#' **Design Decision: 1:1 Mapping**
#' For the current scope of the Shiny integration, each supported engine string
#' maps intentionally to exactly one canonical model specification. For example,
#' `"nnet"` strictly maps to `parsnip::mlp()`. This explicit 1:1 mapping keeps
#' the Shiny dropdown interface simple and focuses on the most common models
#' used for exploring decision boundaries, rather than exposing every possible
#' model/engine combination simultaneously.
#'
#' @noRd
.tidymodels_registry <- list(
  "rpart" = function() {
    rlang::check_installed("parsnip")
    parsnip::set_engine(parsnip::decision_tree(mode = "classification"), "rpart")
  },
  "randomForest" = function() {
    rlang::check_installed("parsnip")
    parsnip::set_engine(parsnip::rand_forest(mode = "classification"), "randomForest")
  },
  "kernlab" = function() {
    rlang::check_installed("parsnip")
    parsnip::set_engine(parsnip::svm_rbf(mode = "classification"), "kernlab")
  },
  "nnet" = function() {
    rlang::check_installed("parsnip")
    parsnip::set_engine(parsnip::mlp(mode = "classification"), "nnet")
  }
)

#' Tidymodels Shiny Bridge
#'
#' An internal helper that connects Shiny UI text inputs to the `tidymodels`
#' `workflow_set` engine.
#'
#' @param data A data frame containing the features and target.
#' @param response A string representing the name of the target column.
#' @param models A character vector of engine names (e.g., `c("rpart", "nnet")`).
#' @param range Optional list of axis limits.
#' @param resolution Numeric grid resolution.
#'
#' @return A classbound object containing a multi-model grid.
#' @noRd
tidymodels_bridge <- function(data, response, models, range = NULL, resolution = 100) {
  rlang::check_installed(c("parsnip", "workflowsets"))

  if (length(models) == 0) {
    rlang::abort("The `models` argument cannot be empty.")
  }

  if (length(models) != length(unique(models))) {
    rlang::abort("Duplicate models detected. The `models` argument must contain unique engine names.")
  }

  # Validate against registry
  valid_keys <- names(.tidymodels_registry)
  invalid_models <- setdiff(models, valid_keys)
  if (length(invalid_models) > 0) {
    rlang::abort(
      paste0(
        "Unsupported models requested: ", paste(invalid_models, collapse = ", "), ". ",
        "Supported models are: ", paste(valid_keys, collapse = ", "), "."
      )
    )
  }

  # Fetch specs from registry
  specs <- list()
  for (m in models) {
    specs[[m]] <- .tidymodels_registry[[m]]()
  }

  # Create generic formula
  f <- stats::as.formula(paste(response, "~ ."))

  # Construct workflow set
  wf_set <- workflowsets::workflow_set(
    preproc = list(base = f),
    models = specs
  )

  # Execute boundary computation
  boundary_workflow_set(
    wf_set,
    data = data,
    response = response,
    range = range,
  )
}

#' Find Workspace Models
#'
#' Scans an environment for objects that inherit from `workflow`, `model_fit`, or `model_spec`.
#'
#' @param env The environment to scan, defaults to `.GlobalEnv`.
#' @return A character vector of object names.
#' @noRd
find_workspace_models <- function(env = .GlobalEnv) {
  objs <- ls(envir = env)
  if (length(objs) == 0) return(character(0))
  
  is_model <- vapply(objs, function(x) {
    obj <- get(x, envir = env)
    inherits(obj, c("workflow", "model_fit", "model_spec"))
  }, logical(1))
  
  objs[is_model]
}

#' @export
fit_model.workflow <- function(data, formula, classifier, ...) {
  rlang::check_installed(c("parsnip", "workflows"))
  
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  orig_class_levels <- if (is.factor(y)) levels(y) else sort(unique(as.character(y)))
  
  processed <- preprocess_data(data, y)
  data <- processed$data
  
  # A fresh workflow from workflow() %>% add_model() lacks a preprocessor (formula).
  # If it lacks a preprocessor, we add the provided formula.
  has_pre <- tryCatch({ workflows::extract_preprocessor(classifier); TRUE }, error = function(e) FALSE)
  if (!has_pre) {
    classifier <- workflows::add_formula(classifier, formula)
  }
  
  model_fit <- parsnip::fit(classifier, data = data)
  
  response_var <- all.vars(formula[[2]])
  predictors_df <- data[, setdiff(colnames(data), response_var), drop = FALSE]
  feature_meta <- extract_feature_metadata(predictors_df)
  
  structure(
    list(
      fit = model_fit,
      metadata = list(features = feature_meta, class_levels = orig_class_levels),
      boundary_data = NULL
    ),
    class = "classbound"
  )
}

#' @export
fit_model.model_spec <- function(data, formula, classifier, ...) {
  rlang::check_installed("parsnip")
  
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  orig_class_levels <- if (is.factor(y)) levels(y) else sort(unique(as.character(y)))
  
  processed <- preprocess_data(data, y)
  data <- processed$data
  
  model_fit <- parsnip::fit(classifier, formula, data = data)
  
  response_var <- all.vars(formula[[2]])
  predictors_df <- data[, setdiff(colnames(data), response_var), drop = FALSE]
  feature_meta <- extract_feature_metadata(predictors_df)
  
  structure(
    list(
      fit = model_fit,
      metadata = list(features = feature_meta, class_levels = orig_class_levels),
      boundary_data = NULL
    ),
    class = "classbound"
  )
}

#' @export
fit_model.model_fit <- function(data, formula, classifier, ...) {
  # Reconstruct model frame to retrieve true labels
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  orig_class_levels <- if (is.factor(y)) levels(y) else sort(unique(as.character(y)))
  
  response_var <- all.vars(formula[[2]])
  predictors_df <- data[, setdiff(colnames(data), response_var), drop = FALSE]
  feature_meta <- extract_feature_metadata(predictors_df)
  
  structure(
    list(
      fit = classifier,
      metadata = list(features = feature_meta, class_levels = orig_class_levels),
      boundary_data = NULL
    ),
    class = "classbound"
  )
}

#' @export
fit_model.workflow <- function(data, formula, classifier, ...) {
  rlang::check_installed("workflows")
  
  # Workflows embed a fixed preprocessing formula. 
  # Extract the underlying parsnip spec to allow refitting with the internal canvas formula (Sim ~ .).
  spec <- workflows::extract_spec_parsnip(classifier)
  
  fit_model(data, formula, spec, ...)
}
