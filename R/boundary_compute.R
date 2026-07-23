#' Compute 2D the classification decision boundary 
#'
#' @description Computes the decision boundary for a fitted classifier over a specified feature space in 2D.
#'
#' @param model A fitted model object of class `classbound`.
#' @param range An optional named list defining the feature ranges (e.g., `list(V1 = c(4, 8), V2 = c(2, 5))`),
#'   or a character vector of two feature names to automatically compute bounds from the training data.
#'   If `NULL` and the model has exactly two numeric features, it will automatically use their ranges.
#' @param resolution An integer specifying the grid resolution.
#' @param predfun A custom function to generate predictions for non-standard models.
#'   The function must accept at least two arguments: \code{model} (the fitted native model) 
#'   and \code{newdata} (a data frame of new observations). It should return either a vector/factor 
#'   of predicted classes, or a list containing \code{class} (predicted labels) and \code{probs} (a probability matrix).
#' @param projection An optional list defining a high-dimensional projection.
#'   Must contain `basis` (a numeric matrix defining the projection frame, which MUST be orthonormal).
#'   Optionally contains `center` and `scale` (numeric vectors) to reverse visual normalization.
#'   Projection is only supported for models trained exclusively on numeric features.
#'   If provided, `range` defines the limits of the 2D projected space, and the grid
#'   is inversely projected back to the original feature space before prediction.
#' @param reference An optional named list of fixed values to use for features not specified in `range`.
#'   If NULL, missing numeric features are imputed with their median, and missing categorical features with their mode.
#' @param ... Additional computation parameters.
#'
#' @return A modified `classbound` model object containing the 2D grid points and predictions in `$boundary_data`.
#' @examples
#' \donttest{
#' data(data69_1)
#' data69_1$Y <- as.factor(data69_1$Y)
#' model <- fit_model(data69_1, Y ~ V1 + V2, rpart::rpart)
#' # Auto-compute ranges for the two numeric features
#' model <- boundary_compute(model)
#' 
#' # Palmer penguins example with projection and resolution
#' library(palmerpenguins)
#' data(penguins)
#' penguins <- na.omit(
#'   penguins[, c("species", "bill_length_mm", "bill_depth_mm", "flipper_length_mm")]
#' )
#' m_peng <- fit_model(penguins, species ~ ., rpart::rpart)
#' 
#' # Create a simple orthonormal projection basis (e.g., PCA or manual)
#' basis <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3)[, 1:2]
#' rownames(basis) <- c("bill_length_mm", "bill_depth_mm", "flipper_length_mm")
#' 
#' m_peng_proj <- boundary_compute(
#'   m_peng, 
#'   range = list(PC1 = c(30, 60), PC2 = c(10, 25)), 
#'   resolution = 50,
#'   projection = list(basis = basis)
#' )
#' }
#' @export
boundary_compute <- function(model, range = NULL, resolution = 100, predfun = NULL, projection = NULL, reference = NULL, ...) {
  is_multi <- FALSE
  is_bare_list <- is.list(model) && (class(model)[1] == "list")

  if (is_bare_list) {
    if (!all(vapply(model, inherits, "classbound", FUN.VALUE = logical(1)))) {
      stop("If model is a list, all elements must be 'classbound' objects.", call. = FALSE)
    }
    is_multi <- TRUE
    model_list <- model
  } else if (inherits(model, "classbound")) {
    model_list <- list(model)
  } else {
    stop("model must be a 'classbound' object or a list of 'classbound' objects.", call. = FALSE)
  }

  first_model <- model_list[[1]]

  if (is.null(first_model$metadata$features)) {
    stop("Model metadata is missing. Please retrain the model.", call. = FALSE)
  }

  if (is.null(range) || is.character(range)) {
    avail_features <- first_model$metadata$features$names
    avail_numeric <- avail_features[
      vapply(first_model$metadata$features$types, function(x) any(c("numeric", "integer", "double") %in% x), logical(1))
    ]
    
    if (is.character(range)) {
      if (length(range) != 2) stop("If `range` is a character vector, it must specify exactly 2 feature names.", call. = FALSE)
      target_features <- range
    } else {
      if (length(avail_numeric) == 2) {
        target_features <- avail_numeric
      } else {
        stop("Cannot auto-compute `range`: model does not have exactly 2 numeric features. Please provide `range` explicitly.", call. = FALSE)
      }
    }
    
    range <- list()
    for (f in target_features) {
      rng <- first_model$metadata$features$range[[f]]
      if (is.null(rng)) stop(sprintf("Could not find numeric range for feature '%s'.", f), call. = FALSE)
      range[[f]] <- rng
    }
  }

  if (!is.list(range) || length(range) != 2 || is.null(names(range))) {
    stop("range must be a named list of length 2 specifying feature ranges.", call. = FALSE)
  }

  var_names <- names(range)

  if (any(duplicated(var_names))) {
    stop("Duplicate feature names found in `range`.", call. = FALSE)
  }

  expected_names <- first_model$metadata$features$names
  expected_classes <- first_model$metadata$class_levels

  if (is_multi) {
    for (m in model_list) {
      if (!identical(m$metadata$features$names, expected_names)) {
        stop("All models in a multi-model list must be trained on the exact same features (names and order).", call. = FALSE)
      }
      if (!identical(m$metadata$class_levels, expected_classes)) {
        stop("All models in a multi-model list must share the exact same class levels.", call. = FALSE)
      }
      # If we have feature types (e.g. from preprocessing), validate them too
      if (!is.null(first_model$metadata$features$types) && !is.null(m$metadata$features$types)) {
        if (!identical(m$metadata$features$types, first_model$metadata$features$types)) {
          stop("All models in a multi-model list must share identical feature types.", call. = FALSE)
        }
      }
    }
  }

  if (is.null(projection)) {
    if (length(var_names) != 2) {
      stop(
        sprintf("Expected exactly 2 boundary variables, but got %d.", length(var_names)),
        call. = FALSE
      )
    }

    if (!is.null(reference)) {
      if (!is.list(reference) || is.null(names(reference))) {
        stop("`reference` must be a named list of fixed values.", call. = FALSE)
      }
      ref_invalid <- setdiff(names(reference), expected_names)
      if (length(ref_invalid) > 0) {
        stop(paste0("Names in `reference` do not match training features. Invalid features: ", paste(ref_invalid, collapse = ", ")), call. = FALSE)
      }
    }
  } else {
    if (length(var_names) != 2) {
      stop(
        sprintf("Expected exactly 2 boundary variables for the projection plane, but got %d.", length(var_names)),
        call. = FALSE
      )
    }
  }

  missing_names <- setdiff(expected_names, var_names)
  invalid_names <- setdiff(var_names, expected_names)

  if (is.null(projection)) {
    if (length(invalid_names) > 0) {
      msg <- "Names in `range` do not match the training features."
      if (length(invalid_names) > 0) {
        msg <- paste0(msg, "\nInvalid features: ", paste(invalid_names, collapse = ", "))
      }
      stop(msg, call. = FALSE)
    }
  } else {
    # Validate that all features are numeric
    if (!all(first_model$metadata$features$type %in% c("numeric", "integer", "double"))) {
      stop("High-dimensional projection is only supported for models trained exclusively on numeric features.", call. = FALSE)
    }

    # Validate projection object
    if (!is.list(projection) || is.null(projection$basis)) {
      stop("`projection` must be a list containing at least a `basis` matrix.", call. = FALSE)
    }
    basis <- projection$basis
    if (!is.matrix(basis) || !is.numeric(basis) || anyNA(basis)) {
      stop("`projection$basis` must be a numeric matrix without missing values.", call. = FALSE)
    }
    if (nrow(basis) != length(expected_names) || ncol(basis) != 2) {
      stop(
        sprintf("`projection$basis` must be a %d x 2 matrix to match training features.", length(expected_names)),
        call. = FALSE
      )
    }
    if (!isTRUE(all.equal(crossprod(basis), diag(2), check.attributes = FALSE, tolerance = 1e-6))) {
      stop("The projection basis is not orthonormal. Reconstructing the grid requires an orthonormal basis.", call. = FALSE)
    }
    if (!is.null(rownames(basis)) && !identical(rownames(basis), expected_names)) {
      stop(
        "Row names of `projection$basis` must match the exact feature names and ordering of the training data.",
        call. = FALSE
      )
    }
    if (!is.null(projection$center)) {
      if (!is.numeric(projection$center) || length(projection$center) != length(expected_names) || anyNA(projection$center)) {
        stop(
          sprintf("`center` must be a numeric vector of length %d without missing values.", length(expected_names)),
          call. = FALSE
        )
      }
    }
    if (!is.null(projection$scale)) {
      if (!is.numeric(projection$scale) || length(projection$scale) != length(expected_names) || anyNA(projection$scale)) {
        stop(
          sprintf("`scale` must be a numeric vector of length %d without missing values.", length(expected_names)),
          call. = FALSE
        )
      }
    }
  }

  # Ensure requested boundary features are numeric
  if (is.null(projection)) {
    feature_types <- first_model$metadata$features$types[var_names]
  } else {
    feature_types <- first_model$metadata$features$types
  }
  is_numeric <- vapply(feature_types, function(x) any(c("numeric", "integer") %in% x), logical(1))

  if (any(!is_numeric)) {
    stop("Boundary generation requires numeric features. Categorical ranges are not supported.", call. = FALSE)
  }
  # Generate grid points for each dimension
  seq_x <- seq(from = min(range[[1]]), to = max(range[[1]]), length.out = resolution)
  seq_y <- seq(from = min(range[[2]]), to = max(range[[2]]), length.out = resolution)

  # Create grid
  grid_df <- expand.grid(seq_x, seq_y)
  colnames(grid_df) <- var_names

  # Reconstruct high-dimensional feature space if projection is provided
  if (!is.null(projection)) {
    z_matrix <- as.matrix(grid_df)
    x_matrix <- z_matrix %*% t(projection$basis)
    if (!is.null(projection$scale)) {
      x_matrix <- sweep(x_matrix, 2, projection$scale, "*")
    }
    if (!is.null(projection$center)) {
      x_matrix <- sweep(x_matrix, 2, projection$center, "+")
    }
    predict_df <- as.data.frame(x_matrix)
    colnames(predict_df) <- expected_names
  } else {
    predict_df <- grid_df

    missing_feats <- setdiff(expected_names, colnames(predict_df))
    for (feat in missing_feats) {
      if (!is.null(reference) && feat %in% names(reference)) {
        val <- reference[[feat]]
        expected_levels <- first_model$metadata$features$levels[[feat]]
        if (!is.null(expected_levels)) {
          if (!as.character(val) %in% expected_levels) {
            stop(sprintf("Reference value '%s' for feature '%s' is not a valid level.", val, feat), call. = FALSE)
          }
          val <- factor(val, levels = expected_levels)
        }
        predict_df[[feat]] <- val
      } else {
        predict_df[[feat]] <- first_model$metadata$features$imputation[[feat]]
      }
    }

    # Reorder columns to exactly match expected_names
    predict_df <- predict_df[, expected_names, drop = FALSE]
  }

  # Predict over all models
  res_list <- lapply(seq_along(model_list), function(i) {
    m <- model_list[[i]]
    preds <- predict_model(m, newdata = predict_df, predfun = predfun, ...)

    df <- data.frame(
      x = grid_df[[1]],
      y = grid_df[[2]],
      prediction = preds$class
    )

    if (!is.null(preds$probs)) {
      df <- cbind(df, as.data.frame(preds$probs))
    }

    if (is_multi) {
      model_name <- names(model_list)[i]
      if (is.null(model_name) || model_name == "") model_name <- paste0("Model_", i)
      df <- cbind(model = model_name, df)
    }

    df
  })

  final_boundary <- do.call(rbind, res_list)

  if (is_multi) {
    structure(
      list(
        fit = NULL,
        fits = lapply(model_list, function(x) x$fit),
        metadata = list(
          features = first_model$metadata$features,
          class_levels = expected_classes
        ),
        boundary_data = final_boundary,
        projection = projection,
        boundary_features = var_names
      ),
      class = c("classbound_multi", "classbound")
    )
  } else {
    first_model$boundary_data <- final_boundary
    first_model$projection <- projection
    first_model$boundary_features <- var_names
    first_model
  }
}
