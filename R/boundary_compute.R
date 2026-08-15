#' Compute the classification decision boundary
#'
#' @description
#' Generates a 2D grid of predicted class labels (and optionally class probabilities)
#' for a fitted `classbound` model. The resulting boundary data is stored in the returned
#' object and consumed directly by `plot_boundary()`.
#'
#' @details
#' ## 2D slice (two-feature visualization)
#'
#' When `projection` is `NULL`, `boundary_compute()` generates a regular grid over the
#' two features named in `feature_range`. If the model was trained on more than two
#' features, all remaining numeric features are fixed at their training-set median and
#' categorical features at their training-set mode — unless you supply explicit values
#' via `reference`.
#'
#' This is a **2D slice** of the full multivariate decision boundary. Two observations
#' that appear in the same region may still be separated in a dimension that is held
#' fixed. Use this mode when you want to examine how two specific features interact,
#' or when your model was trained on exactly two features.
#'
#' ## Projection (high-dimensional visualization)
#'
#' When `projection` is provided, the 2D grid is generated in the projected space and
#' then **inverse-projected** back to the original feature space before prediction.
#' This means the model uses all of its training features; only the visualization is
#' collapsed to two dimensions.
#'
#' The projection `$basis` must be a numeric matrix with:
#' - row count equal to the number of training features, row names matching feature names
#' - exactly two columns (the projected axes)
#' - orthonormal columns (`crossprod(basis)` must equal `diag(2)`)
#'
#' Suitable bases can be obtained from `prcomp()` or the `tourr` package.
#'
#' ## Multi-model comparison
#'
#' Pass a named list of `classbound` objects (all trained on the same features and class
#' levels) to compute boundaries for multiple models at once. The result contains a
#' `model` column, suitable for use with `plot_boundary(facet_col = "model")`.
#'
#' @param model A `classbound` model object returned by `fit_model()` or `as_classbound()`,
#'   or a named list of `classbound` objects for multi-model comparison. All models in a
#'   list must share the same training features and class levels.
#' @param feature_range A named list of length 2 specifying the axis ranges, e.g.,
#'   `list(bill_length_mm = c(30, 60), bill_depth_mm = c(10, 25))`. Alternatively, a
#'   character vector of exactly two feature names (ranges computed from training data).
#'   If `NULL` and the model has exactly two numeric features, ranges are auto-detected.
#'   When `projection` is provided, this defines the limits of the 2D projected space.
#' @param resolution An integer >= 2 specifying the number of grid points per axis.
#'   Higher values produce smoother boundaries. Default is 100; use 50 for faster
#'   interactive exploration.
#' @param predfun An optional custom prediction function for non-standard classifiers.
#'   Must accept `(model, newdata, ...)` and return a vector/factor of predicted classes,
#'   or a list with `$class` (factor) and `$probs` (probability matrix or `NULL`).
#' @param projection An optional named list defining a high-dimensional projection.
#'   Must contain `$basis` (a numeric matrix, rows = features, columns = 2 axes,
#'   must be orthonormal). Optionally contains `$center` and `$scale` (numeric vectors
#'   of length equal to the number of features) to reverse pre-projection standardization.
#'   Only supported for models trained on numeric features.
#' @param reference An optional named list of fixed reference values for features not
#'   specified in `feature_range` (2D slice mode only). Names must match training feature
#'   names. If `NULL`, numeric features are imputed at their median and categorical
#'   features at their mode.
#' @param ... Additional arguments passed to `predict_model()`.
#'
#' @return A modified `classbound` object with the additional class `"classbound_boundary"`.
#'   The boundary grid is stored in `$boundary_data` (columns: `x`, `y`, `prediction`,
#'   and per-class probability columns when available). For multi-model input, also
#'   contains a `model` column.
#'
#' @examples
#' \donttest{
#' library(palmerpenguins)
#' data(penguins)
#' peng_data <- na.omit(penguins[, c("species", "bill_length_mm", "bill_depth_mm")])
#'
#' # Fit and compute 2D boundary with auto-detected ranges
#' m <- fit_model(peng_data, species ~ ., rpart::rpart)
#' m <- boundary_compute(m, resolution = 50)
#' head(m$boundary_data)
#'
#' # Explicit axis ranges
#' m <- boundary_compute(m,
#'   feature_range = list(bill_length_mm = c(30, 60), bill_depth_mm = c(10, 25)),
#'   resolution = 100
#' )
#'
#' # 3-feature model visualized as a 2D slice
#' peng3d <- na.omit(penguins[, c(
#'   "species", "bill_length_mm",
#'   "bill_depth_mm", "flipper_length_mm"
#' )])
#' m3d <- fit_model(peng3d, species ~ ., rpart::rpart)
#' m3d_slice <- boundary_compute(m3d,
#'   feature_range = list(bill_length_mm = c(30, 60), bill_depth_mm = c(10, 25)),
#'   resolution = 50
#' )
#'
#' # 3-feature model visualized via PCA projection
#' feat_cols <- c("bill_length_mm", "bill_depth_mm", "flipper_length_mm")
#' pca <- prcomp(peng3d[, feat_cols], scale. = TRUE)
#' basis <- pca$rotation[, 1:2]
#' m3d_proj <- boundary_compute(m3d,
#'   feature_range = list(PC1 = c(-4, 4), PC2 = c(-3, 3)),
#'   resolution = 50,
#'   projection = list(basis = basis, center = pca$center, scale = pca$scale)
#' )
#' }
#' @seealso [fit_model()], [plot_boundary()], [boundary_workflow_set()]
#' @export
boundary_compute <- function(model, feature_range = NULL, resolution = 100, predfun = NULL, projection = NULL, reference = NULL, ...) {
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

  if (is.null(feature_range) || is.character(feature_range)) {
    avail_features <- first_model$metadata$features$names
    avail_numeric <- avail_features[
      vapply(first_model$metadata$features$types, function(x) any(c("numeric", "integer", "double") %in% x), logical(1))
    ]

    if (is.character(feature_range)) {
      if (length(feature_range) != 2) stop("If `feature_range` is a character vector, it must specify exactly 2 feature names.", call. = FALSE)
      target_features <- feature_range
    } else {
      if (length(avail_numeric) == 2) {
        target_features <- avail_numeric
      } else {
        stop("Cannot auto-compute `feature_range`: model does not have exactly 2 numeric features. Please provide `feature_range` explicitly.", call. = FALSE)
      }
    }

    feature_range <- list()
    for (f in target_features) {
      rng <- first_model$metadata$features$range[[f]]
      if (is.null(rng)) stop(sprintf("Could not find numeric range for feature '%s'.", f), call. = FALSE)
      feature_range[[f]] <- rng
    }
  }

  if (!is.list(feature_range) || length(feature_range) != 2 || is.null(names(feature_range))) {
    stop("feature_range must be a named list of length 2 specifying feature ranges.", call. = FALSE)
  }

  if (!is.numeric(resolution) || resolution < 2) {
    stop("`resolution` must be an integer >= 2 to generate a valid rendering grid.", call. = FALSE)
  }

  var_names <- names(feature_range)

  if (any(duplicated(var_names))) {
    stop("Duplicate feature names found in `feature_range`.", call. = FALSE)
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
      # Validate feature types.
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

  invalid_names <- setdiff(var_names, expected_names)

  if (is.null(projection)) {
    if (length(invalid_names) > 0) {
      msg <- "Names in `feature_range` do not match the training features."
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

  if (is.null(projection)) {
    feature_types <- first_model$metadata$features$types[var_names]
  } else {
    feature_types <- first_model$metadata$features$types
  }
  is_numeric <- vapply(feature_types, function(x) any(c("numeric", "integer") %in% x), logical(1))

  if (any(!is_numeric)) {
    stop("Boundary generation requires numeric features. Categorical ranges are not supported.", call. = FALSE)
  }
  seq_x <- seq(from = min(feature_range[[1]]), to = max(feature_range[[1]]), length.out = resolution)
  seq_y <- seq(from = min(feature_range[[2]]), to = max(feature_range[[2]]), length.out = resolution)

  # Generate grid.
  grid_df <- expand.grid(seq_x, seq_y)
  colnames(grid_df) <- var_names

  # Apply projection matrix.
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
    imputed_silently <- c()
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
        imputed_silently <- c(imputed_silently, feat)
      }
    }

    if (length(imputed_silently) > 0) {
      warning(
        sprintf("The following high-dimensional features were not provided and were automatically imputed with their median/mode: %s", paste(imputed_silently, collapse = ", ")),
        call. = FALSE
      )
    }

    # Restore expected column ordering.
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
  # Pad missing columns with NA.
  all_cols <- unique(unlist(lapply(res_list, colnames)))
  res_list <- lapply(res_list, function(df) {
    missing_cols <- setdiff(all_cols, colnames(df))
    if (length(missing_cols) > 0) {
      df[missing_cols] <- NA
    }
    df[, all_cols, drop = FALSE]
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
      class = c("classbound_boundary", "classbound_multi", "classbound")
    )
  } else {
    first_model$boundary_data <- final_boundary
    first_model$projection <- projection
    first_model$boundary_features <- var_names
    class(first_model) <- c("classbound_boundary", class(first_model))
    first_model
  }
}
