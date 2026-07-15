#' Compute the classification decision boundary
#'
#' @description Computes the decision boundary for a fitted classifier over a specified feature space.
#'
#' @param model A fitted model object of class `classbound`.
#' @param range A named list defining the feature ranges for a 2D space. Each element should be a
#'   numeric vector of length 2 (e.g., `list(Sepal.Length = c(4, 8), Sepal.Width = c(2, 5))`).
#' @param resolution An integer specifying the grid resolution.
#' @param predfun A custom function to generate predictions for non-standard models.
#' @param projection An optional list defining a high-dimensional projection.
#'   Must contain `basis` (a numeric matrix defining the projection frame, which MUST be orthonormal).
#'   Optionally contains `center` and `scale` (numeric vectors) to reverse visual normalization.
#'   Projection is only supported for models trained exclusively on numeric features.
#'   If provided, `range` defines the limits of the 2D projected space, and the grid
#'   is inversely projected back to the original feature space before prediction.
#' @param ... Additional computation parameters.
#'
#' @return A modified `classbound` model object containing the 2D grid points and predictions in `$boundary_data`.
#' @examples
#' \donttest{
#' data(data69_1)
#' data69_1$Y <- as.factor(data69_1$Y)
#' model <- fit_model(data69_1, Y ~ V1 + V2, rpart::rpart)
#' model <- boundary_compute(model, list(V1 = c(-1, 1), V2 = c(-1, 1)))
#' }
#' @export
boundary_compute <- function(model, range, resolution = 100, predfun = NULL, projection = NULL, ...) {
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

  if (!is.list(range) || length(range) != 2 || is.null(names(range))) {
    stop("range must be a named list of length 2 specifying feature ranges.", call. = FALSE)
  }

  var_names <- names(range)

  if (any(duplicated(var_names))) {
    stop("Duplicate feature names found in `range`.", call. = FALSE)
  }

  if (is.null(first_model$metadata$features)) {
    stop("Model metadata is missing. Please retrain the model.", call. = FALSE)
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
    if (length(expected_names) > 2) {
      stop(
        paste0(
          "Visualizing models with >2 features without a projection requires fixed-value slicing, ",
          "which is currently unsupported. Please supply a projection."
        ),
        call. = FALSE
      )
    }
    if (length(var_names) != 2) {
      stop(
        sprintf("Expected exactly 2 boundary variables, but got %d.", length(var_names)),
        call. = FALSE
      )
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
    if (length(missing_names) > 0 || length(invalid_names) > 0) {
      msg <- "Names in `range` do not match the training features."
      if (length(missing_names) > 0) {
        msg <- paste0(msg, "\nMissing features: ", paste(missing_names, collapse = ", "))
      }
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
        projection = projection
      ),
      class = "classbound"
    )
  } else {
    first_model$boundary_data <- final_boundary
    first_model$projection <- projection
    first_model
  }
}
