#' Explore high-dimensional classification boundaries
#'
#' @description Generates a high-dimensional uniform domain and identifies points near the decision boundaries of multiple models.
#'
#' @param models A named list of fitted `classbound` model objects.
#' @param data A data frame representing the original training data. The bounding box of this data determines the simulation domain.
#' @param n_points Integer. Number of points to simulate in the high-dimensional space (default: 10000).
#' @param threshold Numeric. A percentile threshold (0 to 1) for determining boundary proximity. Points below this advantage percentile are considered boundary-near. Default: 0.05.
#'
#' @return A data frame containing the simulated coordinates, the predicted classes for each model (prefixed with `Y_`), and a logical `is_boundary` column indicating if the point is near the boundary for ANY of the provided models.
#' @export
boundary_explore_nd <- function(models, data, n_points = 10000, threshold = 0.05) {
  if (!is.list(models) || length(models) == 0) {
    stop("Please provide a named list of one or more classbound models.", call. = FALSE)
  }

  if (is.null(names(models))) {
    names(models) <- paste0("Model_", seq_along(models))
  }

  # Ensure all models are trained on the same features
  first_model <- models[[1]]
  expected_names <- first_model$metadata$features$names

  for (m in models) {
    if (!identical(m$metadata$features$names, expected_names)) {
      stop("All models must be trained on the exact same features.", call. = FALSE)
    }
  }

  # Extract bounding box from data
  if (!all(expected_names %in% colnames(data))) {
    stop("The provided data does not contain all features expected by the models.", call. = FALSE)
  }

  data_features <- data[, expected_names, drop = FALSE]

  # Generate non-aligned grid (uniform random)
  sim_list <- lapply(data_features, function(col) {
    if (is.numeric(col)) {
      stats::runif(n_points, min = min(col, na.rm = TRUE), max = max(col, na.rm = TRUE))
    } else {
      # Sample from categorical levels
      factor(sample(levels(col), n_points, replace = TRUE), levels = levels(col))
    }
  })

  sim_df <- as.data.frame(sim_list)

  # Predict all models
  has_class <- requireNamespace("class", quietly = TRUE)

  is_boundary_global <- rep(FALSE, n_points)

  for (m_name in names(models)) {
    m <- models[[m_name]]
    preds <- predict_model(m, sim_df)

    col_name <- paste0("Y_", m_name)
    sim_df[[col_name]] <- preds$class

    # Calculate advantage
    adv <- rep(1, n_points)

    if (!is.null(preds$probs) && ncol(preds$probs) >= 2) {
      # Model supports probabilities. Advantage = top1 - top2
      prob_mat <- as.matrix(preds$probs)
      adv <- apply(prob_mat, 1, function(x) {
        sx <- sort(x, decreasing = TRUE)
        sx[1] - sx[2]
      })
    } else {
      # No probabilities. Use KNN to approximate advantage
      if (!has_class) {
        warning(paste("Model", m_name, "does not output probabilities and the 'class' package is not installed. Boundary detection will fallback to all points."), call. = FALSE)
        adv <- rep(0, n_points)
      } else {
        # Using KNN: for each point, check its k=10 neighbors. Advantage = proportion of neighbors with same class
        # class::knn internally standardizes or we should? class::knn expects numeric matrices
        numeric_sim <- sim_df[, vapply(sim_df, is.numeric, logical(1)), drop = FALSE]
        if (ncol(numeric_sim) > 0) {
          # KNN is slow for 10000 points. We use a smaller K (e.g. 5)
          knn_res <- class::knn(train = numeric_sim, test = numeric_sim, cl = preds$class, k = 5, prob = TRUE)
          # class::knn returns the winning class and the proportion of votes (prob)
          # The advantage is related to the proportion of votes for the winning class
          # If prob is 1.0 (all 5 neighbors agree), advantage is high. If prob is 0.4, advantage is low.
          adv <- attr(knn_res, "prob")
        }
      }
    }

    # Identify boundary points (bottom X percentile of advantage)
    thresh_val <- stats::quantile(adv, probs = threshold, na.rm = TRUE)
    is_bnd <- adv <= thresh_val
    is_boundary_global <- is_boundary_global | is_bnd
  }

  sim_df$is_boundary <- is_boundary_global

  # Subset to boundary points only
  bnd_df <- sim_df[sim_df$is_boundary, , drop = FALSE]

  return(bnd_df)
}
