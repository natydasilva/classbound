#' Visualize the classification decision boundary
#'
#' @description Renders a 2D plot for boundary exploration.
#'
#' @importFrom rlang .data
#'
#' @details
#' If the `model` object contains a `projection` (e.g. from `tourr` or `prcomp`),
#' `plot_boundary()` will automatically apply forward projection to any provided `obs_data`.
#' It will mathematically map the high-dimensional points onto the 2D plane and calculate
#' their orthogonal distances to the plane. Points exactly on the plane are rendered fully opaque
#' (alpha = 1.0), while points further away gradually fade to greater transparency (alpha = 0.2).
#'
#' Any provided training observations (`obs_data`) are rendered with a white border
#' to improve contrast against the underlying colored decision regions.
#'
#' @param model A `classbound` model object returned by `boundary_compute()`.
#'   Must contain boundary data in `$boundary_data`.
#' @param obs_data Optional data frame of training observations to overlay.
#' @param x_col Column name in `obs_data` for the x-axis feature. Also used as the x-axis label.
#'   Required if `obs_data` is provided.
#' @param y_col Column name in `obs_data` for the y-axis feature. Also used as the y-axis label.
#'   Required if `obs_data` is provided.
#' @param true_label Column name in `obs_data` representing the true class labels. Required if `obs_data` is provided.
#' @param facet_col Optional string naming a column in `boundary` to facet the plot by.
#'   Useful for comparing multiple models (e.g. from `boundary_workflow_set()`).
#' @param type The type of visualization to generate. Supported: '2D' or 'disagreement'.
#' @param show_gradient Logical. If \code{TRUE}, classification probabilities (if available)
#'   will be mapped to the transparency (alpha) of the decision regions, creating a gradient
#'   surface. Defaults to \code{FALSE} (renders a flat decision boundary).
#' @param agree_color Color used for areas where all models agree (only for type='disagreement').
#' @param disagree_color Color used for areas where models disagree (only for type='disagreement').
#' @param obs_alpha Numeric transparency level for overlaid observation points (0.0 to 1.0).
#' @param obs_size Numeric size for overlaid observation points.
#' @param render Character string specifying the rendering method for the decision region. Options are \code{"raster"} (high performance, default) and \code{"tile"} (slower, but fully compatible with interactive graphics like Plotly).
#' @param ... Additional visualization parameters.
#'
#' @return A `ggplot2` object.
#' @examples
#' \donttest{
#' library(palmerpenguins)
#' data(penguins)
#' peng_data <- na.omit(penguins[, c("species", "bill_length_mm", "bill_depth_mm")])
#' m_rpart <- fit_model(peng_data, species ~ ., rpart::rpart)
#' 
#' # Auto-compute bounds and generate 2D grid
#' m_rpart <- boundary_compute(m_rpart, resolution = 50)
#' 
#' # Plot the decision boundary with observations overlaid
#' plot_boundary(m_rpart, obs_data = peng_data, 
#'               x_col = "bill_length_mm", y_col = "bill_depth_mm", 
#'               true_label = "species")
#' }
#' @export
plot_boundary <- function(model, obs_data = NULL, x_col = NULL, y_col = NULL,
                          true_label = NULL, facet_col = NULL, type = "2D",
                          show_gradient = FALSE,
                          agree_color = "#006666", disagree_color = "#FF8000",
                          obs_alpha = 1.0, obs_size = 2.5, render = c("raster", "tile"), ...) {
  render <- match.arg(render)

  if (!type %in% c("2D", "disagreement")) {
    stop("Only type='2D' and type='disagreement' are currently supported.", call. = FALSE)
  }

  if (!inherits(model, "classbound")) {
    stop("model must be a 'classbound' object returned by boundary_compute().", call. = FALSE)
  }

  boundary <- model$boundary_data
  if (is.null(boundary)) {
    stop("model does not contain boundary data. Please run boundary_compute() first.", call. = FALSE)
  }

  req_cols <- c("x", "y", "prediction")
  if (!all(req_cols %in% colnames(boundary))) {
    stop("boundary must contain 'x', 'y', and 'prediction' columns.", call. = FALSE)
  }

  default_x <- if (!is.null(model$boundary_features)) model$boundary_features[1] else "Feature 1"
  default_y <- if (!is.null(model$boundary_features)) model$boundary_features[2] else "Feature 2"

  if (type == "disagreement") {
    if (!inherits(model, "classbound_multi")) {
      stop("Disagreement plots require a multi-model comparison (created by passing a list of models to boundary_compute).", call. = FALSE)
    }

    # Aggregate predictions to find disagreements
    disagree_data <- stats::aggregate(
      prediction ~ x + y,
      data = boundary,
      FUN = function(p) if (length(unique(p)) > 1) "Disagree" else "Agree"
    )
    # Ensure it's a factor for discrete coloring
    disagree_data$prediction <- factor(disagree_data$prediction, levels = c("Agree", "Disagree"))

    p <- ggplot2::ggplot(disagree_data, ggplot2::aes(x = .data$x, y = .data$y))
    if (render == "raster") {
      p <- p + ggplot2::geom_raster(ggplot2::aes(fill = .data$prediction))
    } else {
      p <- p + ggplot2::geom_tile(ggplot2::aes(fill = .data$prediction))
    }
    p <- p +
      ggplot2::scale_fill_manual(
        values = c("Agree" = agree_color, "Disagree" = disagree_color),
        drop = FALSE
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = if (!is.null(x_col)) x_col else default_x,
        y = if (!is.null(y_col)) y_col else default_y,
        fill = "Consensus"
      )
  } else {
    if (show_gradient) {
      # Safely extract probability for the predicted class, leaving NA if the model doesn't support it
      boundary$probability <- NA_real_
      for (lvl in levels(boundary$prediction)) {
        if (lvl %in% colnames(boundary)) {
          idx <- which(boundary$prediction == lvl)
          boundary$probability[idx] <- boundary[[lvl]][idx]
        }
      }
      
      has_any_probs <- any(!is.na(boundary$probability))
      
      if (!has_any_probs) {
        warning(
          "show_gradient = TRUE was requested, but no models contain class probabilities. Falling back to a flat decision boundary.",
          call. = FALSE
        )
      }
    } else {
      has_any_probs <- FALSE
    }

    if (has_any_probs) {
      # Basic plot with alpha mapped to probability (missing probs fallback to 0.3 transparency)
      p <- ggplot2::ggplot(boundary, ggplot2::aes(x = .data$x, y = .data$y))
      if (render == "raster") {
        p <- p + ggplot2::geom_raster(ggplot2::aes(fill = .data$prediction, alpha = .data$probability))
      } else {
        p <- p + ggplot2::geom_tile(ggplot2::aes(fill = .data$prediction, alpha = .data$probability))
      }
      p <- p +
        ggplot2::scale_alpha_continuous(limits = c(0, 1), range = c(0, 1), na.value = 0.3, guide = "none") +
        ggplot2::theme_minimal() +
        ggplot2::scale_fill_discrete(drop = FALSE) +
        ggplot2::labs(
          x = if (!is.null(x_col)) x_col else default_x,
          y = if (!is.null(y_col)) y_col else default_y,
          fill = "Prediction"
        )
    } else {
      # Basic plot with hardcoded alpha (no probabilities)
      p <- ggplot2::ggplot(boundary, ggplot2::aes(x = .data$x, y = .data$y))
      if (render == "raster") {
        p <- p + ggplot2::geom_raster(ggplot2::aes(fill = .data$prediction), alpha = 0.3)
      } else {
        p <- p + ggplot2::geom_tile(ggplot2::aes(fill = .data$prediction), alpha = 0.3)
      }
      p <- p +
        ggplot2::theme_minimal() +
        ggplot2::scale_fill_discrete(drop = FALSE) +
        ggplot2::labs(
          x = if (!is.null(x_col)) x_col else default_x,
          y = if (!is.null(y_col)) y_col else default_y,
          fill = "Prediction"
        )
    }
  }

  # Overlay training observations if provided
  if (!is.null(obs_data)) {
    if (!inherits(obs_data, "data.frame")) {
      stop("obs_data must be a data.frame.", call. = FALSE)
    }

    if (is.null(true_label)) {
      stop("When providing obs_data, you must specify true_label.", call. = FALSE)
    }

    if (!true_label %in% colnames(obs_data)) {
      stop(sprintf("obs_data must contain column '%s'.", true_label), call. = FALSE)
    }

    # Projection workflow
    if (!is.null(model$projection)) {
      proj <- model$projection
      features <- rownames(proj$basis)

      if (!all(features %in% colnames(obs_data))) {
        stop("obs_data is missing features required for the projection basis.", call. = FALSE)
      }

      # Extract numeric matrix for projection
      x_mat <- as.matrix(obs_data[, features])

      # Apply standardization
      if (!is.null(proj$scale)) {
        x_mat <- sweep(x_mat, 2, proj$scale, "/")
      }
      if (!is.null(proj$center)) {
        x_mat <- sweep(x_mat, 2, proj$center, "-")
      }

      # Forward project to 2D coordinates
      z_mat <- x_mat %*% proj$basis

      obs_df <- data.frame(
        x_val = z_mat[, 1],
        y_val = z_mat[, 2],
        true_class = obs_data[[true_label]]
      )

      # Depth fading (opacity scaling)
      # Reconstruct high-dimensional points from the 2D plane
      x_recon <- z_mat %*% t(proj$basis)

      # Calculate orthogonal distance from original points to the 2D plane
      dists <- sqrt(rowSums((x_mat - x_recon)^2))
      max_dist <- max(dists, na.rm = TRUE)

      # Map distance to opacity (1.0 = on plane, fades to 0.2)
      if (max_dist > 0) {
        obs_df$alpha_val <- obs_alpha * (1 - 0.8 * (dists / max_dist))
      } else {
        obs_df$alpha_val <- obs_alpha
      }

      p <- p +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_point(
          data = obs_df,
          ggplot2::aes(x = .data$x_val, y = .data$y_val, fill = .data$true_class, alpha = .data$alpha_val),
          size = obs_size, shape = 21, color = "white", stroke = 0.5
        ) +
        ggplot2::scale_fill_discrete(drop = FALSE) +
        ggplot2::scale_alpha_identity() +
        ggplot2::labs(fill = "True Class")
    } else {
      # Standard 2D workflow
      if (is.null(x_col) || is.null(y_col)) {
        stop("When providing obs_data without a projection, you must specify x_col and y_col.", call. = FALSE)
      }

      if (!all(c(x_col, y_col) %in% colnames(obs_data))) {
        stop(sprintf("obs_data must contain columns '%s' and '%s'.", x_col, y_col), call. = FALSE)
      }

      obs_df <- obs_data[, c(x_col, y_col, true_label)]
      colnames(obs_df) <- c("x_val", "y_val", "true_class")

      p <- p +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_point(
          data = obs_df,
          ggplot2::aes(x = .data$x_val, y = .data$y_val, fill = .data$true_class),
          size = obs_size, alpha = obs_alpha, shape = 21, color = "white", stroke = 0.5
        ) +
        ggplot2::scale_fill_discrete(drop = FALSE) +
        ggplot2::labs(fill = "True Class")
    }
  }

  if (is.null(facet_col) && inherits(model, "classbound_multi") && type == "2D") {
    warning(
      paste0(
        "This model contains multiple boundaries. ",
        "Automatically setting facet_col = 'model' to prevent overlapping plots."
      ),
      call. = FALSE
    )
    facet_col <- "model"
  }

  if (!is.null(facet_col)) {
    if (!facet_col %in% colnames(boundary)) {
      stop(sprintf("facet_col '%s' not found in boundary data.", facet_col), call. = FALSE)
    }
    p <- p + ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym(facet_col)))
  }

  p
}
