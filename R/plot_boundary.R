#' Visualize the classification decision boundary
#'
#' @description
#' Renders a 2D ggplot of a classification decision boundary, optionally overlaying
#' observed training data. Input must be a `classbound` object that has already been
#' passed through `boundary_compute()`.
#'
#' @details
#' ## Probability surface (gradient)
#'
#' When `show_gradient = TRUE`, decision regions are shaded by the predicted class
#' probability: deep, saturated regions indicate high model confidence, while faded
#' regions indicate uncertainty near the boundary. Probability shading is only possible
#' when the underlying classifier returns class probabilities. Classifiers that return
#' only class labels (e.g., standard SVMs or PPtree models) produce a flat boundary
#' regardless of `show_gradient`.
#'
#' ## High-dimensional projections and depth fading
#'
#' When the `classbound` object contains a projection (from `boundary_compute(..., projection = ...)`),
#' `plot_boundary()` automatically forward-projects any `obs_data` observations onto the
#' 2D plane. It also computes each point's orthogonal distance from the projection plane
#' and maps this distance to opacity: points lying exactly on the plane are fully opaque
#' (`alpha = 1.0`), while points further away in the original feature space gradually
#' fade toward `alpha = 0.2`. This depth-fading provides visual cues about how faithfully
#' each point's position is captured by the current projection.
#'
#' ## Rendering backend and plotly compatibility
#'
#' The default `render = "raster"` uses `ggplot2::geom_raster()`, which is fast and
#' produces high-quality static output. However, `geom_raster()` is not supported by
#' `plotly::ggplotly()` and will produce a blank interactive plot. To convert a boundary
#' plot to an interactive plotly figure, use `render = "tile"` instead:
#'
#' ```r
#' p <- plot_boundary(model, render = "tile")
#' plotly::ggplotly(p)
#' ```
#'
#' ## Colors and palettes
#'
#' By default, `plot_boundary()` uses `classbound_palette()`, a curated 20-color palette
#' with deterministic (alphabetical) class-to-color assignment, ensuring consistent colors
#' across multiple plots. Supply a `palette` name (e.g., `"Dark2"`) to use an RColorBrewer
#' palette, or supply `colors` as a named vector for explicit control. If an RColorBrewer
#' palette cannot accommodate the number of classes, it falls back to `classbound_palette()`.
#'
#' @param model A `classbound` object returned by `boundary_compute()`. Must contain
#'   boundary data in `$boundary_data`.
#' @param obs_data An optional data frame of observations to overlay on the boundary plot.
#'   Typically the training data. If provided, `x_col`, `y_col`, and `true_label` are
#'   required.
#' @param x_col Column name in `obs_data` for the x-axis feature. Used as the x-axis label.
#'   Required when `obs_data` is provided and the model has no projection.
#' @param y_col Column name in `obs_data` for the y-axis feature. Used as the y-axis label.
#'   Required when `obs_data` is provided and the model has no projection.
#' @param true_label Column name in `obs_data` containing the true class labels.
#'   Required when `obs_data` is provided.
#' @param facet_col Optional string naming a column in the boundary data to facet the
#'   plot by. Set to `"model"` automatically when the input is a multi-model comparison
#'   object. Also compatible with any column from `boundary_workflow_set()` output.
#' @param type The visualization type. `"2D"` (default) renders the predicted class
#'   regions. `"disagreement"` renders a binary map showing where models agree vs.
#'   disagree; requires a multi-model input.
#' @param show_gradient Logical. If `TRUE`, decision regions are shaded by the predicted
#'   class probability (requires a classifier that provides probabilities). Defaults to
#'   `FALSE`.
#' @param agree_color Color for regions where all models agree (only for
#'   `type = "disagreement"`).
#' @param disagree_color Color for regions where models disagree (only for
#'   `type = "disagreement"`).
#' @param obs_alpha Numeric transparency for overlaid observation points (0.0–1.0).
#'   When a projection is active, this is treated as the maximum opacity; actual alpha
#'   varies by depth-fading.
#' @param obs_size Numeric point size for overlaid observations.
#' @param render Rendering method for decision regions: `"raster"` (default, fast, not
#'   plotly-compatible) or `"tile"` (slower, compatible with `plotly::ggplotly()`).
#' @param colors Optional named character vector mapping class labels to colors, e.g.
#'   `c("Adelie" = "#E6194B", "Chinstrap" = "#3CB44B", "Gentoo" = "#4363D8")`. Overrides
#'   both `palette` and the default `classbound_palette()`.
#' @param highlight_outliers Logical. If `TRUE`, observations with `is_outlier == TRUE`
#'   in `obs_data` are rendered as diamonds (shape 23) instead of circles.
#' @param palette Optional RColorBrewer palette name (e.g., `"Dark2"`, `"Set1"`) to
#'   override the default colors. Falls back to `classbound_palette()` if the palette
#'   cannot support the number of classes.
#' @param xlim Optional length-2 numeric vector to set x-axis limits.
#' @param ylim Optional length-2 numeric vector to set y-axis limits.
#' @param ... Additional arguments (currently unused).
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' \donttest{
#' library(palmerpenguins)
#' data(penguins)
#' peng_data <- na.omit(penguins[, c("species", "bill_length_mm", "bill_depth_mm")])
#'
#' m <- fit_model(peng_data, species ~ ., rpart::rpart)
#' m <- boundary_compute(m, resolution = 50)
#'
#' # Basic boundary plot with observations
#' plot_boundary(m,
#'   obs_data   = peng_data,
#'   x_col      = "bill_length_mm",
#'   y_col      = "bill_depth_mm",
#'   true_label = "species"
#' )
#'
#' # Probability gradient (rpart supports probabilities)
#' plot_boundary(m,
#'   obs_data      = peng_data,
#'   x_col         = "bill_length_mm",
#'   y_col         = "bill_depth_mm",
#'   true_label    = "species",
#'   show_gradient = TRUE
#' )
#'
#' # Plotly-compatible rendering
#' p <- plot_boundary(m,
#'   obs_data   = peng_data,
#'   x_col      = "bill_length_mm",
#'   y_col      = "bill_depth_mm",
#'   true_label = "species",
#'   render     = "tile"
#' )
#' # plotly::ggplotly(p)  # uncomment to convert to interactive
#' }
#' @seealso [boundary_compute()], [classbound()], [classbound_palette()]
#' @importFrom rlang .data
#' @export
plot_boundary <- function(model, obs_data = NULL, x_col = NULL, y_col = NULL,
                          true_label = NULL, facet_col = NULL, type = "2D",
                          show_gradient = FALSE, agree_color = "#006666", disagree_color = "#FF8000",
                          obs_alpha = 1.0, obs_size = 2.5, render = c("raster", "tile"), colors = NULL,
                          palette = NULL, highlight_outliers = FALSE, xlim = NULL, ylim = NULL, ...) {
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

  levs <- if (is.factor(boundary$prediction)) levels(boundary$prediction) else unique(as.character(boundary$prediction))
  if (type != "disagreement") {
    if (!is.null(colors)) {
      if (is.null(names(colors))) {
        stop("The 'colors' argument must be a named vector (e.g., c('Class 1' = 'red')).", call. = FALSE)
      }
      missing_classes <- setdiff(levs, names(colors))
      if (length(missing_classes) > 0) {
        stop(sprintf("The 'colors' vector is missing colors for the following classes: %s", paste(missing_classes, collapse = ", ")), call. = FALSE)
      }
      unused_classes <- setdiff(names(colors), levs)
      if (length(unused_classes) > 0) {
        warning(sprintf("The 'colors' vector contains unused classes not present in the data: %s", paste(unused_classes, collapse = ", ")), call. = FALSE)
      }
      color_scale <- ggplot2::scale_fill_manual(values = colors, drop = FALSE)
    } else if (!is.null(palette)) {
      if (!palette %in% rownames(RColorBrewer::brewer.pal.info)) {
        stop(sprintf("Invalid palette name: '%s'", palette), call. = FALSE)
      }
      max_colors <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
      if (length(levs) > max_colors) {
        warning(sprintf("Palette '%s' only supports %d classes but data has %d. Falling back to classbound_palette().", palette, max_colors, length(levs)), call. = FALSE)
        color_scale <- ggplot2::scale_fill_manual(values = classbound_palette(levs), drop = FALSE)
      } else {
        color_scale <- ggplot2::scale_fill_brewer(palette = palette, drop = FALSE)
      }
    } else {
      color_scale <- ggplot2::scale_fill_manual(values = classbound_palette(levs), drop = FALSE)
    }
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
      ggplot2::theme(aspect.ratio = 1) +
      ggplot2::labs(
        x = if (!is.null(x_col)) x_col else default_x,
        y = if (!is.null(y_col)) y_col else default_y,
        fill = "Consensus"
      )
  } else {
    if (show_gradient) {
      # Extract probability for the predicted class (NA if unsupported)
      boundary$probability <- NA_real_
      for (lvl in levels(boundary$prediction)) {
        if (lvl %in% colnames(boundary)) {
          idx <- which(boundary$prediction == lvl)
          boundary$probability[idx] <- boundary[[lvl]][idx]
        }
      }

      has_any_probs <- any(!is.na(boundary$probability))
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
        ggplot2::scale_alpha_continuous(limits = c(0, 1), range = c(0.1, 1), na.value = 0.3, guide = "none") +
        ggplot2::theme_minimal() +
        ggplot2::theme(aspect.ratio = 1) +
        color_scale +
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
        ggplot2::theme(aspect.ratio = 1) +
        color_scale +
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
      if (!is.null(proj$center)) {
        x_mat <- sweep(x_mat, 2, proj$center, "-")
      }
      if (!is.null(proj$scale)) {
        x_mat <- sweep(x_mat, 2, proj$scale, "/")
      }

      # Forward project to 2D coordinates
      z_mat <- x_mat %*% proj$basis

      obs_df <- data.frame(
        x_val = z_mat[, 1],
        y_val = z_mat[, 2],
        true_class = obs_data[[true_label]]
      )

      if (highlight_outliers) {
        obs_df$is_outlier <- if ("is_outlier" %in% colnames(obs_data)) obs_data$is_outlier else FALSE
      } else {
        obs_df$is_outlier <- FALSE
      }

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
        ggnewscale::new_scale("alpha") +
        ggplot2::geom_point(
          data = obs_df[!obs_df$is_outlier, , drop = FALSE],
          ggplot2::aes(x = .data$x_val, y = .data$y_val, fill = .data$true_class, alpha = .data$alpha_val),
          size = obs_size, shape = 21, color = "white", stroke = 0.5
        ) +
        ggplot2::geom_point(
          data = obs_df[obs_df$is_outlier, , drop = FALSE],
          ggplot2::aes(x = .data$x_val, y = .data$y_val, fill = .data$true_class, alpha = .data$alpha_val),
          size = obs_size, shape = 23, color = "black", stroke = 1.2, show.legend = FALSE
        ) +
        color_scale +
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

      if (!is.null(model$metadata$features$names) && length(model$metadata$features$names) > 2) {
        warning("Plotting obs_data on a >2D slice projects points flat, which can be visually misleading. Consider using a projection matrix (e.g., PCA) via boundary_compute() instead.", call. = FALSE)
      }

      obs_df <- obs_data[, c(x_col, y_col, true_label)]
      colnames(obs_df) <- c("x_val", "y_val", "true_class")

      if (highlight_outliers) {
        obs_df$is_outlier <- if ("is_outlier" %in% colnames(obs_data)) obs_data$is_outlier else FALSE
      } else {
        obs_df$is_outlier <- FALSE
      }

      p <- p +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_point(
          data = obs_df[!obs_df$is_outlier, , drop = FALSE],
          ggplot2::aes(x = .data$x_val, y = .data$y_val, fill = .data$true_class),
          size = obs_size, alpha = obs_alpha, shape = 21, color = "white", stroke = 0.5
        ) +
        ggplot2::geom_point(
          data = obs_df[obs_df$is_outlier, , drop = FALSE],
          ggplot2::aes(x = .data$x_val, y = .data$y_val, fill = .data$true_class),
          size = obs_size, alpha = 1, shape = 23, color = "black", stroke = 1.2, show.legend = FALSE
        ) +
        color_scale +
        ggplot2::labs(fill = "True Class")
    }
  }

  if (is.null(facet_col) && inherits(model, "classbound_multi") && type == "2D") {
    message("This model contains multiple boundaries. Automatically setting facet_col = 'model' to prevent overlapping plots.")
    facet_col <- "model"
  }

  if (!is.null(facet_col)) {
    if (!facet_col %in% colnames(boundary)) {
      stop(sprintf("facet_col '%s' not found in boundary data.", facet_col), call. = FALSE)
    }
    p <- p + ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym(facet_col)))
  }
  p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE)
  p
}
