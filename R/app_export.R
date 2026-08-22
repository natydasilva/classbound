#' Export data to CSV
#'
#' @param data A data frame to export.
#' @param file The file path to write to.
#' @return The file path (invisibly).
#' @export
export_data_csv <- function(data, file) {
  utils::write.csv(data, file = file, row.names = FALSE)
  invisible(file)
}

#' Export fitted models to RDS
#'
#' @param models A list of fitted model objects.
#' @param file The file path to write to.
#' @return The file path (invisibly).
#' @export
export_models_rds <- function(models, file) {
  # Save models as a named list
  saveRDS(models, file = file)
  invisible(file)
}

#' Export comparison metrics to CSV
#'
#' @param metrics A data frame of metrics.
#' @param file The file path to write to.
#' @return The file path (invisibly).
#' @export
export_metrics_csv <- function(metrics, file) {
  utils::write.csv(metrics, file = file, row.names = FALSE)
  invisible(file)
}

#' Export ggplot objects to a multi-page PDF
#'
#' @param plots A list of ggplot objects.
#' @param file The file path to write to.
#' @param width Width of the PDF in inches.
#' @param height Height of the PDF in inches.
#' @return The file path (invisibly).
#' @export
export_plots_pdf <- function(plots, file, width = 8, height = 8) {
  if (requireNamespace("grDevices", quietly = TRUE)) {
    grDevices::pdf(file = file, width = width, height = height)
    on.exit(grDevices::dev.off())
    for (p in plots) {
      suppressWarnings(print(p))
    }
  } else {
    stop("Package 'grDevices' is required for PDF export.", call. = FALSE)
  }
  invisible(file)
}

#' Export ggplot objects to individual PNG files in a directory
#'
#' @param plots A list of ggplot objects (named or unnamed).
#' @param dir The directory path to write the PNGs to.
#' @param dpi Resolution of the PNG files.
#' @param width Width of each PNG in inches.
#' @param height Height of each PNG in inches.
#' @return The directory path (invisibly).
#' @export
export_plots_png <- function(plots, dir, dpi = 300, width = 8, height = 8) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    plot_names <- names(plots)
    if (is.null(plot_names)) {
      plot_names <- paste0("plot_", seq_along(plots))
    }

    for (i in seq_along(plots)) {
      # Sanitize file name
      safe_name <- gsub("[^A-Za-z0-9_.-]", "_", plot_names[i])
      file_path <- file.path(dir, paste0(safe_name, ".png"))
      suppressWarnings(
        ggplot2::ggsave(
          filename = file_path, plot = plots[[i]],
          width = width, height = height, dpi = as.numeric(dpi)
        )
      )
    }
  } else {
    stop("Package 'ggplot2' is required for PNG export.", call. = FALSE)
  }
  invisible(dir)
}

#' Export grid predictions to CSV
#'
#' Saves the boundary grid data (X, Y, predicted class, probabilities) for
#' a fitted classbound model so results can be reproduced outside the Shiny app.
#'
#' @param cb_mod A classbound model object.
#' @param data The training data (used to compute axis ranges).
#' @param resolution Grid resolution (points per axis).
#' @param proj_matrix Optional projection matrix for high-dimensional data.
#' @param proj_info Optional list with center/scale for projection.
#' @param predict_args Additional arguments passed to predict.
#' @param file The file path to write the CSV.
#' @return The file path (invisibly).
#' @keywords internal
export_grid_csv <- function(cb_mod, data, resolution, proj_matrix = NULL, proj_info = NULL, predict_args = list(), slice_x = NULL, slice_y = NULL, file) {
  tryCatch(
    {
      if (!is.null(proj_matrix)) {
        proj_list <- list(basis = proj_matrix, center = proj_info$center, scale = proj_info$scale)
        x_mat <- as.matrix(data[, rownames(proj_matrix)])
        if (!is.null(proj_info$center)) x_mat <- sweep(x_mat, 2, proj_info$center, "-")
        if (!is.null(proj_info$scale)) x_mat <- sweep(x_mat, 2, proj_info$scale, "/")
        z_mat <- x_mat %*% proj_matrix
        range_list <- list(
          PC1 = range(z_mat[, 1]) + c(-0.5, 0.5),
          PC2 = range(z_mat[, 2]) + c(-0.5, 0.5)
        )
        cb_bound <- boundary_compute(cb_mod, feature_range = range_list, resolution = resolution, projection = proj_list, predict_args = predict_args)
      } else {
        feat_cols <- setdiff(colnames(data), "Sim")
        x_name <- if (!is.null(slice_x) && slice_x %in% feat_cols) slice_x else feat_cols[1]
        y_name <- if (!is.null(slice_y) && slice_y %in% feat_cols) slice_y else feat_cols[2]
        range_list <- list()
        range_list[[x_name]] <- range(data[[x_name]]) + c(-0.5, 0.5)
        range_list[[y_name]] <- range(data[[y_name]]) + c(-0.5, 0.5)
        cb_bound <- boundary_compute(cb_mod, feature_range = range_list, resolution = resolution, predict_args = predict_args)
      }
      utils::write.csv(cb_bound$boundary_data, file = file, row.names = FALSE)
      invisible(file)
    },
    error = function(e) {
      warning(sprintf("Grid export failed for model: %s", conditionMessage(e)), call. = FALSE)
      invisible(NULL)
    }
  )
}

#' Export current app configuration to JSON
#'
#' Captures the current state of all Shiny UI inputs as a JSON file for
#' reproducibility. Includes data mode, simulation parameters, grid resolution,
#' and selected classifiers.
#'
#' @param input_list A named list of input values.
#' @param file The file path to write the JSON.
#' @return The file path (invisibly).
#' @keywords internal
#' @importFrom utils str
export_config_json <- function(input_list, file) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    # If jsonlite isn't available, fallback to str()
    utils::capture.output(utils::str(input_list), file = file)
    return(invisible(file))
  }
  jsonlite::write_json(input_list, path = file, pretty = TRUE, auto_unbox = TRUE)
  invisible(file)
}

#' Export a reproducibility R script
#'
#' Writes a self-contained R script that loads the exported data and models
#' and regenerates the boundary plots. Handles both 2D and high-dimensional
#' (projection-based) datasets.
#'
#' @param model_names Character vector of model names included in the export.
#' @param has_projection Logical, whether a projection.rds file was exported.
#' @param file The file path to write the script.
#' @return The file path (invisibly).
#' @keywords internal
export_reproduce_script <- function(model_names, has_projection = FALSE, slice_x = NULL, slice_y = NULL, resolution = 100, show_probs = FALSE, zoom_x = NULL, zoom_y = NULL, file) {
  header <- c(
    "# classbound - Reproduce Exported Results",
    "# This script was auto-generated by classbound's Export Wizard.",
    "# It loads the exported data and models from this ZIP and regenerates the plots.",
    "#",
    "# Note: Ensure this script and the data/model files remain in the same folder.",
    "#",
    "# To run this script in RStudio:",
    "# 1. Open this file in RStudio.",
    "# 2. Go to Session -> Set Working Directory -> To Source File Location.",
    "# 3. Click the 'Source' button in the top right corner of the editor.",
    "",
    "library(classbound)",
    "",
    "# Load Data",
    "if (file.exists(\"data.csv\")) {",
    "  dat <- read.csv(\"data.csv\", stringsAsFactors = TRUE)",
    "  cat(\"Loaded data:\", nrow(dat), \"observations,\", ncol(dat) - 1, \"features\\n\")",
    "} else {",
    "  stop(\"data.csv not found. Make sure 'Data' was included in the export.\")",
    "}",
    "",
    "# Load Models",
    "if (file.exists(\"models.rds\")) {",
    "  models <- readRDS(\"models.rds\")",
    "  cat(\"Loaded models:\", paste(names(models), collapse = \", \"), \"\\n\")",
    "} else {",
    "  stop(\"models.rds not found. Make sure 'Fitted Models' was included in the export.\")",
    "}",
    ""
  )

  zoom_str <- if (!is.null(zoom_x) && !is.null(zoom_y)) {
    sprintf(" + ggplot2::coord_cartesian(xlim = c(%f, %f), ylim = c(%f, %f), expand = FALSE)", zoom_x[1], zoom_x[2], zoom_y[1], zoom_y[2])
  } else {
    ""
  }

  if (has_projection) {
    plot_block <- c(
      "# Load Projection (high-dimensional data was projected to 2D)",
      "if (file.exists(\"projection.rds\")) {",
      "  proj <- readRDS(\"projection.rds\")",
      "  cat(\"Loaded projection matrix:\", nrow(proj$basis), \"features -> 2D\\n\")",
      "} else {",
      "  stop(\"projection.rds not found. This export used a projection for high-dimensional data.\")",
      "}",
      "",
      "# Regenerate Boundary Plots",
      "x_mat <- as.matrix(dat[, rownames(proj$basis)])",
      "if (!is.null(proj$center)) x_mat <- sweep(x_mat, 2, proj$center, \"-\")",
      "if (!is.null(proj$scale)) x_mat <- sweep(x_mat, 2, proj$scale, \"/\")",
      "z_mat <- x_mat %*% proj$basis",
      "range_list <- list(",
      "  PC1 = range(z_mat[, 1]) + c(-0.5, 0.5),",
      "  PC2 = range(z_mat[, 2]) + c(-0.5, 0.5)",
      ")",
      "",
      "for (model_name in names(models)) {",
      "  cat(\"\\nPlotting boundary for:\", model_name, \"\\n\")",
      "  cb_mod <- models[[model_name]]",
      sprintf("  cb_bound <- boundary_compute(cb_mod, feature_range = range_list, resolution = %d, projection = proj)", resolution),
      sprintf("  p <- plot_boundary(cb_bound, obs_data = dat, x_col = \"PC1\", y_col = \"PC2\", true_label = \"Sim\", show_gradient = %s) +", ifelse(show_probs, "TRUE", "FALSE")),
      "    ggplot2::ggtitle(model_name)",
      if (zoom_str != "") sprintf("  p <- suppressMessages(p%s)", zoom_str) else "",
      "  print(p)",
      "}"
    )
  } else {
    slice_x_str <- if (!is.null(slice_x)) sprintf("\"%s\"", slice_x) else "feat_cols[1]"
    slice_y_str <- if (!is.null(slice_y)) sprintf("\"%s\"", slice_y) else "feat_cols[2]"

    plot_block <- c(
      "# Regenerate Boundary Plots",
      "feat_cols <- setdiff(colnames(dat), \"Sim\")",
      sprintf("x_name <- if (%s %%in%% feat_cols) %s else feat_cols[1]", slice_x_str, slice_x_str),
      sprintf("y_name <- if (%s %%in%% feat_cols) %s else feat_cols[2]", slice_y_str, slice_y_str),
      "range_list <- list()",
      "range_list[[x_name]] <- range(dat[[x_name]]) + c(-0.5, 0.5)",
      "range_list[[y_name]] <- range(dat[[y_name]]) + c(-0.5, 0.5)",
      "",
      "for (model_name in names(models)) {",
      "  cat(\"\\nPlotting boundary for:\", model_name, \"\\n\")",
      "  cb_mod <- models[[model_name]]",
      sprintf("  cb_bound <- boundary_compute(cb_mod, feature_range = range_list, resolution = %d)", resolution),
      sprintf("  p <- plot_boundary(cb_bound, obs_data = dat, x_col = x_name, y_col = y_name, true_label = \"Sim\", show_gradient = %s) +", ifelse(show_probs, "TRUE", "FALSE")),
      "    ggplot2::ggtitle(model_name)",
      if (zoom_str != "") sprintf("  p <- suppressMessages(p%s)", zoom_str) else "",
      "  print(p)",
      "}"
    )
  }

  footer <- c(
    "",
    "# Load Metrics (if available)",
    "if (file.exists(\"metrics.csv\")) {",
    "  metrics <- read.csv(\"metrics.csv\")",
    "  cat(\"\\nPerformance Metrics:\\n\")",
    "  print(metrics)",
    "}",
    ""
  )

  writeLines(c(header, plot_block, footer), con = file)
  invisible(file)
}
