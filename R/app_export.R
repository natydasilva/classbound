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
      print(p)
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
      ggplot2::ggsave(filename = file_path, plot = plots[[i]], 
                      width = width, height = height, dpi = as.numeric(dpi))
    }
  } else {
    stop("Package 'ggplot2' is required for PNG export.", call. = FALSE)
  }
  invisible(dir)
}
