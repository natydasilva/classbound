#' Print a classbound model
#'
#' @description Prints a clean summary of a classbound model object, hiding the internal wrapper
#'   structure. For single models, it delegates to the native model's print method. For multi-model
#'   objects (`classbound_multi`), it prints a comparison summary of the included models.
#'
#' @param x A `classbound` or `classbound_multi` object.
#' @param ... Additional arguments passed to the native model's print method.
#' @return The object, invisibly.
#' @export
print.classbound <- function(x, ...) {
  cat("=== Classbound Model Pipeline ===\n")
  cat("Features: ", paste(x$metadata$features$names, collapse = ", "), "\n")
  cat("Classes:  ", paste(x$metadata$class_levels, collapse = ", "), "\n")

  if (!is.null(x$boundary_data)) {
    cat("Boundary: Computed (", nrow(x$boundary_data), " grid points)\n\n", sep = "")
  } else {
    cat("Boundary: Not yet computed (Run `boundary_compute()`)\n\n")
  }

  cat("-- Native Model --\n")
  cat(utils::capture.output(x$fit, ...), sep = "\\n")
  invisible(x)
}

#' @rdname print.classbound
#' @export
print.classbound_multi <- function(x, ...) {
  cat("=== Classbound Model Pipeline ===\n")
  cat("Features: ", paste(x$metadata$features$names, collapse = ", "), "\n")
  cat("Classes:  ", paste(x$metadata$class_levels, collapse = ", "), "\n")

  if (!is.null(x$boundary_data)) {
    cat("Boundary: Computed (", nrow(x$boundary_data), " grid points)\n\n", sep = "")
  } else {
    cat("Boundary: Not yet computed (Run `boundary_compute()`)\n\n")
  }

  cat("-- Multi-Model Comparison --\n")
  cat(length(x$fits), "models compared:\n")
  model_names <- names(x$fits)
  if (is.null(model_names)) model_names <- paste("Model", seq_along(x$fits))
  cat(paste("  -", model_names), sep = "\n")
  invisible(x)
}

#' Plot a classbound model boundary
#'
#' @description A standard S3 plot method that delegates to `plot_boundary()`.
#'
#' @param x A `classbound` model object.
#' @param ... Additional arguments passed to `plot_boundary()`.
#' @return A `ggplot2` object.
#' @export
plot.classbound <- function(x, ...) {
  plot_boundary(model = x, ...)
}

#' Plot a computed boundary
#'
#' @description An S3 plot method that delegates to `plot_boundary()` for
#' objects that have had boundaries computed via `boundary_compute()`.
#'
#' @param x A `classbound_boundary` object.
#' @param ... Additional arguments passed to `plot_boundary()`.
#' @return A `ggplot2` object.
#' @export
plot.classbound_boundary <- function(x, ...) {
  plot_boundary(model = x, ...)
}

#' Summarize a classbound model
#'
#' @description For single models, delegates the summary function directly to the native model
#'   wrapped inside the classbound object. For multi-model objects (`classbound_multi`), it
#'   returns a summary of the contained models.
#'
#' @param object A `classbound` or `classbound_multi` object.
#' @param ... Additional arguments passed to the native model's summary method.
#' @return The summary of the native model.
#' @export
summary.classbound <- function(object, ...) {
  summary(object$fit, ...)
}

#' @rdname summary.classbound
#' @export
summary.classbound_multi <- function(object, ...) {
  cat("Multi-Model Comparison Object containing", length(object$fits), "models.\n")
  invisible(object$fits)
}
