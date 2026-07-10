#' Print a classbound model
#'
#' @description Prints a clean summary of a classbound model object, hiding the internal wrapper
#'   structure, and delegates to the native model's print method.
#'
#' @param x A `classbound` object.
#' @param ... Additional arguments passed to the native model's print method.
#' @return The object, invisibly.
#' @export
print.classbound <- function(x, ...) {
  cat("=== Classbound Model Pipeline ===\n")
  cat("Features: ", paste(x$metadata$features$names, collapse = ", "), "\n")
  cat("Classes:  ", paste(x$metadata$class_levels, collapse = ", "), "\n\n")
  cat("-- Native Model --\n")
  print(x$fit, ...)
  invisible(x)
}

#' Summarize a classbound model
#'
#' @description Delegates the summary function directly to the native model wrapped inside the classbound object.
#'
#' @param object A `classbound` object.
#' @param ... Additional arguments passed to the native model's summary method.
#' @return The summary of the native model.
#' @export
summary.classbound <- function(object, ...) {
  summary(object$fit, ...)
}
