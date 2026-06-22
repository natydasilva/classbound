#' Classifier Adapter Interface
#'
#' @description This file defines the contract for implementing a classifier adapter.
#' Each adapter must implement its own `fit_*` and `predict_*` functions
#' to map specific classifier package implementations (e.g. randomForest, e1071)
#' into the unified `classbound` structure.
#'
#' @details
#' The `predict_*` function must return a list with two elements:
#' \itemize{
#'   \item \code{class}: A factor vector of predicted class labels.
#'   \item \code{probs}: A probability matrix, or strictly \code{NULL} if the classifier does not support probabilities.
#'   Downstream pipeline functions (e.g., \code{boundary_compute()}) will handle \code{probs = NULL} natively.
#' }
#'
#' @name adapter_interface
NULL

# Example adapter skeleton:
#
# fit_rf_adapter <- function(data, labels, ...) {
#   randomForest::randomForest(x = data, y = labels, ...)
# }
#
# predict_rf_adapter <- function(model, newdata, ...) {
#   preds <- predict(model, newdata, type = "response")
#
#   # probs must be either a probability matrix or NULL if unsupported
#   probs <- predict(model, newdata, type = "prob")
#
#   list(class = preds, probs = probs)
# }
