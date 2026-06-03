#' Classifier Adapter Interface
#'
#' @description This file defines the contract for implementing a classifier adapter.
#' Each adapter must implement its own `fit_*` and `predict_*` functions
#' to map specific classifier package implementations (e.g. randomForest, e1071)
#' into the unified `classbound` structure.
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
#   probs <- predict(model, newdata, type = "prob")
#   list(predictions = preds, probabilities = probs)
# }
