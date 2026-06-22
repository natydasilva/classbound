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
#' @section Data Validation Contract:
#' Prior to calling an adapter's `fit_*` method, the core pipeline (via \code{preprocess_data()})
#' guarantees that the inputs are fully validated. Specifically:
#' \itemize{
#'   \item \code{data} is a valid data frame with no missing values.
#'   \item \code{labels} is a valid factor vector with no missing values and no unused levels.
#'   \item All categorical predictors (character columns) in \code{data} have been converted to factors with unused levels dropped.
#' }
#' Because of this strict contract, individual adapter implementations should **assume validated inputs**
#' and should **not duplicate basic validation logic** (e.g., checking for NAs, coercing characters to factors, or dropping levels).
#' Adapters should focus solely on the model-specific logic required to fit and predict.
#'
#' @section Implementation Convention:
#' Adapters that use the formula interface combine features and labels via a reserved
#' column name \code{.label.}: \code{train_data$.label. <- labels}. The core pipeline
#' validates that user data does not contain a column with this reserved name.
#'
#' Example S3 adapter skeleton for a custom classifier "mymodel":
#' ```R
#' fit_adapter.classbound_mymodel <- function(object, data, labels, ...) {
#'   train_data <- data
#'   train_data$.label. <- labels
#'   mypackage::myfit(.label. ~ ., data = train_data, ...)
#' }
#' 
#' predict_model.classbound_mymodel <- function(model, newdata, ...) {
#'   raw_model <- model$model
#'   preds <- as.factor(predict(raw_model, newdata, ...))
#' 
#'   # probs must be either a probability matrix or NULL if unsupported
#'   probs <- NULL
#' 
#'   list(class = preds, probs = probs)
#' }
#' ```
#'
#' @name adapter_interface
NULL
