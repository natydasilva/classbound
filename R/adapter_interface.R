#' Classifier Adapter Interface
#'
#' @description This file defines the contract for implementing a classifier adapter.
#' Each adapter must implement its own `predict_adapter` function
#' to map specific classifier package implementations (e.g. randomForest, e1071)
#' into the unified `classbound` structure.
#'
#' @details
#' The `predict_adapter` function must return a list with two elements:
#' \itemize{
#'   \item \code{class}: A factor vector of predicted class labels.
#'   \item \code{probs}: A probability matrix, or strictly \code{NULL} if the classifier does not support probabilities.
#'   Downstream pipeline functions (e.g., \code{boundary_compute()}) will handle \code{probs = NULL} natively.
#' }
#'
#' @section Data Validation Contract:
#' Prior to calling a classifier's fitting logic, the core pipeline (via \code{preprocess_data()})
#' guarantees that the inputs are fully validated. Specifically:
#' \itemize{
#'   \item The data is a valid data frame with no missing values.
#'   \item The labels are a valid factor vector with no missing values and no unused levels.
#'   \item All categorical predictors (character columns) in the data have been converted to factors with unused levels dropped.
#' }
#'
#' @section User-Defined Classifiers:
#' Because `classbound` uses standard S3 dispatch for predictions, you can natively support
#' any custom classifier without modifying the package source code. You only need to pass 
#' your classifier function to \code{fit_model()}, and then define a \code{predict_adapter} 
#' S3 method for your model's native class.
#' 
#' Example S3 adapter skeleton for a custom classifier that returns class "mymodel":
#' ```R
#' # 1. Define the prediction logic for the native model class
#' predict_adapter.mymodel <- function(model, newdata, ...) {
#'   # Get predictions from the custom package
#'   raw_preds <- mypackage::predict(model, newdata, ...)
#' 
#'   # probs must be either a probability matrix or NULL if unsupported
#'   probs <- NULL
#' 
#'   list(class = as.factor(raw_preds), probs = probs)
#' }
#' 
#' # 2. Use it natively in the pipeline
#' # my_model <- fit_model(data, formula, classifier = mypackage::myfit)
#' ```
#'
#' @name adapter_interface
NULL
