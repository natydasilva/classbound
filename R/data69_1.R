#' Waveform Dataset (UCI)
#'
#' @description
#' A subset of the UCI Waveform Database Generator dataset (Version 1). Contains
#' 5000 observations with 21 numeric features and a 3-class response variable.
#' All 21 features contribute to class membership, making this a useful dataset for
#' demonstrating high-dimensional boundary visualization, including 2D slice and
#' projection-based approaches.
#'
#' @details
#' Each observation is generated from one of three waveform classes. All features
#' are numeric; the response variable `Y` is a factor with levels `"1"`, `"2"`, and
#' `"3"`. Load the dataset with `data(data69_1)`.
#'
#' This dataset is included primarily to demonstrate high-dimensional boundary
#' visualization with tools like PCA projection and `tourr`. For introductory
#' 2D examples, the `palmerpenguins::penguins` dataset provides a more
#' accessible alternative.
#'
#' @format A data frame with 5000 rows and 22 columns:
#' \describe{
#' \item{Y}{Class label: a factor with 3 levels (`"1"`, `"2"`, `"3"`)}
#' \item{V1}{Numeric feature variable}
#' \item{V2}{Numeric feature variable}
#' \item{V3}{Numeric feature variable}
#' \item{V4}{Numeric feature variable}
#' \item{V5}{Numeric feature variable}
#' \item{V6}{Numeric feature variable}
#' \item{V7}{Numeric feature variable}
#' \item{V8}{Numeric feature variable}
#' \item{V9}{Numeric feature variable}
#' \item{V10}{Numeric feature variable}
#' \item{V11}{Numeric feature variable}
#' \item{V12}{Numeric feature variable}
#' \item{V13}{Numeric feature variable}
#' \item{V14}{Numeric feature variable}
#' \item{V15}{Numeric feature variable}
#' \item{V16}{Numeric feature variable}
#' \item{V17}{Numeric feature variable}
#' \item{V18}{Numeric feature variable}
#' \item{V19}{Numeric feature variable}
#' \item{V20}{Numeric feature variable}
#' \item{V21}{Numeric feature variable}
#' }
#' @docType data
#' @keywords datasets
#' @name data69_1
#' @usage data(data69_1)
#' @source \url{https://archive.ics.uci.edu/dataset/107/waveform+database+generator+version+1}
#' @examples
#' data(data69_1)
#' dim(data69_1) # 5000 x 22
#' levels(data69_1$Y) # "1" "2" "3"
#' head(data69_1[, 1:5])
NULL
