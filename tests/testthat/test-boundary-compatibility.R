test_that("boundary_compute works correctly with pptree model", {
  skip_if_not_installed("PPtreeExt")
  suppressWarnings(suppressPackageStartupMessages(library(palmerpenguins)))
  penguins <- as.data.frame(na.omit(penguins[, -c(2, 7, 8)]))
  train_data <- penguins[, c("bill_length_mm", "bill_depth_mm", "species")]

  model <- fit_model(data = train_data, formula = species ~ ., classifier = PPtreeExt::PPtreeExtclass)

  feature_range <- list(
    bill_length_mm = c(30.0, 60.0),
    bill_depth_mm = c(10.0, 25.0)
  )

  res_model <- boundary_compute(model, feature_range = feature_range, resolution = 10)
  res <- res_model$boundary_data

  expect_s3_class(res, "data.frame")
  expect_true(all(c("x", "y", "prediction") %in% colnames(res)))
  expect_equal(nrow(res), 100)
  expect_true(is.factor(res$prediction))

  # PPtreeExt does not provide probabilities natively in our adapter
  expect_false("Adelie" %in% colnames(res))
})

test_that("boundary_compute works correctly with randomForest model", {
  skip_if_not_installed("randomForest")
  library(palmerpenguins)
  penguins <- as.data.frame(na.omit(penguins[, -c(2, 7, 8)]))
  train_data <- penguins[, c("bill_length_mm", "bill_depth_mm", "species")]

  model <- fit_model(data = train_data, formula = species ~ ., classifier = randomForest::randomForest)

  feature_range <- list(
    bill_length_mm = c(30.0, 60.0),
    bill_depth_mm = c(10.0, 25.0)
  )

  res_model <- boundary_compute(model, feature_range = feature_range, resolution = 10)
  res <- res_model$boundary_data

  expect_s3_class(res, "data.frame")
  expect_true(all(c("x", "y", "prediction") %in% colnames(res)))
  expect_equal(nrow(res), 100)
  expect_true(is.factor(res$prediction))

  # Check probabilities exist
  expect_true(all(c("Adelie", "Chinstrap", "Gentoo") %in% colnames(res)))
  probs_sum <- rowSums(res[, c("Adelie", "Chinstrap", "Gentoo")])
  expect_true(all(abs(probs_sum - 1) < 1e-6))
})

test_that("boundary_compute works correctly with ppforest2 model", {
  skip_if_not_installed("ppforest2")
  library(palmerpenguins)
  penguins <- as.data.frame(na.omit(penguins[, -c(2, 7, 8)]))
  train_data <- penguins[, c("bill_length_mm", "bill_depth_mm", "species")]

  model <- fit_model(
    data = train_data,
    formula = species ~ bill_length_mm + bill_depth_mm,
    classifier = ppforest2::pprf
  )

  feature_range <- list(
    bill_length_mm = c(30.0, 60.0),
    bill_depth_mm = c(10.0, 25.0)
  )

  res_model <- boundary_compute(model, feature_range = feature_range, resolution = 10)
  res <- res_model$boundary_data

  expect_s3_class(res, "data.frame")
  expect_true(all(c("x", "y", "prediction") %in% colnames(res)))
  expect_equal(nrow(res), 100)
  expect_true(is.factor(res$prediction))

  # ppforest2 adapter provides probabilities natively
  expect_true(all(c("Adelie", "Chinstrap", "Gentoo") %in% colnames(res)))
  probs_sum <- rowSums(res[, c("Adelie", "Chinstrap", "Gentoo")])
  expect_true(all(abs(probs_sum - 1) < 1e-6))
})

test_that("boundary_compute works correctly with PPtreeViz model", {
  skip_if_not_installed("PPtreeViz")
  library(palmerpenguins)
  penguins <- as.data.frame(na.omit(penguins[, -c(2, 7, 8)]))
  train_data <- penguins[, c("bill_length_mm", "bill_depth_mm", "species")]

  model <- fit_model(data = train_data, formula = species ~ ., classifier = PPtreeViz::PPTreeclass)

  feature_range <- list(
    bill_length_mm = c(30.0, 60.0),
    bill_depth_mm = c(10.0, 25.0)
  )

  res_model <- boundary_compute(model, feature_range = feature_range, resolution = 10)
  res <- res_model$boundary_data

  expect_s3_class(res, "data.frame")
  expect_true(all(c("x", "y", "prediction") %in% colnames(res)))
  expect_equal(nrow(res), 100)
  expect_true(is.factor(res$prediction))

  # PPtreeViz adapter does not provide probabilities
  expect_false("Adelie" %in% colnames(res))
})
