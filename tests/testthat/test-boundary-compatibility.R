test_that("boundary_compute works correctly with pptree model", {
  skip_if_not_installed("PPtreeExt")
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins[, c("bill_length_mm", "bill_depth_mm")]
  train_labels <- penguins$species
  
  model <- fit_model(data = train_data, labels = train_labels, method = "pptree")
  
  feature_range <- list(
    bill_length_mm = c(30.0, 60.0),
    bill_depth_mm = c(10.0, 25.0)
  )
  
  res <- boundary_compute(model, range = feature_range, resolution = 10)
  
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
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins[, c("bill_length_mm", "bill_depth_mm")]
  train_labels <- penguins$species
  
  model <- fit_model(data = train_data, labels = train_labels, method = "randomForest")
  
  feature_range <- list(
    bill_length_mm = c(30.0, 60.0),
    bill_depth_mm = c(10.0, 25.0)
  )
  
  res <- boundary_compute(model, range = feature_range, resolution = 10)
  
  expect_s3_class(res, "data.frame")
  expect_true(all(c("x", "y", "prediction") %in% colnames(res)))
  expect_equal(nrow(res), 100)
  expect_true(is.factor(res$prediction))
  
  # Check probabilities exist
  expect_true(all(c("Adelie", "Chinstrap", "Gentoo") %in% colnames(res)))
  probs_sum <- rowSums(res[, c("Adelie", "Chinstrap", "Gentoo")])
  expect_true(all(abs(probs_sum - 1) < 1e-6))
})

test_that("boundary_compute works correctly with PPforest model", {
  skip_if_not_installed("PPforest")
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- as.data.frame(penguins[, c("bill_length_mm", "bill_depth_mm")])
  train_labels <- penguins$species
  
  model <- fit_model(data = train_data, labels = train_labels, method = "PPforest", m = 10)
  
  feature_range <- list(
    bill_length_mm = c(30.0, 60.0),
    bill_depth_mm = c(10.0, 25.0)
  )
  
  res <- boundary_compute(model, range = feature_range, resolution = 10)
  
  expect_s3_class(res, "data.frame")
  expect_true(all(c("x", "y", "prediction") %in% colnames(res)))
  expect_equal(nrow(res), 100)
  expect_true(is.factor(res$prediction))
  
  # PPforest adapter does not provide probabilities
  expect_false("Adelie" %in% colnames(res))
})

test_that("boundary_compute works correctly with PPtreeViz model", {
  skip_if_not_installed("PPtreeViz")
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- as.data.frame(penguins[, c("bill_length_mm", "bill_depth_mm")])
  train_labels <- penguins$species
  
  model <- fit_model(data = train_data, labels = train_labels, method = "PPtreeViz")
  
  feature_range <- list(
    bill_length_mm = c(30.0, 60.0),
    bill_depth_mm = c(10.0, 25.0)
  )
  
  res <- boundary_compute(model, range = feature_range, resolution = 10)
  
  expect_s3_class(res, "data.frame")
  expect_true(all(c("x", "y", "prediction") %in% colnames(res)))
  expect_equal(nrow(res), 100)
  expect_true(is.factor(res$prediction))
  
  # PPtreeViz adapter does not provide probabilities
  expect_false("Adelie" %in% colnames(res))
})
