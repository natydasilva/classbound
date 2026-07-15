test_that("boundary_compute slicing works with imputation", {
  skip_if_not_installed("rpart")
  library(palmerpenguins)
  penguins <- na.omit(penguins)

  # Train on 4 numeric features
  train_data <- penguins[, c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g", "species")]
  model <- fit_model(data = train_data, formula = species ~ ., classifier = rpart::rpart)

  # Request 2D boundary on bill_length and bill_depth
  boundary_res <- boundary_compute(
    model,
    range = list(bill_length_mm = c(30, 60), bill_depth_mm = c(10, 25)),
    resolution = 10
  )

  # Verify boundary features are stored
  expect_equal(boundary_res$boundary_features, c("bill_length_mm", "bill_depth_mm"))

  # The other two features should be imputed with their medians
  expected_flipper <- median(train_data$flipper_length_mm, na.rm = TRUE)
  expected_mass <- median(train_data$body_mass_g, na.rm = TRUE)

  expect_equal(model$metadata$features$imputation$flipper_length_mm, expected_flipper)
  expect_equal(model$metadata$features$imputation$body_mass_g, expected_mass)

  # Predictions should exist
  expect_true("prediction" %in% colnames(boundary_res$boundary_data))
})

test_that("boundary_compute slicing works with explicit reference", {
  skip_if_not_installed("rpart")
  library(palmerpenguins)
  penguins <- na.omit(penguins)

  train_data <- penguins[, c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g", "island", "species")]
  model <- fit_model(data = train_data, formula = species ~ ., classifier = rpart::rpart)

  # Request 2D boundary with references
  boundary_res <- boundary_compute(
    model,
    range = list(bill_length_mm = c(30, 60), bill_depth_mm = c(10, 25)),
    reference = list(flipper_length_mm = 200, body_mass_g = 4000, island = "Biscoe"),
    resolution = 10
  )

  expect_equal(boundary_res$boundary_features, c("bill_length_mm", "bill_depth_mm"))
  expect_true("prediction" %in% colnames(boundary_res$boundary_data))
})

test_that("boundary_compute slicing validation catches errors", {
  skip_if_not_installed("rpart")
  library(palmerpenguins)
  penguins <- na.omit(penguins)

  train_data <- penguins[, c("bill_length_mm", "bill_depth_mm", "island", "species")]
  model <- fit_model(data = train_data, formula = species ~ ., classifier = rpart::rpart)

  # Invalid reference feature name
  expect_error(
    boundary_compute(model, range = list(bill_length_mm = c(30, 60), bill_depth_mm = c(10, 25)), reference = list(fake_feat = 10)),
    "Names in `reference` do not match training features."
  )

  # Invalid factor level
  expect_error(
    boundary_compute(model, range = list(bill_length_mm = c(30, 60), bill_depth_mm = c(10, 25)), reference = list(island = "FakeIsland")),
    "Reference value 'FakeIsland' for feature 'island' is not a valid level."
  )
})
