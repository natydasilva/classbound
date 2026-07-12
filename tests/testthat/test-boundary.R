test_that("boundary_compute works correctly with rpart model", {
  skip_if_not_installed("rpart")
  # Use a 2D subset of penguins to avoid missing columns in predict
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins[, c("bill_length_mm", "bill_depth_mm", "species")]

  model <- fit_model(
    data = train_data,
    formula = species ~ .,
    classifier = rpart::rpart
  )

  # Define range
  feature_range <- list(
    bill_length_mm = c(30.0, 60.0),
    bill_depth_mm = c(10.0, 25.0)
  )

  # Compute boundary
  res_model <- boundary_compute(model, range = feature_range, resolution = 10)
  res <- res_model$boundary_data

  # Check structure
  expect_s3_class(res, "data.frame")
  expect_true(all(c("x", "y", "prediction") %in% colnames(res)))

  # 10 * 10 grid = 100 points
  expect_equal(nrow(res), 100)
  expect_true(is.factor(res$prediction))

  # Check probabilities exist (rpart provides them)
  # levels of species are Adelie, Chinstrap, Gentoo
  expect_true(all(c("Adelie", "Chinstrap", "Gentoo") %in% colnames(res)))

  # Check probabilities sum to 1 row-wise
  probs_sum <- rowSums(res[, c("Adelie", "Chinstrap", "Gentoo")])
  expect_true(all(abs(probs_sum - 1) < 1e-6))
})

test_that("boundary_compute handles errors gracefully", {
  skip_if_not_installed("rpart")
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins[, c("bill_length_mm", "bill_depth_mm", "species")]
  model <- fit_model(data = train_data, formula = species ~ ., classifier = rpart::rpart)

  # Not a classbound object
  expect_error(
    boundary_compute(model$fit, range = list(x = c(1, 2), y = c(1, 2)), resolution = 10),
    "model must be a 'classbound' object"
  )

  # Invalid range
  expect_error(
    boundary_compute(model, range = c(1, 2, 3), resolution = 10),
    "range must be a named list of length 2"
  )

  # Unnamed list
  expect_error(
    boundary_compute(model, range = list(c(1, 2), c(1, 2)), resolution = 10),
    "range must be a named list of length 2"
  )
})

test_that("boundary_compute handles metadata validation gracefully", {
  skip_if_not_installed("rpart")
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins[, c("bill_length_mm", "bill_depth_mm", "species")]
  model <- fit_model(data = train_data, formula = species ~ ., classifier = rpart::rpart)

  # Duplicate range names
  expect_error(
    boundary_compute(model, range = list(bill_length_mm = c(1, 2), bill_length_mm = c(1, 2))),
    "Duplicate feature names found in `range`."
  )

  # Invalid name & Missing name
  expect_error(
    boundary_compute(model, range = list(bill_length_mm = c(1, 2), wrong_name = c(1, 2))),
    "Names in `range` do not match the training features."
  )

  expect_error(
    boundary_compute(model, range = list(bill_length_mm = c(1, 2), wrong_name = c(1, 2))),
    "Invalid features: wrong_name"
  )

  # Model trained on 3 features, range provides 2
  train_data_3 <- penguins[, c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "species")]
  model_3 <- fit_model(data = train_data_3, formula = species ~ ., classifier = rpart::rpart)

  expect_error(
    boundary_compute(model_3, range = list(bill_length_mm = c(1, 2), bill_depth_mm = c(1, 2))),
    "Visualizing models with >2 features without a projection requires fixed-value slicing"
  )
})

test_that("boundary_compute rejects categorical features", {
  skip_if_not_installed("rpart")
  library(palmerpenguins)
  penguins <- na.omit(penguins)

  # Train model on a numeric and a character/factor column
  train_data <- penguins[, c("bill_length_mm", "island", "species")]
  model <- fit_model(data = train_data, formula = species ~ bill_length_mm + island, classifier = rpart::rpart)

  expect_error(
    boundary_compute(model, range = list(bill_length_mm = c(30, 60), island = c(1, 2)), resolution = 10),
    "Boundary generation requires numeric features. Categorical ranges are not supported."
  )
})
