test_that("tidymodels_bridge creates a multi-model grid for valid engines", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("workflowsets")
  skip_if_not_installed("rpart")

  data <- palmerpenguins::penguins
  data <- na.omit(data[, c("bill_length_mm", "bill_depth_mm", "species")])

  # Test valid engine names
  models <- c("rpart")
  res <- tidymodels_bridge(
    data = data,
    response = "species",
    models = models,
    range = list(bill_length_mm = c(30, 60), bill_depth_mm = c(10, 25)),
    resolution = 10
  )

  expect_s3_class(res, "classbound")
  expect_s3_class(res, "classbound_multi")
  expect_true("model" %in% colnames(res$boundary_data))
})

test_that("tidymodels_bridge rejects empty models list", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("workflowsets")

  data <- palmerpenguins::penguins
  data <- na.omit(data[, c("bill_length_mm", "bill_depth_mm", "species")])

  expect_error(
    tidymodels_bridge(data, "species", character(0)),
    "cannot be empty"
  )
})

test_that("tidymodels_bridge rejects duplicates", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("workflowsets")

  data <- palmerpenguins::penguins
  data <- na.omit(data[, c("bill_length_mm", "bill_depth_mm", "species")])

  expect_error(
    tidymodels_bridge(data, "species", c("rpart", "rpart")),
    "Duplicate models detected"
  )
})

test_that("tidymodels_bridge rejects unsupported models", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("workflowsets")

  data <- palmerpenguins::penguins
  data <- na.omit(data[, -c(2, 7, 8)])

  expect_error(
    tidymodels_bridge(data, "species", c("rpart", "fake_engine")),
    "Unsupported models requested: fake_engine"
  )
})

test_that("tidymodels_bridge auto-computes range robustly", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("workflowsets")
  skip_if_not_installed("rpart")

  # 1. Successful auto-computation
  data <- palmerpenguins::penguins
  data <- na.omit(data[, c("bill_length_mm", "bill_depth_mm", "species")])
  
  res <- tidymodels_bridge(data, "species", c("rpart"), resolution = 5)
  expect_true(is.list(res$boundary_data))
  
  # 2. Rejects > 2 predictors without explicit range
  data_3d <- na.omit(palmerpenguins::penguins[, c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "species")])
  expect_error(
    tidymodels_bridge(data_3d, "species", c("rpart")),
    "exactly 2 numeric features"
  )
  
  # 3. Rejects categorical predictors in auto-computation
  data_cat <- data.frame(
    x1 = as.factor(c("A", "B", "A", "B")),
    x2 = c(1, 2, 3, 4),
    y = as.factor(c("Yes", "No", "Yes", "No"))
  )
  expect_error(
    tidymodels_bridge(data_cat, "y", c("rpart")),
    "exactly 2 numeric features"
  )
  
  # 4. Rejects NA/Inf only data
  data_na <- data.frame(
    x1 = as.numeric(c(NA, NA, NA)),
    x2 = c(1, 2, 3),
    y = as.factor(c("Yes", "No", "Yes"))
  )
  expect_error(
    tidymodels_bridge(data_na, "y", c("rpart")),
    "Could not find numeric range"
  )
})
