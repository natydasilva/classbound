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
