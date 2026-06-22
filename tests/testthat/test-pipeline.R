test_that("fit_model handles invalid inputs and unsupported methods gracefully", {
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins[, 2:5]
  train_labels <- penguins$species

  # Missing method
  expect_error(
    fit_model(data = train_data, labels = train_labels),
    "Please specify a classification method"
  )

  # Unsupported method
  expect_error(
    fit_model(data = train_data, labels = train_labels, method = "unsupported_magic_model"),
    "not supported"
  )
})

test_that("predict_model correctly enforces structural constraints", {
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  raw_rpart <- rpart::rpart(species ~ ., data = penguins)

  expect_error(
    predict_model(raw_rpart, newdata = penguins),
    "must be a 'classbound_model'"
  )
})

test_that("boundary_compute handles classifiers with probs = NULL gracefully", {
  # Mock a simple classifier adapter that returns probs = NULL
  fit_adapter.classbound_mock_noprobs <- function(object, data, labels, ...) {
    list(dummy_model = TRUE)
  }
  
  predict_model.classbound_mock_noprobs <- function(model, newdata, ...) {
    # Always predict the first level for simplicity
    class_preds <- factor(rep("A", nrow(newdata)), levels = c("A", "B"))
    list(class = class_preds, probs = NULL)
  }
  
  # Register the S3 methods temporarily
  local({
    registerS3method("fit_adapter", "classbound_mock_noprobs", fit_adapter.classbound_mock_noprobs, envir = parent.frame())
    registerS3method("predict_model", "classbound_mock_noprobs", predict_model.classbound_mock_noprobs, envir = parent.frame())
  })
  
  # Create some dummy data
  data <- data.frame(X1 = c(1, 2, 3), X2 = c(3, 2, 1))
  labels <- factor(c("A", "B", "A"))
  
  # Fit the mock model
  model <- fit_model(data, labels, "mock_noprobs")
  
  # Compute boundary
  grid <- boundary_compute(model, list(X1 = c(0, 4), X2 = c(0, 4)), resolution = 10)
  
  # Validation
  expect_true(is.data.frame(grid))
  expect_equal(nrow(grid), 100)
  expect_true(all(c("x", "y", "prediction") %in% colnames(grid)))
  
  # Ensure no probability columns exist
  expect_false(any(c("A", "B") %in% colnames(grid)))
})
