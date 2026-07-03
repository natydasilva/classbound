test_that("user-defined custom adapters are supported natively via S3 dispatch", {
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins[, c("bill_length_mm", "bill_depth_mm")]
  train_labels <- penguins$species

  # 1. Define a custom fitting adapter in the current environment
  fit_adapter.classbound_my_custom_model <- function(object, data, labels, ...) {
    # Store just the labels in a dummy model structure for validation
    list(dummy_labels = labels, fitted = TRUE)
  }

  # 2. Define a custom prediction adapter in the current environment
  predict_adapter.classbound_my_custom_model <- function(model, newdata, ...) {
    # Just predict the first level for all rows
    n <- nrow(newdata)
    level_1 <- levels(model$model$dummy_labels)[1]
    preds <- factor(rep(level_1, n), levels = levels(model$model$dummy_labels))
    list(class = preds, probs = NULL)
  }

  # 3. Register the S3 methods locally so testthat environment can see them
  local({
    registerS3method("fit_adapter", "classbound_my_custom_model", fit_adapter.classbound_my_custom_model, envir = parent.frame())
    registerS3method("predict_adapter", "classbound_my_custom_model", predict_adapter.classbound_my_custom_model, envir = parent.frame())
  })

  # 4. Use the core pipeline to fit the model using the custom method string
  my_model <- fit_model(train_data, train_labels, method = "my_custom_model")
  
  # Validate that the model object was properly constructed using the custom adapter
  expect_s3_class(my_model, "classbound_my_custom_model")
  expect_s3_class(my_model, "classbound_model")
  expect_true(my_model$model$fitted)
  
  # 5. Use the core pipeline to generate a boundary (which tests predict_model)
  grid <- boundary_compute(my_model, list(bill_length_mm = c(30, 60), bill_depth_mm = c(10, 20)), resolution = 10)
  
  # Validate that boundary generation succeeded using the custom adapter
  expect_s3_class(grid, "data.frame")
  expect_equal(nrow(grid), 100)
  expect_true("prediction" %in% colnames(grid))
  expect_true(is.factor(grid$prediction))
})
