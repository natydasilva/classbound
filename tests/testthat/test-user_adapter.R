test_that("user-defined custom adapters are supported natively via S3 dispatch", {
  library(palmerpenguins)
  train_data <- na.omit(palmerpenguins::penguins[, c("bill_length_mm", "bill_depth_mm", "species")])

  # 1. Define a custom classifier function
  my_custom_classifier <- function(formula, data, ...) {
    structure(
      list(formula = formula, dummy_labels = data[[all.vars(formula[[2]])]], fitted = TRUE),
      class = "my_model"
    )
  }

  # 2. Define a custom prediction adapter for that class
  predict_adapter.my_model <- function(model, newdata, ...) { # nolint: object_name_linter.
    # Just predict the first level for all rows
    n <- nrow(newdata)
    level_1 <- levels(model$dummy_labels)[1]
    preds <- factor(rep(level_1, n), levels = levels(model$dummy_labels))
    list(class = preds, probs = NULL)
  }

  # 3. Register the S3 method locally so testthat environment can see them
  local({
    registerS3method("predict_adapter", "my_model", predict_adapter.my_model, envir = parent.frame())
  })

  # 4. Use the core pipeline to fit the model using the custom function
  my_fit <- fit_model(train_data, species ~ ., classifier = my_custom_classifier)

  # Validate that the model object was properly constructed
  expect_s3_class(my_fit, "classbound")
  expect_s3_class(my_fit$fit, "my_model")
  expect_true(my_fit$fit$fitted)

  # 5. Use the core pipeline to generate a boundary (which tests predict)
  grid <- boundary_compute(my_fit, list(bill_length_mm = c(30, 60), bill_depth_mm = c(10, 20)), resolution = 10)

  # Validate that boundary generation succeeded using the custom predict adapter
  expect_s3_class(grid, "data.frame")
  expect_equal(nrow(grid), 100)
  expect_true("prediction" %in% colnames(grid))
  expect_true(is.factor(grid$prediction))
})
