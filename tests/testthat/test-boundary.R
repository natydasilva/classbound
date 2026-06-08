test_that("boundary_compute works correctly with rpart model", {
  # Use a 2D subset of iris to avoid missing columns in predict
  train_data <- iris[, c("Sepal.Length", "Sepal.Width")]
  train_labels <- iris$Species
  
  model <- fit_model(
    data = train_data,
    labels = train_labels,
    method = "rpart"
  )
  
  # Define range
  feature_range <- list(
    Sepal.Length = c(4.0, 8.0),
    Sepal.Width = c(2.0, 5.0)
  )
  
  # Compute boundary
  res <- boundary_compute(model, range = feature_range, resolution = 10)
  
  # Check structure
  expect_s3_class(res, "data.frame")
  expect_true(all(c("x", "y", "prediction") %in% colnames(res)))
  
  # 10 * 10 grid = 100 points
  expect_equal(nrow(res), 100)
  expect_true(is.factor(res$prediction))
  
  # Check probabilities exist (rpart provides them)
  # levels of Species are setosa, versicolor, virginica
  expect_true(all(c("setosa", "versicolor", "virginica") %in% colnames(res)))
  
  # Check probabilities sum to 1 row-wise
  probs_sum <- rowSums(res[, c("setosa", "versicolor", "virginica")])
  expect_true(all(abs(probs_sum - 1) < 1e-6))
})

test_that("boundary_compute handles errors gracefully", {
  train_data <- iris[, c("Sepal.Length", "Sepal.Width")]
  train_labels <- iris$Species
  model <- fit_model(data = train_data, labels = train_labels, method = "rpart")
  
  # Not a classbound_model
  expect_error(
    boundary_compute(model$model, range = list(x=c(1,2), y=c(1,2)), resolution = 10),
    "model must be a 'classbound_model' object"
  )
  
  # Invalid range
  expect_error(
    boundary_compute(model, range = c(1,2,3), resolution = 10),
    "range must be a named list of length 2"
  )
  
  # Unnamed list
  expect_error(
    boundary_compute(model, range = list(c(1,2), c(1,2)), resolution = 10),
    "range must be a named list of length 2"
  )
})
