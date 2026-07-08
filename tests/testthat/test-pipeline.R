test_that("fit_model handles missing inputs gracefully", {
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins

  # Missing classifier
  expect_error(
    fit_model(data = train_data, formula = species ~ .),
    "Please specify a classifier function"
  )

  # Missing formula
  expect_error(
    fit_model(data = train_data, classifier = rpart::rpart),
    "Please specify a formula"
  )
})

test_that("predict correctly enforces structural constraints", {
  skip_if_not_installed("rpart")
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  raw_rpart <- rpart::rpart(species ~ ., data = penguins)

  expect_error(
    predict_model(raw_rpart, newdata = penguins),
    "must be a 'classbound'"
  )
})

test_that("predict defaults to error on complex objects, but succeeds with predfun", {
  # Mock a model with no predict_adapter
  mock_model <- structure(
    list(
      fit = structure(list(), class = "unsupported_magic_model"),
      features = list(),
      class_levels = c("A", "B")
    ),
    class = "classbound"
  )

  # Register a dummy predict method for this test that returns a list
  predict.unsupported_magic_model <<- function(object, ...) list(a = 1, b = 2)
  on.exit(rm("predict.unsupported_magic_model", envir = globalenv()))

  expect_error(
    predict(mock_model, newdata = data.frame()),
    "Please supply `predfun` to extract the desired class predictions"
  )

  # Now test that predfun fixes it
  preds <- predict(mock_model, newdata = data.frame(), predfun = function(m, nd, ...) c("A", "B"))
  expect_equal(as.character(preds$class), c("A", "B"))
})

test_that("predict ensures correct factor levels", {
  skip_if_not_installed("rpart")
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins

  model <- fit_model(train_data, species ~ ., classifier = rpart::rpart)
  preds <- predict(model, train_data)
  expect_equal(levels(preds$class), levels(train_data$species))
})
