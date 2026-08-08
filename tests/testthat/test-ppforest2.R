test_that("ppforest2 predict_adapter extracts classes and probs correctly", {
  skip_if_not_installed("ppforest2")
  skip_if_not_installed("palmerpenguins")

  library(ppforest2)
  data(penguins, package = "palmerpenguins")
  df <- na.omit(penguins[, c("bill_length_mm", "bill_depth_mm", "species")])

  model <- pprf(species ~ bill_length_mm + bill_depth_mm, data = df, size = 10)

  preds <- predict_adapter.pprf_classification(model, df)

  expect_type(preds, "list")
  expect_named(preds, c("class", "probs"))
  expect_s3_class(preds$class, "factor")
  expect_true(is.matrix(preds$probs))

  # Verify columns align precisely with class levels in the correct order
  expect_equal(colnames(preds$probs), levels(preds$class))

  expect_equal(nrow(preds$probs), nrow(df))
  expect_equal(length(preds$class), nrow(df))

  b_grid <- expand.grid(
    bill_length_mm = seq(30, 60, length.out = 5),
    bill_depth_mm = seq(10, 25, length.out = 5)
  )

  b_preds <- predict_adapter.pprf_classification(model, b_grid)
  expect_equal(nrow(b_preds$probs), nrow(b_grid))
})

test_that("ppforest2 tidymodels integration works through classbound architecture", {
  skip_if_not_installed("ppforest2")
  skip_if_not_installed("palmerpenguins")
  skip_if_not_installed("parsnip")
  skip_if_not_installed("workflows")

  data(penguins, package = "palmerpenguins")
  df <- na.omit(penguins[, -c(2, 7, 8)])

  spec <- parsnip::set_engine(ppforest2::pp_rand_forest(mode = "classification"), "ppforest2")

  # Test fit_model pipeline
  fitted <- fit_model(data = df, formula = species ~ ., classifier = spec)

  # Test predict_model pipeline
  preds <- predict_model(fitted, df)

  expect_type(preds, "list")
  expect_named(preds, c("class", "probs"))
  expect_s3_class(preds$class, "factor")
  expect_true(is.matrix(preds$probs))

  # Ensure columns align perfectly
  expect_equal(colnames(preds$probs), levels(preds$class))
  expect_equal(nrow(preds$probs), nrow(df))
})
