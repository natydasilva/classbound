test_that("preprocess_data validates inputs correctly", {
  # Non-dataframe
  expect_error(
    preprocess_data(c(1, 2, 3)),
    "must be a data.frame"
  )

  # Length mismatch
  expect_error(
    preprocess_data(data.frame(x = 1:5), labels = c("A", "B")),
    "Length of `labels` must match"
  )

  # Missing values in data
  expect_error(
    preprocess_data(data.frame(x = c(1, NA, 3))),
    "Missing values are not supported"
  )

  # Missing values in labels
  expect_error(
    preprocess_data(data.frame(x = 1:3), labels = c("A", NA, "B")),
    "Missing values are not supported"
  )
})

test_that("preprocess_data handles factors correctly", {
  df <- data.frame(
    num = 1:3,
    char = c("red", "blue", "red"),
    fct = factor(c("small", "large", "small"), levels = c("small", "medium", "large")),
    stringsAsFactors = FALSE
  )
  labels <- factor(c("yes", "no", "yes"), levels = c("yes", "no", "maybe"))

  processed <- preprocess_data(df, labels)

  # Labels should be factor with 'maybe' dropped
  expect_true(is.factor(processed$labels))
  expect_equal(levels(processed$labels), c("yes", "no"))

  # char column should be converted to factor
  expect_true(is.factor(processed$data$char))
  expect_equal(levels(processed$data$char), c("blue", "red"))

  # fct column should have 'medium' dropped
  expect_true(is.factor(processed$data$fct))
  expect_equal(levels(processed$data$fct), c("small", "large"))

  # num column should remain numeric
  expect_true(is.numeric(processed$data$num))
})

test_that("preprocess_data handles scaling correctly", {
  df <- data.frame(
    num1 = c(10, 20, 30),
    num2 = c(100, 200, 300),
    char = c("A", "B", "C")
  )

  processed <- preprocess_data(df, scale = TRUE)

  # numeric columns should be scaled
  expect_equal(mean(processed$data$num1), 0)
  expect_equal(sd(processed$data$num1), 1)
  expect_equal(mean(processed$data$num2), 0)
  expect_equal(sd(processed$data$num2), 1)

  # character column should be factor but unaffected by scaling
  expect_true(is.factor(processed$data$char))
})
