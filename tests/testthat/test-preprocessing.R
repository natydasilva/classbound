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


test_that("preprocess_data rejects data frames with fewer than 2 rows", {
  expect_error(
    preprocess_data(data.frame()),
    "must have at least two rows"
  )

  # Test a data.frame with columns but 0 rows
  expect_error(
    preprocess_data(data.frame(x = numeric(0), y = numeric(0))),
    "must have at least two rows"
  )
  
  # Test a single-row data.frame
  expect_error(
    preprocess_data(data.frame(x = 1, y = 2)),
    "must have at least two rows"
  )
})

test_that("preprocess_data rejects infinite values", {
  expect_error(
    preprocess_data(data.frame(x = c(1, Inf, 3))),
    "Infinite values are not supported"
  )
  
  expect_error(
    preprocess_data(data.frame(x = c(1, -Inf, 3))),
    "Infinite values are not supported"
  )

})

test_that("preprocess_data rejects duplicate column names", {
  df <- data.frame(a = 1:3, b = 4:6)
  colnames(df) <- c("x", "x")
  expect_error(
    preprocess_data(df),
    "Duplicate column names"
  )
})
