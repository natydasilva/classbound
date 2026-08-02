test_that("generate_extreme_outlier creates deterministic extreme outliers", {
  # Simple 2D numeric dataset
  df <- data.frame(
    Sim = c("A", "A", "B", "B"),
    X1 = c(1, 2, 8, 9),
    X2 = c(1, 2, 8, 9)
  )
  
  outlier <- generate_extreme_outlier(df, "C", magnitude = 2, target_col = "Sim")
  
  expect_equal(nrow(outlier), 1)
  expect_equal(outlier$Sim, "C")
  
  # Max of X1 is 9, range is 8. 
  # Class C is index 3. sub-step = 2 * 0.03 * 2 = 0.12. 
  # Magnitude = 2 + 0.12 = 2.12
  # Expansion = 2.12 * 8 * 0.1 = 1.696
  # Expected X1 = 9 + 1.696 = 10.696
  expect_equal(outlier$X1, 10.696)
  expect_equal(outlier$X2, 10.696)
})

test_that("generate_extreme_outlier handles high-dimensional and categorical data", {
  df <- data.frame(
    Sim = c("A", "B", "A", "B", "A"),
    X1 = c(1, 2, 3, 4, 5), # Numeric, range 4
    X2 = c(10, 20, 30, 40, 50), # Numeric, range 40
    X3 = c(100, 200, 300, 400, 500), # Numeric > 2nd dimension
    Cat = factor(c("Yes", "No", "Yes", "Yes", "No"), levels = c("Yes", "No", "Maybe"))
  )
  
  outlier <- generate_extreme_outlier(df, "B", magnitude = 5, target_col = "Sim")
  
  # B is class index 2. sub-step = 1 * 0.03 * 5 = 0.15. Magnitude = 5.15
  # X1: max 5, range 4 -> 5 + (5.15 * 4 * 0.1) = 7.06
  expect_equal(outlier$X1, 7.06)
  
  # X2: min 10, range 40, magnitude 5.15 -> 10 - (5.15 * 40 * 0.1) = -10.6
  expect_equal(outlier$X2, -10.6)
  
  # X3: > 2nd dimension, should use median (300)
  expect_equal(outlier$X3, 300)
  
  # Cat: Mode is "Yes", should remain a factor with same levels
  expect_equal(as.character(outlier$Cat), "Yes")
  expect_true(is.factor(outlier$Cat))
  expect_equal(levels(outlier$Cat), c("Yes", "No", "Maybe"))
})

test_that("generate_extreme_outlier handles zero variance gracefully", {
  df <- data.frame(
    Sim = c("A", "B"),
    X1 = c(5, 5),
    X2 = c(10, 10)
  )
  
  outlier <- generate_extreme_outlier(df, "A", magnitude = 1)
  # Range is 0, so it uses 1. Min is 5, magnitude 1 -> 5 - (1 * 1 * 0.1) = 4.9
  expect_equal(outlier$X1, 4.9)
})

test_that("create_boundary_plot generates a ggplot object for 2D data", {
  df <- data.frame(
    Sim = factor(c("A", "B", "A", "B")),
    X1 = c(1, 2, 8, 9),
    X2 = c(1, 2, 8, 9)
  )
  m <- fit_model(Sim ~ ., data = df, classifier = rpart::rpart)
  
  p <- create_boundary_plot(m, df, "Test Plot", resolution = 10)
  
  expect_true(inherits(p, "ggplot"))
  expect_equal(p$labels$title, "Test Plot")
})
