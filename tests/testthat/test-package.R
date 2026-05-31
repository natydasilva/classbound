test_that("classbound package loads", {
  expect_true(requireNamespace("classbound", quietly = TRUE))
})

test_that("classbound has required dependencies", {
  expect_true(requireNamespace("ggplot2", quietly = TRUE))
  expect_true(requireNamespace("shiny", quietly = TRUE))
  expect_true(requireNamespace("MASS", quietly = TRUE))
})
