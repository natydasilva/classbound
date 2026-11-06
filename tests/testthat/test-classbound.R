test_that("classbound() fits, computes and plots in one step", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("palmerpenguins")

  d <- as.data.frame(na.omit(palmerpenguins::penguins[, c("species", "bill_length_mm", "bill_depth_mm")]))
  p <- classbound(d, species ~ bill_length_mm + bill_depth_mm, rpart::rpart, resolution = 12)

  expect_true(inherits(p, "gg"))
  expect_equal(nrow(p$data), 12 * 12)
  expect_equal(range(p$data$x), range(d$bill_length_mm))
  expect_equal(range(p$data$y), range(d$bill_depth_mm))
  expect_setequal(levels(p$data$prediction), levels(d$species))
})

test_that("classbound() uses a projection for more than two predictors", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("palmerpenguins")

  d <- as.data.frame(na.omit(palmerpenguins::penguins[, c("species", "bill_length_mm", "bill_depth_mm", "flipper_length_mm")]))
  expect_error(
    classbound(d, species ~ ., rpart::rpart, resolution = 5),
    "projection"
  )

  pca <- stats::prcomp(d[, -1], scale. = TRUE)
  proj <- list(basis = pca$rotation[, 1:2], center = pca$center, scale = pca$scale)
  p <- classbound(d, species ~ ., rpart::rpart, projection = proj, resolution = 8)
  expect_equal(nrow(p$data), 8 * 8)
  expect_equal(range(p$data$x), range(pca$x[, 1]), tolerance = 1e-8)
})
