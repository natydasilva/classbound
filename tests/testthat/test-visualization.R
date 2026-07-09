test_that("plot_boundary works correctly with basic boundary", {
  boundary_data <- data.frame(
    x = rep(1:10, 10),
    y = rep(1:10, each = 10),
    prediction = factor(sample(c("A", "B"), 100, replace = TRUE))
  )

  p <- plot_boundary(boundary_data)
  expect_s3_class(p, "ggplot")

  # Check if geom_raster is present
  layers <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomRaster" %in% layers)
})

test_that("plot_boundary works correctly with observations", {
  boundary_data <- data.frame(
    x = rep(1:10, 10),
    y = rep(1:10, each = 10),
    prediction = factor(sample(c("A", "B"), 100, replace = TRUE))
  )

  obs_data <- data.frame(
    feat1 = runif(10, 1, 10),
    feat2 = runif(10, 1, 10),
    label = factor(sample(c("A", "B"), 10, replace = TRUE))
  )

  p <- plot_boundary(
    boundary_data,
    obs_data = obs_data,
    x_col = "feat1",
    y_col = "feat2",
    true_label = "label"
  )

  expect_s3_class(p, "ggplot")

  layers <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomPoint" %in% layers)
  expect_true("GeomRaster" %in% layers)
})

test_that("plot_boundary handles errors gracefully", {
  boundary_data <- data.frame(x = 1:5, y = 1:5, prediction = factor(1:5))

  # Type not 2D
  expect_error(plot_boundary(boundary_data, type = "tour"), "Only type='2D' is currently supported")

  # Missing prediction column
  bad_boundary <- data.frame(x = 1:5, y = 1:5)
  expect_error(plot_boundary(bad_boundary), "boundary must contain 'x', 'y', and 'prediction' columns")

  # Missing obs_data columns mapping
  expect_error(
    plot_boundary(boundary_data, obs_data = data.frame(a = 1, b = 2)),
    "you must also specify x_col, y_col, and true_label"
  )

  # Missing obs_data columns actual
  expect_error(
    plot_boundary(boundary_data, obs_data = data.frame(a = 1, b = 2), x_col = "x", y_col = "y", true_label = "z"),
    "obs_data must contain columns"
  )
})
