test_that("plot_boundary works correctly with basic boundary", {
  boundary_data <- data.frame(
    x = rep(1:10, 10),
    y = rep(1:10, each = 10),
    prediction = factor(sample(c("A", "B"), 100, replace = TRUE))
  )

  mock_model <- structure(list(boundary_data = boundary_data), class = "classbound")
  p <- plot_boundary(mock_model)
  expect_s3_class(p, "ggplot")

  # Verify the plot builds successfully without error
  expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
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

  mock_model <- structure(list(boundary_data = boundary_data), class = "classbound")
  p <- plot_boundary(
    mock_model,
    obs_data = obs_data,
    x_col = "feat1",
    y_col = "feat2",
    true_label = "label"
  )

  expect_s3_class(p, "ggplot")

  # Verify the plot builds successfully without error
  expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
})

test_that("plot_boundary handles errors gracefully", {
  boundary_data <- data.frame(x = 1:5, y = 1:5, prediction = factor(1:5))

  mock_model <- structure(list(boundary_data = boundary_data), class = "classbound")

  # Type not 2D
  expect_error(plot_boundary(mock_model, type = "tour"), "Only type='2D' and type='disagreement' are currently supported")

  # Missing prediction column
  bad_boundary <- data.frame(x = 1:5, y = 1:5)
  bad_model <- structure(list(boundary_data = bad_boundary), class = "classbound")
  expect_error(plot_boundary(bad_model), "boundary must contain 'x', 'y', and 'prediction' columns")

  # Missing obs_data columns mapping
  expect_error(
    plot_boundary(mock_model, obs_data = data.frame(a = 1, b = 2)),
    "you must specify true_label"
  )

  # Missing obs_data columns actual
  expect_error(
    plot_boundary(mock_model, obs_data = data.frame(a = 1, b = 2), x_col = "x", y_col = "y", true_label = "z"),
    "obs_data must contain column"
  )
})

test_that("plot.classbound_boundary S3 method delegates to plot_boundary", {
  boundary_data <- data.frame(
    x = rep(1:10, 10),
    y = rep(1:10, each = 10),
    prediction = factor(sample(c("A", "B"), 100, replace = TRUE))
  )

  mock_model <- structure(list(boundary_data = boundary_data), class = c("classbound_boundary", "classbound"))


  # Call plot() using S3 dispatch
  p <- plot(mock_model)

  expect_s3_class(p, "ggplot")

  # Verify the plot builds successfully without error
  expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))
})
