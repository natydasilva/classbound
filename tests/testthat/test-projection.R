test_that("boundary_compute correctly saves projection", {
  # Simple mock projection
  basis <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 3, ncol = 2)
  rownames(basis) <- c("V1", "V2", "V3")
  proj <- list(
    basis = basis,
    center = c(0, 0, 0),
    scale = c(1, 1, 1)
  )

  # Mock model
  data <- data.frame(V1 = rnorm(10), V2 = rnorm(10), V3 = rnorm(10), Y = factor(sample(c("A", "B"), 10, replace = TRUE)))
  model <- fit_model(data, Y ~ V1 + V2 + V3, rpart::rpart)

  # Compute boundary with projection
  res <- boundary_compute(model, range = list(PC1 = c(-2, 2), PC2 = c(-2, 2)), projection = proj)

  # Check that projection is saved at the top level
  expect_false(is.null(res$projection))
  expect_equal(res$projection, proj)
})

test_that("plot_boundary handles projection forward mapping correctly", {
  basis <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 3, ncol = 2)
  rownames(basis) <- c("V1", "V2", "V3")
  proj <- list(
    basis = basis,
    center = c(0, 0, 0),
    scale = c(1, 1, 1)
  )

  data <- data.frame(V1 = rnorm(10), V2 = rnorm(10), V3 = rnorm(10), Y = factor(sample(c("A", "B"), 10, replace = TRUE)))
  model <- fit_model(data, Y ~ V1 + V2 + V3, rpart::rpart)
  model <- boundary_compute(model, range = list(PC1 = c(-2, 2), PC2 = c(-2, 2)), projection = proj)

  # Plot with obs_data
  p <- plot_boundary(model, obs_data = data, true_label = "Y")

  # Verify ggplot layer aesthetics
  expect_s3_class(p, "ggplot")

  # The second layer should be the points (first is raster)
  point_data <- p$layers[[2]]$data

  # Check if alpha_val exists and is within [0, 1]
  expect_true("alpha_val" %in% colnames(point_data))
  expect_true(all(round(point_data$alpha_val, 5) >= 0.2 & round(point_data$alpha_val, 5) <= 1.0))

  # Check if x_val and y_val were synthesized
  expect_true("x_val" %in% colnames(point_data))
  expect_true("y_val" %in% colnames(point_data))
})

test_that("plot_boundary fails gracefully if projection features are missing", {
  basis <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 3, ncol = 2)
  rownames(basis) <- c("V1", "V2", "V3")
  proj <- list(
    basis = basis,
    center = c(0, 0, 0),
    scale = c(1, 1, 1)
  )

  data <- data.frame(V1 = rnorm(10), V2 = rnorm(10), V3 = rnorm(10), Y = factor(sample(c("A", "B"), 10, replace = TRUE)))
  model <- fit_model(data, Y ~ V1 + V2 + V3, rpart::rpart)
  model <- boundary_compute(model, range = list(PC1 = c(-2, 2), PC2 = c(-2, 2)), projection = proj)

  bad_data <- data.frame(V1 = rnorm(10), Y = factor(sample(c("A", "B"), 10, replace = TRUE)))

  expect_error(
    plot_boundary(model, obs_data = bad_data, true_label = "Y"),
    "missing features required for the projection basis"
  )
})
