library(testthat)

# Use palmerpenguins for tests, drop NAs, and keep only numeric features for projection
penguins_data <- na.omit(palmerpenguins::penguins)
# Select Bill Length, Bill Depth, Flipper Length, Body Mass (4 numeric features)
features <- c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g")
penguins_data <- penguins_data[, c(features, "species")]

test_that("boundary_compute validates projection object correctly", {
  skip_if_not_installed("rpart")
  model <- fit_model(penguins_data, species ~ ., rpart::rpart)

  # Missing projection for high-dimensional model
  expect_error(
    boundary_compute(model, list(z1 = c(-1, 1), z2 = c(-1, 1))),
    "Visualizing models with >2 features without a projection requires fixed-value slicing"
  )

  # Model with a factor feature
  penguins_with_factor <- na.omit(palmerpenguins::penguins)[, c(features, "sex", "species")]
  model_factor <- fit_model(penguins_with_factor, species ~ ., rpart::rpart)
  expect_error(
    boundary_compute(model_factor, list(z1 = c(-1, 1), z2 = c(-1, 1)), projection = list(basis = matrix(1))),
    "High-dimensional projection is only supported for models trained exclusively on numeric features."
  )

  # Invalid projection type
  expect_error(
    boundary_compute(model, list(z1 = c(-1, 1), z2 = c(-1, 1)), projection = "not a list"),
    "`projection` must be a list containing at least a `basis` matrix"
  )

  # Missing basis
  expect_error(
    boundary_compute(model, list(z1 = c(-1, 1), z2 = c(-1, 1)), projection = list(center = rep(0, 4))),
    "`projection` must be a list containing at least a `basis` matrix"
  )

  # Basis not a matrix
  expect_error(
    boundary_compute(model, list(z1 = c(-1, 1), z2 = c(-1, 1)), projection = list(basis = 1:8)),
    "`projection\\$basis` must be a numeric matrix"
  )

  # Basis wrong dimensions (e.g. 3x2 instead of 4x2)
  bad_basis <- matrix(runif(6), nrow = 3, ncol = 2)
  expect_error(
    boundary_compute(model, list(z1 = c(-1, 1), z2 = c(-1, 1)), projection = list(basis = bad_basis)),
    "`projection\\$basis` must be a 4 x 2 matrix to match training features"
  )

  # Center wrong length
  good_basis <- qr.Q(qr(matrix(runif(8), nrow = 4, ncol = 2)))
  expect_error(
    boundary_compute(
      model, list(z1 = c(-1, 1), z2 = c(-1, 1)),
      projection = list(basis = good_basis, center = c(1, 2, 3))
    ),
    "`center` must be a numeric vector of length 4"
  )

  # NAs in basis
  basis_na <- good_basis
  basis_na[1, 1] <- NA
  expect_error(
    boundary_compute(model, list(z1 = c(-1, 1), z2 = c(-1, 1)), projection = list(basis = basis_na)),
    "`projection\\$basis` must be a numeric matrix without missing values"
  )

  # NAs in center
  expect_error(
    boundary_compute(
      model, list(z1 = c(-1, 1), z2 = c(-1, 1)),
      projection = list(basis = good_basis, center = c(1, 2, NA, 4))
    ),
    "`center` must be a numeric vector of length 4 without missing values"
  )

  # Scale wrong length
  expect_error(
    boundary_compute(
      model, list(z1 = c(-1, 1), z2 = c(-1, 1)),
      projection = list(basis = good_basis, scale = c(1, 2, 3))
    ),
    "`scale` must be a numeric vector of length 4"
  )

  # NAs in scale
  expect_error(
    boundary_compute(
      model, list(z1 = c(-1, 1), z2 = c(-1, 1)),
      projection = list(basis = good_basis, scale = c(1, 2, NA, 4))
    ),
    "`scale` must be a numeric vector of length 4 without missing values"
  )

  # Incorrect rownames in basis
  basis_wrong_names <- good_basis
  rownames(basis_wrong_names) <- c("a", "b", "c", "d")
  expect_error(
    boundary_compute(model, list(z1 = c(-1, 1), z2 = c(-1, 1)), projection = list(basis = basis_wrong_names)),
    "Row names of `projection\\$basis` must match the exact feature names and ordering"
  )
  # Non-orthonormal basis error
  non_ortho_basis <- matrix(runif(8), nrow = 4, ncol = 2)
  expect_error(
    boundary_compute(model, list(z1 = c(-1, 1), z2 = c(-1, 1)), projection = list(basis = non_ortho_basis)),
    "The projection basis is not orthonormal"
  )
})

test_that("boundary_compute works with a valid projection", {
  skip_if_not_installed("rpart")
  model <- fit_model(penguins_data, species ~ ., rpart::rpart)

  # Create a mock orthonormal basis (e.g., from PCA)
  # Simple basis: feature 1 and 2 mapped to z1, feature 3 and 4 mapped to z2
  basis <- matrix(c(
    1, 0, 0, 0,
    0, 1, 0, 0
  ), nrow = 4, ncol = 2)
  center <- colMeans(penguins_data[, features])

  projection <- list(basis = basis, center = center)
  range_z <- list(z1 = c(-10, 10), z2 = c(-5, 5))

  res_model <- boundary_compute(model, range = range_z, resolution = 10, projection = projection)
  res <- res_model$boundary_data

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 100) # 10x10 resolution
  expect_true(all(c("x", "y", "prediction") %in% colnames(res)))

  # Check if probability columns are present
  expect_true(all(levels(penguins_data$species) %in% colnames(res)))

  # Test without center
  res_no_center_model <- boundary_compute(model, range = range_z, resolution = 10, projection = list(basis = basis))
  res_no_center <- res_no_center_model$boundary_data
  expect_s3_class(res_no_center, "data.frame")
  expect_equal(nrow(res_no_center), 100)
})

test_that("projection reduces mathematically to 2D workflow", {
  skip_if_not_installed("rpart")
  # Train a 2D model
  penguins_2d <- na.omit(palmerpenguins::penguins)[, c("bill_length_mm", "bill_depth_mm", "species")]
  model_2d <- fit_model(penguins_2d, species ~ ., rpart::rpart)

  range_2d <- list(
    bill_length_mm = c(30, 60),
    bill_depth_mm = c(10, 25)
  )

  # Compute standard 2D boundary
  res_standard_model <- boundary_compute(model_2d, range = range_2d, resolution = 10)
  res_standard <- res_standard_model$boundary_data

  # Compute boundary via inverse projection identity mapping
  projection_identity <- list(
    basis = diag(2),
    center = c(0, 0)
  )
  # Rownames must match exactly if provided, or be NULL
  rownames(projection_identity$basis) <- c("bill_length_mm", "bill_depth_mm")

  res_projected_model <- boundary_compute(
    model_2d,
    range = range_2d,
    resolution = 10,
    projection = projection_identity
  )
  res_projected <- res_projected_model$boundary_data

  # Ensure mathematically identical predictions (ignoring column names since projection outputs expected_names)
  # Actually, both should output a dataframe with columns (x, y, prediction, ...)
  # The projected dataframe will have x and y as expected_names if projection was not NULL, wait no,
  # wait: if projection is NULL, x=seq_x, y=seq_y. If projection is not NULL, x=seq_x, y=seq_y.
  expect_equal(res_standard$prediction, res_projected$prediction)
  expect_equal(res_standard$x, res_projected$x)
  expect_equal(res_standard$y, res_projected$y)
})
