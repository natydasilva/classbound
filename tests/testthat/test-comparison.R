test_that(
  "comparison workflow enforces class levels and structures correctly",
  {
    skip_if_not_installed("rpart")
    skip_if_not_installed("randomForest")

    # Prepare a subset that excludes one class entirely
    # penguins has 3 classes: Adelie, Chinstrap, Gentoo
    penguins <- palmerpenguins::penguins
    penguins <- na.omit(penguins[, -c(2, 7, 8)])

    # Drop "Gentoo" entirely from the subset
    subset_idx <- penguins$species != "Gentoo"
    subset_data <- penguins[
      subset_idx, c("bill_length_mm", "bill_depth_mm", "species")
    ]

    models <- list(
      rpart = fit_model(subset_data, species ~ ., classifier = rpart::rpart),
      randomForest = fit_model(
        subset_data, species ~ .,
        classifier = randomForest::randomForest
      )
    )

    # Validate list-of-models structure
    expect_type(models, "list")
    expect_s3_class(models$rpart, "classbound")
    expect_s3_class(models$randomForest, "classbound")

    # Validate that original global levels are preserved in metadata
    expect_equal(models$rpart$metadata$class_levels, c("Adelie", "Chinstrap", "Gentoo"))

    # Generate comparison multi-boundary grid
    range <- list(
      bill_length_mm = c(30, 60),
      bill_depth_mm = c(10, 25)
    )

    grids <- lapply(names(models), function(m_name) {
      grid_model <- boundary_compute(models[[m_name]], range, resolution = 10)
      grid <- grid_model$boundary_data
      grid$model <- m_name
      grid
    })
    multi_grid <- do.call(rbind, grids)

    # Validate combined grid structure
    expect_s3_class(multi_grid, "data.frame")
    expect_true(all(c("x", "y", "prediction", "model") %in% colnames(multi_grid)))

    # Validate that class predictions maintain the global 3 levels
    # This guarantees consistent colors across models even if one was trained
    # on a subset
    expect_true(is.factor(multi_grid$prediction))
    expect_equal(
      levels(multi_grid$prediction), c("Adelie", "Chinstrap", "Gentoo")
    )
  }
)
