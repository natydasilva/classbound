test_that("explorapp requires at least two numeric features", {
  bad <- data.frame(y = c("a", "b"), x = c(1, 2), f = c("u", "v"))
  expect_error(explorapp(bad, "y"), "2 numeric feature")
})

test_that("imported data keeps rows with missing values only in unused columns", {
  skip_if_not_installed("palmerpenguins")

  pg <- as.data.frame(palmerpenguins::penguins)
  n_complete <- sum(stats::complete.cases(pg[, c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g")]))

  shiny::testServer(explorapp(pg, "species"), {
    expect_equal(nrow(current_data()), n_complete)
    expect_equal(setdiff(colnames(current_data()), "Sim"), c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g", "year"))
  })
})

test_that("the projection plane is recomputed when new data is generated", {
  skip_if_not_installed("MixSim")
  skip_if_not_installed("tourr")

  shiny::testServer(explorapp(), {
    session$setInputs(
      data_mode = "Simulate Data", sim_engine = "mixsim", sim_seed = 1, sim_noise = 0,
      sim_mixsim_k = 3, sim_mixsim_p = 3, sim_mixsim_omega = 0.05, sim_mixsim_n = 200,
      sim_mixsim_do = 1
    )
    expect_equal(nrow(current_projection()), 3)

    # More dimensions: the basis must match the new number of features.
    session$setInputs(sim_mixsim_p = 4, sim_mixsim_do = 2)
    expect_equal(nrow(current_projection()), 4)

    # Same dimensions, new data: the plane must pass through the new data centre.
    session$setInputs(sim_seed = 7, sim_mixsim_do = 3)
    expect_equal(
      unname(current_projection_info()$center),
      unname(colMeans(current_data()[, -1]))
    )
  })
})

test_that("clearing the drawing canvas also clears the undo history", {
  shiny::testServer(explorapp(), {
    session$setInputs(data_mode = "Draw Data", interaction_mode = "Draw Point", draw_class = "Class 1")
    session$setInputs(plot_click = list(x = 1, y = 1))
    session$setInputs(plot_click = list(x = 2, y = 2))
    expect_equal(nrow(current_data()), 2)
    expect_length(undo_history(), 2)

    session$setInputs(confirm_draw_reset = 1)
    expect_equal(nrow(current_data()), 0)
    expect_length(undo_history(), 0)
    expect_equal(class_choices(), c("Class 1", "Class 2", "Class 3"))
  })
})

test_that("custom models can declare probability support", {
  shiny::testServer(
    explorapp(custom_models = list(
      "My model" = list(fn = function(...) NULL, supports_prob = TRUE)
    )),
    {
      session$setInputs(selected_models = "My model")
      html <- as.character(output$prob_surface_ui$html)
      expect_true(grepl('id="show_probs"', html, fixed = TRUE))
      expect_false(grepl("show_probs_disabled", html, fixed = TRUE))
    }
  )
})
