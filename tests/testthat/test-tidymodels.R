test_that("tidymodels adapter works", {
  skip_if_not_installed("tidymodels")

  library(tidymodels)

  # Prepare simple data
  data(data69_1, package = "classbound")
  df <- data69_1[, c("V1", "V2", "Y")]
  df$Y <- as.factor(df$Y)

  # 1. Test parsnip model_fit natively
  spec <- decision_tree(mode = "classification") %>% set_engine("rpart")
  mf_fit <- fit(spec, Y ~ V1 + V2, data = df)

  cb_mf <- as_classbound(mf_fit, data = df, response = "Y")
  expect_s3_class(cb_mf, "classbound")
  expect_equal(cb_mf$metadata$class_levels, c("0", "1", "2"))

  grid_mf_model <- boundary_compute(cb_mf, list(V1 = c(-1, 1), V2 = c(-1, 1)), resolution = 10)
  grid_mf <- grid_mf_model$boundary_data
  expect_true(all(c("x", "y", "prediction") %in% colnames(grid_mf)))

  # 2. Test workflow natively
  rec <- recipe(Y ~ V1 + V2, data = df) %>% step_normalize(all_numeric_predictors())
  wf <- workflow() %>%
    add_recipe(rec) %>%
    add_model(spec)
  wf_fit <- fit(wf, data = df)

  cb_wf <- as_classbound(wf_fit, data = df, response = "Y")
  grid_wf_model <- boundary_compute(cb_wf, list(V1 = c(-1, 1), V2 = c(-1, 1)), resolution = 10)
  grid_wf <- grid_wf_model$boundary_data
  expect_true(all(c("x", "y", "prediction") %in% colnames(grid_wf)))

  # 3. Test workflow_set native helper
  wf_set <- workflow_set(preproc = list(rec = rec), models = list(tree = spec))
  # boundary_workflow_set should auto-fit since this is unfitted
  grid_wfs <- boundary_workflow_set(wf_set, data = df, range = list(V1 = c(-1, 1), V2 = c(-1, 1)), response = "Y", resolution = 10)
  expect_true(all(c("model", "x", "y", "prediction") %in% colnames(grid_wfs)))
  expect_equal(grid_wfs$model[1], "rec_tree")
})
