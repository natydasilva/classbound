peng_2d <- function() {
  as.data.frame(na.omit(palmerpenguins::penguins[, c("species", "bill_length_mm", "bill_depth_mm")]))
}

test_that("print.classbound prints the native model with real line breaks", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("palmerpenguins")

  m <- fit_model(peng_2d(), species ~ ., rpart::rpart)
  out <- utils::capture.output(print(m))

  expect_false(any(grepl("\\n", out, fixed = TRUE)))
  expect_equal(
    out[(which(out == "-- Native Model --") + 1):length(out)],
    utils::capture.output(print(m$fit))
  )
})

test_that("tidymodels_bridge passes `resolution` on to the grid", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("workflowsets")
  skip_if_not_installed("rpart")
  skip_if_not_installed("palmerpenguins")

  res <- tidymodels_bridge(peng_2d(), "species", c("rpart"), resolution = 7)
  expect_equal(nrow(res$boundary_data), 7 * 7)
})

test_that("recomputing a boundary does not duplicate the class attribute", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("palmerpenguins")

  m <- fit_model(peng_2d(), species ~ ., rpart::rpart)
  m <- boundary_compute(m, resolution = 5)
  m <- boundary_compute(m, resolution = 5)
  expect_equal(class(m), c("classbound_boundary", "classbound"))
})

test_that("disagreement plots can overlay observations", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("palmerpenguins")

  d <- peng_2d()
  set.seed(1)
  m1 <- fit_model(d, species ~ ., rpart::rpart)
  m2 <- fit_model(d, species ~ ., randomForest::randomForest, interface = "matrix")
  mb <- boundary_compute(list(tree = m1, rf = m2), resolution = 10)

  p <- plot_boundary(mb,
    obs_data = d, x_col = "bill_length_mm", y_col = "bill_depth_mm",
    true_label = "species", type = "disagreement"
  )
  expect_s3_class(p, "ggplot")
  built <- suppressWarnings(ggplot2::ggplot_build(p))
  expect_equal(nrow(built$data[[2]]), nrow(d))
})

test_that("boundary_explore_nd kNN fallback is not dominated by large-scale features", {
  skip_if_not_installed("e1071")

  set.seed(1)
  n <- 400
  df <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n, 0, 1000))
  df$y <- factor(ifelse(df$x1 > 0.5, "a", "b"))
  # A default e1071::svm fit returns no class probabilities, so the kNN
  # approximation of the advantage is used.
  m <- fit_model(df, y ~ ., e1071::svm)
  expect_null(predict_model(m, df[1:3, ])$probs)

  set.seed(2)
  res <- boundary_explore_nd(list(svm = m), df, n_points = 2000, threshold = 0.02)

  # The true boundary is x1 = 0.5; boundary-near points must concentrate there.
  expect_lt(nrow(res), 200)
  expect_lt(stats::median(abs(res$x1 - 0.5)), 0.05)
  expect_true(all(res$is_boundary))
})

test_that("data69_1 matches its documentation", {
  e <- new.env(parent = emptyenv())
  utils::data("data69_1", package = "classbound", envir = e)
  d <- e$data69_1

  expect_equal(dim(d), c(5000L, 22L))
  expect_equal(sort(unique(d$Y)), c(0, 1, 2))
  expect_true(all(vapply(d[, paste0("V", 1:21)], is.numeric, logical(1))))
})

test_that("fit_model ignores missing values in columns the formula does not use", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("palmerpenguins")

  pg <- as.data.frame(palmerpenguins::penguins)
  pg <- pg[!is.na(pg$bill_length_mm), ] # `sex` still has missing values
  expect_true(anyNA(pg$sex))

  m <- fit_model(pg, species ~ bill_length_mm + bill_depth_mm, rpart::rpart)
  expect_equal(m$metadata$features$names, c("bill_length_mm", "bill_depth_mm"))
  expect_error(
    fit_model(as.data.frame(palmerpenguins::penguins), species ~ bill_length_mm + bill_depth_mm, rpart::rpart),
    "Missing values"
  )
})

test_that("fit_model treats numeric class labels as classes", {
  skip_if_not_installed("rpart")

  e <- new.env(parent = emptyenv())
  utils::data("data69_1", package = "classbound", envir = e)
  d <- as.data.frame(e$data69_1[1:300, c("Y", "V1", "V2")])
  expect_true(is.numeric(d$Y))

  m <- fit_model(d, Y ~ ., rpart::rpart)
  expect_equal(m$fit$method, "class")
  b <- boundary_compute(m, resolution = 5)
  expect_true(is.factor(b$boundary_data$prediction))
  expect_false(anyNA(b$boundary_data$prediction))
  expect_equal(levels(b$boundary_data$prediction), c("0", "1", "2"))
})

test_that("boundary_compute can recompute a multi-model result", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("palmerpenguins")

  m <- fit_model(peng_2d(), species ~ ., rpart::rpart)
  mb <- boundary_compute(list(a = m, b = m), resolution = 5)
  mb2 <- boundary_compute(mb, resolution = 6)
  expect_s3_class(mb2, "classbound_multi")
  expect_equal(nrow(mb2$boundary_data), 2 * 6 * 6)
  expect_equal(unique(mb2$boundary_data$model), c("a", "b"))
})

test_that("boundary_explore_nd keeps feature names that are not syntactic", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("palmerpenguins")

  d <- peng_2d()
  names(d) <- c("species", "bill length", "bill depth")
  m <- fit_model(d, species ~ ., rpart::rpart)
  set.seed(1)
  res <- boundary_explore_nd(list(tree = m), d, n_points = 500)
  expect_true(all(c("bill length", "bill depth", "Y_tree") %in% colnames(res)))
})

test_that("grid export and reproduce script use the same grid as the app plots", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("palmerpenguins")

  d <- peng_2d()
  names(d)[1] <- "Sim"
  m <- fit_model(d, Sim ~ ., rpart::rpart)
  rng <- app_grid_ranges(d)

  f <- tempfile(fileext = ".csv")
  export_grid_csv(m, d, resolution = 10, file = f)
  g <- utils::read.csv(f)
  expect_equal(range(g$x), rng[[1]])
  expect_equal(range(g$y), rng[[2]])

  dir <- tempfile()
  dir.create(dir)
  export_data_csv(d, file.path(dir, "data.csv"))
  export_models_rds(list(tree = m), file.path(dir, "models.rds"))
  export_reproduce_script("tree",
    resolution = 10, predict_args = list(tree = list()),
    file = file.path(dir, "reproduce.R")
  )
  script <- readLines(file.path(dir, "reproduce.R"))
  expect_true(any(grepl("pad_range", script, fixed = TRUE)))
  expect_true(any(grepl("predict_args = get_predict_args(model_name)", script, fixed = TRUE)))

  owd <- setwd(dir)
  on.exit(setwd(owd), add = TRUE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  env <- new.env()
  expect_no_error(utils::capture.output(sys.source("reproduce.R", envir = env)))
  expect_equal(env$range_list[[1]], rng[[1]])
})

test_that("transformations in the formula work with the formula interface", {
  m <- fit_model(iris, Species ~ log(Sepal.Length) + I(Sepal.Width^2), rpart::rpart)
  # The features are the raw variables, so a grid can be passed back through the formula.
  expect_equal(m$metadata$features$names, c("Sepal.Length", "Sepal.Width"))
  b <- boundary_compute(m, resolution = 8)
  expect_equal(nrow(b$boundary_data), 64)
  expect_s3_class(plot_boundary(b), "ggplot")
})

test_that("formula transformations still work with the matrix interface", {
  skip_if_not_installed("randomForest")
  expect_no_error(
    fit_model(iris, Species ~ log(Sepal.Length) + Sepal.Width, randomForest::randomForest,
      interface = "matrix"
    )
  )
})

test_that("fit_model() requires at least two classes", {
  one_class <- iris[iris$Species == "setosa", ]
  expect_error(
    fit_model(one_class, Species ~ Sepal.Length + Sepal.Width, rpart::rpart),
    "at least two classes"
  )
})

test_that("lda and qda are supported without a custom predfun", {
  skip_if_not_installed("MASS")
  for (clf in list(MASS::lda, MASS::qda)) {
    m <- fit_model(iris, Species ~ Sepal.Length + Sepal.Width, clf)
    b <- boundary_compute(m, resolution = 8)
    expect_true(all(levels(iris$Species) %in% colnames(b$boundary_data)))
    expect_equal(unname(rowSums(b$boundary_data[, levels(iris$Species)])), rep(1, 64), tolerance = 1e-8)
  }
})

test_that("a projection basis without row names can be used with obs_data and classbound()", {
  m <- fit_model(iris, Species ~ ., rpart::rpart)
  pc <- stats::prcomp(iris[, 1:4], scale. = TRUE)
  proj <- list(basis = unname(pc$rotation[, 1:2]), center = pc$center, scale = pc$scale)
  b <- boundary_compute(m,
    feature_range = list(PC1 = c(-4, 4), PC2 = c(-3, 3)), resolution = 8,
    projection = proj
  )
  expect_no_error(plot_boundary(b, obs_data = iris, true_label = "Species"))
  expect_no_error(classbound(iris, Species ~ ., rpart::rpart, projection = proj, resolution = 8))
})

test_that("classbound() works with transformations in the formula", {
  expect_s3_class(
    classbound(iris, Species ~ log(Sepal.Length) + Sepal.Width, rpart::rpart, resolution = 8),
    "ggplot"
  )
})

test_that("the export folder can be zipped and contains the files", {
  skip_if_not_installed("zip")
  dir <- tempfile()
  dir.create(file.path(dir, "plots"), recursive = TRUE)
  writeLines("a", file.path(dir, "data.csv"))
  writeLines("b", file.path(dir, "plots", "p.png"))
  z <- tempfile(fileext = ".zip")
  zip_export_dir(z, dir)
  expect_setequal(utils::unzip(z, list = TRUE)$Name, c("data.csv", "plots/p.png"))
})

test_that("flag_outliers() marks the appended rows", {
  d <- data.frame(Sim = c("a", "b", "a", "b"), X1 = 1:4, X2 = 4:1)
  expect_equal(flag_outliers(d, 2)$is_outlier, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(flag_outliers(d, 0)$is_outlier, rep(FALSE, 4))
  expect_equal(flag_outliers(d, 10)$is_outlier, rep(TRUE, 4))
})

test_that("the reproduce script keeps full precision of zoom limits and highlights outliers", {
  f <- tempfile(fileext = ".R")
  export_reproduce_script("tree",
    resolution = 10, zoom_x = c(1.123456789012, 2.2), zoom_y = c(-0.000000123456, 1),
    predict_args = list(tree = list()), file = f
  )
  script <- paste(readLines(f), collapse = "
")
  expect_match(script, "1.123456789012", fixed = TRUE)
  expect_match(script, "-1.23456e-07", fixed = TRUE)
  expect_match(script, 'highlight_outliers = "is_outlier" %in% colnames(dat)', fixed = TRUE)
})

