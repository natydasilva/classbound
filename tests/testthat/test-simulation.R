test_that("simulate_mixsim works correctly", {
  skip_if_not_installed("MixSim")
  
  df <- simulate_mixsim(n = 100, K = 3, p = 2, MaxOmega = 0.1, seed = 42)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 100)
  expect_equal(ncol(df), 3) # Sim + X1 + X2
  expect_equal(colnames(df), c("Sim", "X1", "X2"))
  expect_true(is.factor(df$Sim))
  expect_equal(nlevels(df$Sim), 3)
})

test_that("simulate_mixsim noise_ratio works", {
  skip_if_not_installed("MixSim")
  
  df <- simulate_mixsim(n = 100, K = 3, p = 2, MaxOmega = 0.1, seed = 42, noise_ratio = 0.1)
  expect_equal(nrow(df), 110) # 100 + 10
})

test_that("simulate_mixsim is reproducible", {
  skip_if_not_installed("MixSim")
  
  df1 <- simulate_mixsim(n = 100, K = 3, p = 2, MaxOmega = 0.1, seed = 123)
  df2 <- simulate_mixsim(n = 100, K = 3, p = 2, MaxOmega = 0.1, seed = 123)
  df3 <- simulate_mixsim(n = 100, K = 3, p = 2, MaxOmega = 0.1, seed = 456)
  
  expect_equal(df1, df2)
  expect_false(identical(df1, df3))
})

test_that("simu_n works correctly", {
  skip_if_not_installed("MASS")
  
  means <- list(c(0, 0), c(5, 5))
  covs <- list(diag(2), diag(2))
  ns <- c(50, 50)
  
  df <- simu_n(means, covs, ns, seed = 42)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 100)
  expect_equal(ncol(df), 3)
  expect_equal(colnames(df), c("Sim", "X1", "X2"))
  expect_true(is.factor(df$Sim))
  expect_equal(nlevels(df$Sim), 2)
})

test_that("simu_n noise_ratio works", {
  skip_if_not_installed("MASS")
  
  means <- list(c(0, 0), c(5, 5))
  covs <- list(diag(2), diag(2))
  ns <- c(50, 50)
  
  df <- simu_n(means, covs, ns, seed = 42, noise_ratio = 0.2)
  expect_equal(nrow(df), 120) # 100 + 20
  
  df_zero <- simu_n(means, covs, ns, seed = 42, noise_ratio = 0)
  expect_equal(nrow(df_zero), 100)
})

test_that("simu_n is reproducible", {
  skip_if_not_installed("MASS")
  
  means <- list(c(0, 0), c(5, 5))
  covs <- list(diag(2), diag(2))
  ns <- c(50, 50)
  
  df1 <- simu_n(means, covs, ns, seed = 123)
  df2 <- simu_n(means, covs, ns, seed = 123)
  df3 <- simu_n(means, covs, ns, seed = 456)
  
  expect_equal(df1, df2)
  expect_false(identical(df1, df3))
})

test_that("noise bounds are respected and errors are caught", {
  skip_if_not_installed("MASS")
  
  means <- list(c(0, 0), c(5, 5))
  covs <- list(diag(2), diag(2))
  ns <- c(50, 50)
  
  df_base <- simu_n(means, covs, ns, seed = 42)
  df_noise <- simu_n(means, covs, ns, seed = 42, noise_ratio = 1.0)
  
  # Base points max/mins
  min_x1 <- min(df_base$X1)
  max_x1 <- max(df_base$X1)
  min_x2 <- min(df_base$X2)
  max_x2 <- max(df_base$X2)
  
  # Check if all noise points are strictly within bounds (noise rows are appended at end)
  noise_only <- df_noise[101:200, ]
  
  expect_true(all(noise_only$X1 >= min_x1 & noise_only$X1 <= max_x1))
  expect_true(all(noise_only$X2 >= min_x2 & noise_only$X2 <= max_x2))
  
  # Invalid noise ratios
  expect_error(simu_n(means, covs, ns, noise_ratio = "high"))
})

test_that("global .Random.seed is preserved when using seed argument", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("MixSim")
  
  # Ensure .Random.seed exists
  set.seed(999)
  initial_seed <- globalenv()$.Random.seed
  
  means <- list(c(0, 0), c(5, 5))
  covs <- list(diag(2), diag(2))
  ns <- c(50, 50)
  
  # Calling simu_n with a seed
  invisible(simu_n(means, covs, ns, seed = 42))
  after_simu_n <- globalenv()$.Random.seed
  expect_identical(initial_seed, after_simu_n)
  
  # Calling simulate_mixsim with a seed
  invisible(simulate_mixsim(n = 100, K = 3, p = 2, MaxOmega = 0.1, seed = 123))
  after_mixsim <- globalenv()$.Random.seed
  expect_identical(initial_seed, after_mixsim)
})
