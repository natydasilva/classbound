test_that("generate_outlier uses 2D Mahalanobis for positively correlated data", {
  df2 <- data.frame(
    Sim = rep("A", 4),
    X1 = c(1, 2, 3, 4),
    X2 = c(1.1, 1.9, 3.1, 3.9) # Strongly positively correlated, n=4, full rank
  )
  
  # Means: mu_X1 = 2.5, mu_X2 = 2.5
  outlier2 <- generate_outlier(df2, "A", magnitude = 1, index = 1)
  
  mu <- colMeans(df2[, c("X1", "X2")])
  Sigma <- stats::cov(df2[, c("X1", "X2")])
  pt <- c(outlier2$X1, outlier2$X2)
  dist_sq <- as.numeric(t(pt - mu) %*% solve(Sigma) %*% (pt - mu))
  expect_equal(dist_sq, 1^2, tolerance = 1e-5)
})

test_that("generate_outlier falls back to Tukey for n=2", {
  df <- data.frame(
    Sim = c("A", "A"),
    X1 = c(1, 4),
    X2 = c(2, 6)
  )
  
  # For n=2, covariance is singular, Mahalanobis fails. Fallback to Tukey.
  # Q1=1.75, Q3=3.25, IQR=1.5
  outlier <- generate_outlier(df, "A", magnitude = 1, index = 3) # corner = 0 (+max)
  
  expect_equal(outlier$X1, 3.25 + 1.5 * 1)
})

test_that("generate_outlier falls back to deterministic scaling for zero IQR / zero variance", {
  df <- data.frame(
    Sim = c("A", "A", "A", "B", "B"),
    X1 = c(1, 1, 1, 0, 10), # Class A has 0 IQR. Global range is 10.
    X2 = c(2, 3, 4, 1, 5)   # Class A has valid IQR
  )
  
  outlier <- generate_outlier(df, "A", magnitude = 2, index = 3) # +max
  
  # X1: global scale = 10, base = 1. max -> 1 + 2 * 10 = 21
  expect_equal(outlier$X1, 21)
  
  # X2: IQR = 4 - 2 = 2 (Wait, Q1=2.5, Q3=3.5, IQR=1.0). max -> 3.5 + 2 * 1.0 = 5.5
  expect_equal(outlier$X2, 5.5)
})

test_that("generate_outlier handles non-visualized dimensions and categoricals", {
  df <- data.frame(
    Sim = c("A", "A", "A", "B", "B"),
    X1 = c(1, 2, 3, 4, 5), # visual
    X2 = c(10, 20, 30, 40, 50), # visual
    X3 = c(100, 200, 300, 400, 500), # non-visual
    Cat = factor(c("Yes", "Yes", "No", "No", "No"))
  )
  
  outlier <- generate_outlier(df, "A", magnitude = 1, index = 1)
  
  # X3 should be median of target class (200)
  expect_equal(outlier$X3, 200)
  # Cat should be mode of target class ("Yes")
  expect_equal(as.character(outlier$Cat), "Yes")
})

test_that("generate_outlier keeps spatial separation with identical Mahalanobis distance", {
  df <- data.frame(
    Sim = rep("A", 4),
    X1 = c(1, 2, 1, 2),
    X2 = c(1, 1, 2, 2)
  )
  
  pts <- lapply(1:4, function(i) generate_outlier(df, "A", magnitude = 2, index = i))
  
  mu <- c(1.5, 1.5)
  Sigma <- stats::cov(df[, c("X1", "X2")])
  inv_Sigma <- solve(Sigma)
  
  dists <- sapply(pts, function(p) {
    vec <- c(p$X1, p$X2) - mu
    as.numeric(t(vec) %*% inv_Sigma %*% vec)
  })
  
  # All distances squared should be 2^2 = 4
  expect_true(all(abs(dists - 4) < 1e-5))
  
  # They should be spatially separated (no two points identical)
  # We check uniqueness by summing coordinates
  sums <- sapply(pts, function(p) p$X1 + p$X2 * 10)
  expect_length(unique(sums), 4)
})
