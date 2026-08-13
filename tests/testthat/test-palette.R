test_that("classbound_palette returns empty character vector for empty input", {
  expect_equal(length(classbound_palette(character(0))), 0)
  expect_equal(length(classbound_palette(NULL)), 0)
})

test_that("classbound_palette is deterministic regardless of input order", {
  classes1 <- c("B", "A", "C")
  classes2 <- c("A", "C", "B")
  
  pal1 <- classbound_palette(classes1)
  pal2 <- classbound_palette(classes2)
  
  expect_equal(names(pal1), c("A", "B", "C"))
  expect_equal(names(pal2), c("A", "B", "C"))
  expect_equal(pal1, pal2)
})

test_that("classbound_palette returns curated colors for <= 20 classes", {
  classes <- paste0("Class", 1:5)
  pal <- classbound_palette(classes)
  
  expect_equal(length(pal), 5)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pal)))
  
  # Ensure the first class alphabetically gets the first curated color
  # Classes alphabetically: Class1, Class2, Class3, Class4, Class5
  expect_equal(pal[["Class1"]], "#E6194B")
})

test_that("classbound_palette uses golden-angle fallback for > 20 classes", {
  classes <- paste0("Class", 1:25)
  pal <- classbound_palette(classes)
  
  expect_equal(length(pal), 25)
  expect_equal(length(unique(pal)), 25) # all colors must be unique
  
  # 21st class (alphabetically) gets the first generated color
  # Note: alphabetical sorting of Class1..Class25: Class1, Class10, Class11...
  sorted_classes <- sort(classes)
  expect_equal(names(pal), sorted_classes)
  
  col_21 <- pal[[21]]
  expect_true(grepl("^#[0-9A-Fa-f]{6}$", col_21))
  expect_true(col_21 != "#E6194B")
})

test_that("classbound_palette produces identical colors on repeated calls", {
  classes <- c("Red", "Green", "Blue")
  pal1 <- classbound_palette(classes)
  pal2 <- classbound_palette(classes)
  expect_equal(pal1, pal2)
})

test_that("plot_boundary respects manual colors override over palette", {
  df <- data.frame(
    x = c(1, 2, 3),
    y = c(1, 2, 3),
    prediction = factor(c("A", "B", "C"))
  )
  
  mod <- list(boundary_data = df)
  class(mod) <- "classbound"
  
  manual_cols <- c("A" = "#111111", "B" = "#222222", "C" = "#333333")
  
  # Manual colors should override palette="Dark2"
  p <- plot_boundary(mod, colors = manual_cols, palette = "Dark2")
  
  # Check if scale_fill_manual is in the plot layers/scales
  scales <- p$scales$scales
  fill_scale <- Find(function(s) "fill" %in% s$aesthetics, scales)
  
  expect_false(is.null(fill_scale))
  expect_equal(fill_scale$palette(1), manual_cols)
})

test_that("plot_boundary warns and falls back to classbound_palette when palette limit exceeded", {
  df <- data.frame(
    x = 1:10,
    y = 1:10,
    prediction = factor(paste0("C", 1:10)) # 10 classes
  )
  
  mod <- list(boundary_data = df)
  class(mod) <- "classbound"
  
  expect_warning(
    p <- plot_boundary(mod, palette = "Dark2"),
    "only supports 8 classes"
  )
  
  scales <- p$scales$scales
  fill_scale <- Find(function(s) "fill" %in% s$aesthetics, scales)
  
  # Should fallback to manual scale (classbound_palette)
  expect_true(inherits(fill_scale, "ScaleDiscrete"))
})

test_that("plot_boundary throws error for invalid palette", {
  df <- data.frame(x = 1, y = 1, prediction = factor("A"))
  mod <- list(boundary_data = df)
  class(mod) <- "classbound"
  
  expect_error(
    plot_boundary(mod, palette = "FakePalette"),
    "Invalid palette name"
  )
})
