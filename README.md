
<!-- README.md is generated from README.Rmd. Please edit that file -->

# classbound

<!-- badges: start -->

[![R-CMD-check](https://github.com/natydasilva/classbound/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/natydasilva/classbound/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/natydasilva/classbound/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/natydasilva/classbound/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

The **classbound** package provides tools to explore, visualize, and
analyze classification decision boundaries in R. It offers a unified
interface for fitting machine learning models and rendering 2D
visualizations of their prediction surfaces.

## Installation

You can install the development version of classbound from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("natydasilva/classbound")
```

## Quick Start

### Visualizing a Decision Boundary

Here is a basic example of how to fit a model and visualize its decision
boundary:

``` r
library(classbound)
library(palmerpenguins)

# 1. Load and prepare sample dataset
penguins <- na.omit(penguins[, -c(2, 7, 8)])

# 2. Fit a decision tree (rpart) predicting species
model <- fit_model(
  data = penguins[, c("bill_length_mm", "bill_depth_mm")], 
  labels = penguins$species, 
  method = "rpart"
)

# 3. Compute the boundary
feature_range <- list(
  bill_length_mm = c(30.0, 60.0),
  bill_depth_mm = c(10.0, 25.0)
)
grid_data <- boundary_compute(model, range = feature_range, resolution = 100)

# 4. Plot the decision boundary
plot_boundary(
  boundary = grid_data, 
  obs_data = penguins, 
  x_col = "bill_length_mm", 
  y_col = "bill_depth_mm", 
  true_label = "species"
)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" alt="" width="100%" />

### Interactive Exploration

The package also includes a Shiny app to interactively explore models
and boundaries. You can launch it using:

``` r
explorapp()
```
