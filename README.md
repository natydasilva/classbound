
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

# 1. Load sample dataset
data(data69_1)

# Convert the target to a factor for classification
data69_1$Y <- as.factor(data69_1$Y)

# 2. Fit a decision tree (rpart)
model <- fit_model(data = data69_1[, c("V1", "V2")], labels = data69_1$Y, method = "rpart")

# 3. Compute the boundary
feature_range <- list(V1 = c(min(data69_1$V1), max(data69_1$V1)),
                      V2 = c(min(data69_1$V2), max(data69_1$V2)))
grid_data <- boundary_compute(model, range = feature_range, resolution = 100)

# 4. Plot the decision boundary
plot_boundary(
  boundary = grid_data, 
  obs_data = data69_1, 
  x_col = "V1", 
  y_col = "V2", 
  true_label = "Y"
)
```

### Interactive Exploration

The package also includes a Shiny app to interactively explore models
and boundaries. You can launch it using:

``` r
explorapp()
```
