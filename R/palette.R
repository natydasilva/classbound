#' Generate deterministic colors for classification classes
#'
#' @description
#' Provides a deterministic color palette mapping for a given set of classes. 
#' It uses a highly curated palette of 20 visually distinct colors for the first
#' 20 classes. For any additional classes, it falls back to a mathematical formula 
#' based on the golden angle to generate infinitely many unique, perceptually distinct colors.
#' 
#' @details
#' The mapping is deterministic: classes are first sorted alphabetically before colors are assigned.
#' This ensures that a given class label will always receive the exact same color regardless
#' of the order in which the classes are provided.
#'
#' @param classes A character vector or factor containing class labels.
#' @return A named character vector where the names are the unique, sorted classes and 
#'   the values are the assigned hex color codes.
#' @export
classbound_palette <- function(classes) {
  if (length(classes) == 0) return(stats::setNames(character(0), character(0)))
  
  unique_classes <- sort(unique(as.character(classes)))
  
  curated_colors <- c(
    "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4",
    "#BFEF45", "#F032E6", "#469990", "#9A6324", "#FABED4",
    "#808000", "#2E8B57", "#FF69B4", "#DCBEFF", "#800000",
    "#42D4F4", "#A9A9A9", "#AAFFC3", "#000075", "#FFE119"
  )
  
  get_palette_color <- function(idx) {
    if (idx <= length(curated_colors)) return(curated_colors[idx])
    overflow_idx <- idx - length(curated_colors)
    hue <- (overflow_idx * 0.618033988749895) %% 1
    grDevices::hcl(h = hue * 360, c = 65, l = 55)
  }
  
  cols <- vapply(seq_along(unique_classes), get_palette_color, character(1))
  stats::setNames(cols, unique_classes)
}
