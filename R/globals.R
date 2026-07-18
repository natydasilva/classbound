# Declare global variables to prevent R CMD check notes when using ggplot2 non-standard evaluation
utils::globalVariables(c(
  ".",
  ".data",
  "estimate", 
  "label",
  "X1",
  "X2",
  "pred",
  "ppred",
  "Sim"
))
