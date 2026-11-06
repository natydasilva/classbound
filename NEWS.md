# classbound 1.0.1

## Bug fixes

* `fit_model()` now treats numeric or character class labels as classes; previously
  a numeric response (such as `Y` in `data69_1`) was fitted as a regression.
* `fit_model()` only checks the columns used in the formula, so missing values in
  unrelated columns no longer stop fitting.
* `boundary_compute()` can recompute a multi-model result, and recomputing a
  boundary no longer repeats the `"classbound_boundary"` class.
* `plot_boundary(type = "disagreement")` no longer fails when `obs_data` is supplied.
* `boundary_explore_nd()` rescales features before its k-nearest-neighbour fallback
  for models without class probabilities, as classifly does; previously one
  large-scale feature could dominate and boundary points were close to random. It
  also keeps non-syntactic feature names and accepts character columns.
* `fit_model()` evaluates transformations in the formula (such as `log(x)` or
  `I(x^2)`) with the formula interface; previously the classifier failed with
  "object not found". The features of the model are the original variables.
* `fit_model()` now stops with a clear message when the response has fewer than
  two classes, instead of failing inside the classifier.
* `plot_boundary()` with `obs_data` and `classbound()` accept a projection basis
  without row names (such as one taken directly from `tourr`), assuming the rows
  follow the order of the training features, instead of failing.
* `classbound()` takes the axes from the fitted model, so formulas containing
  transformations work.
* `print()` on a `classbound` object prints the fitted model on separate lines
  instead of a literal `\n`.
* The internal `tidymodels_bridge()` now passes `resolution` on.

## Shiny application

* In Projection mode, the projection plane is recomputed when new data are
  generated; previously the old plane was kept, which broke the plots when the
  number of dimensions changed.
* Imported data keep rows whose missing values are only in columns the app does
  not use, and data with fewer than two numeric features are rejected up front.
* Custom models can declare `supports_prob = TRUE`, and fitted parsnip models from
  the workspace can show a probability surface.
* Clear Canvas also resets the undo history and the number of base classes.
* Exported grid predictions and the reproduce script now use the same plot range
  as the app, and the script passes on the prediction arguments (such as the
  PPtree rule). Each export uses its own temporary folder.
* The export wizard creates its ZIP file with the 'zip' package when it is
  installed, so exporting no longer needs a `zip` program (often missing on Windows).
* The reproduce script writes zoom limits with full precision, and the exported
  `data.csv` has an `is_outlier` column so the script redraws the outliers as
  diamonds, as in the app.
* A zero-size brush no longer causes a division by zero when zooming.
* An empty or invalid Point Density value no longer stops the app while drawing.
* Outliers belonging to a class that has been removed no longer cause ordinary
  observations to be drawn as outliers.
* Only installed models are preselected at start-up.
* The information table reports the 2D Slice mode correctly.

## New features

* `lda` and `qda` models from 'MASS' are supported out of the box, with their
  posterior probabilities available for the probability surface.

## Documentation

* Projection views are now described as what they show: the predicted regions on
  a plane through the data center (a planar section), with the observations
  projected onto that plane.
* `boundary_explore_nd()` credits the 'classifly' package and Wickham, Cook and
  Hofmann (2015).
* The `data69_1` help page now states that `Y` is numeric with values 0, 1, 2 and
  that the data are the full UCI dataset.
* Corrected the return values documented for `predict_adapter()` and
  `boundary_workflow_set()`, and Ignacio Alvarez-Castro's e-mail address.
* The README and the getting-started vignette show installation from CRAN, and
  the README has a CRAN status badge.

## Tests

* Added regression tests for all of the fixes above (the suite now has 96 tests),
  tests for `classbound()`, and tests of the Shiny application's server logic
  with `shiny::testServer()`.

# classbound 1.0.0

* First CRAN release.
