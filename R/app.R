#' Shiny app to visually explore and compare classification boundaries for built-in, custom, and tidymodels algorithms in 2D
#'
#' @param data Optional data frame to import directly into the app.
#' @param target_col Optional string specifying the column in `data` that contains the true class labels.
#' @param custom_models A list of custom models to inject into the app's comparison UI. Each element should be a named list containing at least `fn` (the model fitting function). Optionally, it can contain `args` (a list of arguments to pass to `fn`) and `predict_args` (a function returning a list of arguments for prediction).
#' @return No return value, called for side effects. Shinyapp is launched.
#' @export
#' @examples
#' if (interactive()) {
#'   # Launch with default models
#'   explorapp()
#'
#'   # Launch with a custom SVM model
#'   explorapp(custom_models = list(
#'     "SVM" = list(
#'       fn = e1071::svm,
#'       args = list(kernel = "linear")
#'     )
#'   ))
#' }
#'
explorapp <- function(data = NULL, target_col = NULL, custom_models = list()) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' must be installed to use explorapp().", call. = FALSE)
  }
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("Package 'DT' must be installed to use explorapp().", call. = FALSE)
  }

  if (!is.null(data)) {
    if (is.null(target_col) || !(target_col %in% colnames(data))) {
      stop("Please provide a valid 'target_col' that exists in 'data'.", call. = FALSE)
    }

    feat_cols <- setdiff(colnames(data), target_col)
    if (length(feat_cols) < 2) {
      stop("Classbound requires at least 2 feature columns (excluding the target column) for 2D visualization.", call. = FALSE)
    }
  }


  # UI to Package API Mapping (lazy closures to defer namespace loading)
  app_methods <- list(
    "rpart" = list(
      fn = function(...) {
        if (!requireNamespace("rpart", quietly = TRUE)) {
          stop("Package 'rpart' must be installed to use this model.", call. = FALSE)
        }
        rpart::rpart(...)
      },
      args = list(), supports_prob = TRUE,
      fit_args_fn = function(input) {
        if (!requireNamespace("rpart", quietly = TRUE)) {
          return(list())
        }
        list(control = rpart::rpart.control(cp = input$rpart_cp))
      }
    ),
    "randomForest" = list(
      fn = function(...) {
        if (!requireNamespace("randomForest", quietly = TRUE)) {
          stop("Package 'randomForest' must be installed to use this model.", call. = FALSE)
        }
        randomForest::randomForest(...)
      },
      args = list(), supports_prob = TRUE,
      fit_args_fn = function(input) {
        l <- list(ntree = input$rf_ntree)
        if (!is.na(input$rf_mtry) && input$rf_mtry > 0) l$mtry <- input$rf_mtry
        l
      }
    ),
    "PPtreeViz" = list(
      fn = function(...) {
        if (!requireNamespace("PPtreeViz", quietly = TRUE)) {
          stop("Package 'PPtreeViz' must be installed to use this model.", call. = FALSE)
        }
        PPtreeViz::PPTreeclass(...)
      },
      args = list(), supports_prob = FALSE,
      fit_args_fn = function(input) list(PPmethod = input$pp_method)
    ),
    "PPtreeExtclass" = list(
      fn = function(...) {
        if (!requireNamespace("PPtreeExt", quietly = TRUE)) {
          stop("Package 'PPtreeExt' must be installed to use this model.", call. = FALSE)
        }
        PPtreeExt::PPtreeExtclass(...)
      },
      args = list(), supports_prob = FALSE,
      fit_args_fn = function(input) list(PPmethod = input$pp_method, stop = input$stop)
    ),
    "PPtreeExt_split" = list(
      fn = function(...) {
        if (!requireNamespace("PPtreeExt", quietly = TRUE)) {
          stop("Package 'PPtreeExt' must be installed to use this model.", call. = FALSE)
        }
        PPtreeExt::PPtreeExt_split(...)
      },
      args = list(), supports_prob = FALSE,
      fit_args_fn = function(input) list(PPmethod = input$pp_method)
    ),
    "ppforest2" = list(
      fn = function(...) {
        if (!requireNamespace("ppforest2", quietly = TRUE)) {
          stop("Package 'ppforest2' must be installed to use this model.", call. = FALSE)
        }
        ppforest2::pprf(...)
      },
      args = list(), supports_prob = TRUE,
      fit_args_fn = function(input) {
        if (!requireNamespace("ppforest2", quietly = TRUE)) {
          return(list())
        }
        list(size = input$pprf_size, lambda = input$pprf_lambda, vars = ppforest2::vars_all())
      }
    )
  )

  predict_args <- list(
    "rpart"           = function(...) list(),
    "randomForest"    = function(...) list(),
    "PPtreeViz"       = function(ru) list(Rule = ru),
    "PPtreeExtclass"  = function(...) list(),
    "PPtreeExt_split" = function(ru) list(Rule = ru),
    "ppforest2"       = function(...) list()
  )

  # Dynamically add Tidymodels presets
  if (requireNamespace("parsnip", quietly = TRUE)) {
    tm_friendly <- c(
      "rpart" = "Decision Tree (rpart)",
      "randomForest" = "Random Forest",
      "kernlab" = "SVM (kernlab)",
      "nnet" = "Neural Net (nnet)",
      "ppforest2" = "PP Forest (ppforest2)"
    )
    for (m_key in names(.tidymodels_registry)) {
      friendly_name <- paste0("Tidymodels: ", if (m_key %in% names(tm_friendly)) tm_friendly[[m_key]] else m_key)
      app_methods[[friendly_name]] <- list(
        fn = .tidymodels_registry[[m_key]](),
        args = list(),
        supports_prob = TRUE
      )
      predict_args[[friendly_name]] <- function(...) list()
    }
  }

  # Integrate custom models
  custom_uis <- list()
  for (m_name in names(custom_models)) {
    c_mod <- custom_models[[m_name]]
    if (!is.list(c_mod) || is.null(c_mod$fn)) {
      warning("Custom model '", m_name, "' must be a list containing at least 'fn'. Skipping.")
      next
    }

    app_methods[[m_name]] <- list(
      fn = c_mod$fn,
      args = if (!is.null(c_mod$args)) c_mod$args else list(),
      fit_args_fn = c_mod$fit_args
    )

    if (!is.null(c_mod$ui)) {
      custom_uis[[m_name]] <- c_mod$ui
    }

    predict_args[[m_name]] <- if (!is.null(c_mod$predict_args) && is.function(c_mod$predict_args)) {
      c_mod$predict_args
    } else {
      function(...) list()
    }
  }

  # UI
  ui <- shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        .col-sm-4 { height: 100vh; position: sticky; top: 0; overflow-y: auto; }
        .btn { white-space: normal; word-wrap: break-word; max-width: 100%; }
        body.mode-draw .shiny-plot-output {
          cursor: var(--brush-cursor, crosshair) !important;
        }
        body.mode-draw-point .shiny-plot-output .select-area {
          display: none !important;
        }
      ")),
      shiny::tags$script(shiny::HTML("
        var currentMode = 'Navigate';
        var currentBrushRadius = 12;
        var activeStrokeSvg = null;
        var activeStrokePath = null;
        var activeStrokeD = '';
        var activePlotId = null;

        function clearVisualStroke() {
          if (activeStrokeSvg) {
            activeStrokeSvg.remove();
            activeStrokeSvg = null;
          }
          activeStrokePath = null;
          activeStrokeD = '';
        }

        function getBrushRadius() {
          var slider = document.getElementById('brush_spread');
          if (slider && slider.value && !isNaN(parseFloat(slider.value))) {
            return Math.max(5, parseFloat(slider.value) * 400);
          }
          return currentBrushRadius;
        }

        function updateCursor() {
           var r = currentMode === 'Draw Cluster' ? getBrushRadius() : 5;
           var visualR = Math.min(r, 60); // Clamp to max 120px size (browsers drop cursors >128x128)
           var size = Math.round(visualR * 2 + 2);
           var svg = '<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"' + size + '\" height=\"' + size + '\"><circle cx=\"' + (size/2) + '\" cy=\"' + (size/2) + '\" r=\"' + visualR + '\" stroke=\"black\" stroke-width=\"1\" fill=\"rgba(0,0,0,0.1)\" /></svg>';
           var encoded = btoa(svg);
           document.documentElement.style.setProperty('--brush-cursor', 'url(data:image/svg+xml;base64,' + encoded + ') ' + Math.round(size/2) + ' ' + Math.round(size/2) + ', crosshair');
        }

        $(document).on('shiny:inputchanged', function(event) {
          if (event.name === 'interaction_mode') {
            currentMode = event.value;
            if (currentMode === 'Draw Cluster' || currentMode === 'Draw Point') {
              $('body').addClass('mode-draw');
              if (currentMode === 'Draw Cluster') {
                $('body').addClass('mode-draw-cluster');
                $('body').removeClass('mode-draw-point');
              } else {
                $('body').addClass('mode-draw-point');
                $('body').removeClass('mode-draw-cluster');
              }
              updateCursor();
            } else {
              $('body').removeClass('mode-draw mode-draw-cluster mode-draw-point');
              clearVisualStroke();
            }
          } else if (event.name === 'data_mode') {
             if (event.value !== 'Draw Data') {
               $('body').removeClass('mode-draw mode-draw-cluster mode-draw-point');
               currentMode = 'Navigate';
               clearVisualStroke();
             }
          } else if (event.name === 'brush_spread') {
             var r = Math.max(5, event.value * 400);
             currentBrushRadius = r;
             updateCursor();
          }
        });
        document.addEventListener('mousedown', function(e) {
          var plotEl = e.target.closest('.shiny-plot-output');
          if (plotEl && currentMode === 'Draw Cluster' && $('body').hasClass('mode-draw')) {
            e.stopPropagation();

            activePlotId = plotEl.id;

            if (currentMode === 'Draw Cluster') {
              activeStrokeSvg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
              activeStrokeSvg.style.position = 'fixed';
              activeStrokeSvg.style.top = '0';
              activeStrokeSvg.style.left = '0';
              activeStrokeSvg.style.width = '100vw';
              activeStrokeSvg.style.height = '100vh';
              activeStrokeSvg.style.pointerEvents = 'none';
              activeStrokeSvg.style.zIndex = '9999';

              var strokeW = (getBrushRadius() * 2);
              var opacity = '0.15';

              activeStrokePath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
              activeStrokePath.setAttribute('fill', 'none');
              activeStrokePath.setAttribute('stroke', 'rgba(0, 0, 0, ' + opacity + ')');
              activeStrokePath.setAttribute('stroke-width', strokeW.toString());
              activeStrokePath.setAttribute('stroke-linecap', 'round');
              activeStrokePath.setAttribute('stroke-linejoin', 'round');

              activeStrokeD = 'M ' + e.clientX + ' ' + e.clientY;
              activeStrokePath.setAttribute('d', activeStrokeD);
              activeStrokeSvg.appendChild(activeStrokePath);
              document.body.appendChild(activeStrokeSvg);
            }

            var mapped = null;
            var img = $(plotEl).find('img');
            var coordmap = $(plotEl).data('coordmap') || (img.length ? img.data('coordmap') : null);
            if (coordmap && coordmap.panels && coordmap.panels.length > 0) {
               var p = coordmap.panels[0];
               var rect = (img.length ? img[0] : plotEl).getBoundingClientRect();
               var ox = e.clientX - rect.left;
               var oy = e.clientY - rect.top;
               var dx = p.domain.right - p.domain.left;
               var dy = p.domain.top - p.domain.bottom;
               var rx = p.range.right - p.range.left;
               var ry = p.range.top - p.range.bottom;
               if (rx !== 0 && ry !== 0) {
                  mapped = {
                     x: p.domain.left + ((ox - p.range.left) / rx) * dx,
                     y: p.domain.bottom + ((oy - p.range.bottom) / ry) * dy,
                     domain: p.domain,
                     range: p.range
                  };
               }
            }

            Shiny.setInputValue('draw_stroke_start', {
              plot_id: plotEl.id,
              mapped: mapped
            }, {priority: 'event'});
          }
        }, true);

        document.addEventListener('mousemove', function(e) {
          if (currentMode === 'Draw Cluster' && activePlotId && activeStrokePath) {
             var hoveredPlot = e.target.closest('.shiny-plot-output');
             if (hoveredPlot && hoveredPlot.id === activePlotId) {
               activeStrokeD += ' L ' + e.clientX + ' ' + e.clientY;
               activeStrokePath.setAttribute('d', activeStrokeD);
             }
          }
        });

        document.addEventListener('mouseup', function(e) {
          if (currentMode === 'Draw Cluster') {
            Shiny.setInputValue('draw_stroke_end', Math.random());
            clearVisualStroke();
          }
        });

        var wheelTimeout = null;
        var accumulatedDelta = 0;

        document.addEventListener('wheel', function(e) {
          var plotEl = e.target.closest('.shiny-plot-output');
          if (!plotEl) return;

          if (!e.ctrlKey && !e.metaKey) {
            // Not pressing modifier key, allow normal page scroll
            return;
          }

          e.preventDefault(); // Only prevent default if zooming

          accumulatedDelta += e.deltaY;

          if (wheelTimeout) clearTimeout(wheelTimeout);
          wheelTimeout = setTimeout(function() {
            var dir = accumulatedDelta > 0 ? 1 : -1;
            var steps = Math.min(4, Math.max(1, Math.round(Math.abs(accumulatedDelta) / 100)));

            Shiny.setInputValue('plot_wheel', {
              plot_id: plotEl.id,
              direction: dir * steps,
              nonce: Math.random()
            });
            accumulatedDelta = 0;
          }, 100); // Debounce to smooth out the jitter
        }, { passive: false });
      "))
    ),
    shiny::titlePanel("Classbound Exploration & Comparison"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::wellPanel(
          shiny::radioButtons("data_mode", "Data Mode", choices = if (!is.null(data)) c("Import Data", "Simulate Data", "Draw Data") else c("Simulate Data", "Draw Data"), selected = if (!is.null(data)) "Import Data" else "Simulate Data"),
          shiny::conditionalPanel(
            condition = "input.data_mode == 'Draw Data'",
            shiny::radioButtons("interaction_mode", "Interaction Mode", choices = c("Navigate", "Draw Point", "Draw Cluster")),
            shiny::helpText("Tip: Use Ctrl + Mouse Wheel or crossbar to zoom, and double-click to reset the view.")
          ),
          shiny::conditionalPanel(
            condition = "input.data_mode == 'Simulate Data'",
            shiny::hr(),
            shiny::radioButtons("sim_engine", "Simulation Engine", choices = c("Multivariate Normal (MVN)" = "mvn", "MixSim" = "mixsim")),
            shiny::numericInput("sim_seed", "Random Seed (Optional)", value = 5, min = 1, width = "100%"),
            shiny::helpText("Set a seed to reproduce the same randomly generated dataset. Leave blank to generate a new random dataset each time."),
            shiny::div(
              title = "Adds uniformly distributed random points across the feature space to test how well the classifier handles background contamination.",
              shiny::sliderInput("sim_noise", "Background Noise", min = 0, max = 100, value = 0, step = 5, post = "%", width = "100%")
            ),
            shiny::conditionalPanel(
              condition = "input.sim_engine == 'mvn'",
              shiny::numericInput("sim_n_classes", "Number of Classes", value = 3, min = 2, max = 10),
              shiny::helpText("Tip: Class parameters appear dynamically based on the number of classes chosen above."),
              shiny::uiOutput("sim_params_ui"),
              shiny::div(
                style = "display: flex; gap: 10px; margin-top: 15px; margin-bottom: 10px;",
                shiny::actionButton("sim_do", "Generate Data", class = "btn-primary"),
                shiny::actionButton("sim_reset_mvn", "Clear Canvas", class = "btn-danger")
              ),
              shiny::actionButton("clone_to_draw_sim_mvn", "Clone to Draw Canvas", icon = shiny::icon("copy"), class = "btn-warning")
            ),
            shiny::conditionalPanel(
              condition = "input.sim_engine == 'mixsim'",
              shiny::numericInput("sim_mixsim_k", "Number of Classes (K)", value = 3, min = 2, max = 10),
              shiny::numericInput("sim_mixsim_p", "Dimensions (p)", value = 2, min = 2, max = 10),
              shiny::numericInput("sim_mixsim_omega", "Max Overlap (MaxOmega)", value = 0.05, min = 0.01, max = 0.5, step = 0.01),
              shiny::numericInput("sim_mixsim_n", "Sample Size", value = 300, min = 10, step = 100),
              shiny::div(
                style = "display: flex; gap: 10px; margin-top: 15px; margin-bottom: 10px;",
                shiny::actionButton("sim_mixsim_do", "Generate Data", class = "btn-primary"),
                shiny::actionButton("sim_reset_mix", "Clear Canvas", class = "btn-danger")
              ),
              shiny::actionButton("clone_to_draw_sim_mix", "Clone to Draw Canvas", icon = shiny::icon("copy"), class = "btn-warning")
            ),
            shiny::tags$details(
              style = "margin-top: 15px;",
              shiny::tags$summary("Data Preview", style = "display: list-item; font-weight: bold; cursor: pointer;"),
              shiny::div(
                style = "margin-top: 10px;",
                shiny::tags$style(shiny::HTML(".nav-pills > li > a { padding: 4px 10px !important; font-size: 13px !important; }")),
                shiny::tabsetPanel(
                  type = "pills",
                  shiny::tabPanel("Train", shiny::div(style = "margin-top: 15px;", DT::dataTableOutput("sim_data_table_train"))),
                  shiny::tabPanel("Test", shiny::div(style = "margin-top: 15px;", DT::dataTableOutput("sim_data_table_test")))
                )
              )
            )
          ),
          shiny::conditionalPanel(
            condition = "input.data_mode == 'Draw Data'",
            shiny::hr(),
            shiny::p("Click or brush on any plot to add data."),
            shiny::div(
              title = "Auto-generates numbered classes (Class 1, Class 2, etc.). You can also add custom named classes below.",
              shiny::numericInput("draw_total_classes", "Base Classes", value = 3, min = 1, step = 1)
            ),
            shiny::selectInput("draw_class", "Active Class", choices = c("Class 1", "Class 2", "Class 3"), selected = "Class 1"),
            shiny::div(
              style = "display: flex; gap: 10px; align-items: baseline; margin-bottom: 15px;",
              shiny::textInput("new_class_name", label = NULL, placeholder = "New class name...", width = "100%"),
              shiny::actionButton("add_class", "Add Class", class = "btn-info")
            ),
            shiny::conditionalPanel(
              condition = "input.interaction_mode == 'Draw Cluster'",
              shiny::div(
                title = "Controls how many observations are generated at each position along the path.",
                shiny::numericInput("brush_size", "Point Density", value = 5, min = 1)
              ),
              shiny::div(
                title = "Controls the physical size of the brush used to spread observations around the path. The brush remains visually consistent when zooming.",
                shiny::sliderInput("brush_spread", "Brush Size", min = 0.01, max = 0.20, value = 0.03, step = 0.01)
              )
            ),
            shiny::div(
              style = "display: flex; gap: 10px; margin-bottom: 15px; flex-wrap: wrap;",
              shiny::actionButton("undo_draw", "Undo Last", icon = shiny::icon("undo"), class = "btn-default"),
              shiny::actionButton("clear", "Clear Canvas", class = "btn-danger")
            ),
            shiny::tags$details(
              style = "margin-top: 15px;",
              shiny::tags$summary("Data Preview", style = "display: list-item; font-weight: bold; cursor: pointer;"),
              shiny::div(style = "margin-top: 10px;", DT::dataTableOutput("drawn_points_table"))
            )
          ),
          shiny::conditionalPanel(
            condition = "input.data_mode == 'Import Data'",
            shiny::p("Using dataset provided via console."),
            shiny::hr(),
            shiny::actionButton("clone_to_draw_imp", "Clone to Draw Canvas", icon = shiny::icon("copy"), class = "btn-warning"),
            shiny::tags$details(
              style = "margin-top: 15px;",
              shiny::tags$summary("Data Preview", style = "display: list-item; font-weight: bold; cursor: pointer;"),
              shiny::div(style = "margin-top: 10px;", DT::dataTableOutput("import_data_table"))
            )
          ),
        ),
        shiny::uiOutput("tour_panel"),
        shiny::wellPanel(
          shiny::tags$details(
            shiny::tags$summary("Outlier Injection", style = "display: list-item; font-size: 18px; font-weight: 500; cursor: pointer; margin-bottom: 10px;"),
            shiny::helpText("Add outliers to see how different models react to them."),
            shiny::selectInput("outlier_class", "Outlier Class", choices = c("Random", "Class 1", "Class 2", "Class 3")),
            shiny::sliderInput("outlier_magnitude", "Outlier Magnitude", min = 0, max = 10, value = 1.5, step = 0.5),
            shiny::helpText("Controls how far the outlier is placed from its class distribution.", style = "font-size: 0.85em; margin-top: -10px; margin-bottom: 15px;"),
            shiny::numericInput("outlier_count", "Number of Outliers", value = 1, min = 1, max = 20, step = 1),
            shiny::checkboxInput("highlight_outliers", "Highlight Outliers (Diamonds)", value = TRUE),
            shiny::div(
              style = "display: flex; gap: 10px; margin-bottom: 10px;",
              shiny::actionButton("inject_outlier_btn", "Inject Outliers", class = "btn-danger"),
              shiny::actionButton("clear_outliers_btn", "Clear", class = "btn-default")
            ),
            shiny::htmlOutput("outlier_status_ui")
          )
        ),
        shiny::wellPanel(
          shiny::tags$details(
            shiny::tags$summary("Import Workspace Models", style = "display: list-item; font-size: 18px; font-weight: 500; cursor: pointer; margin-bottom: 10px;"),
            shiny::helpText("This feature imports `workflow`, `model_fit`, or `model_spec` objects from your R Global Environment."),
            shiny::uiOutput("workspace_import_ui")
          )
        ),
        shiny::wellPanel(
          shiny::tags$details(
            shiny::tags$summary("Model Configuration", style = "display: list-item; font-size: 18px; font-weight: 500; cursor: pointer; margin-bottom: 10px;"),
            shiny::checkboxGroupInput(
              "selected_models",
              "Models to Compare",
              choices = names(app_methods),
              selected = c("randomForest", "ppforest2")
            ),
            shiny::hr(),
            shiny::sliderInput("grid_resolution", "Grid Resolution", min = 50, max = 300, value = 100, step = 25),
            shiny::conditionalPanel(
              condition = "input.selected_models && (input.selected_models.indexOf('PPtreeViz') > -1 || input.selected_models.indexOf('PPtreeExtclass') > -1 || input.selected_models.indexOf('PPtreeExt_split') > -1)",
              shiny::hr(),
              shiny::div(
                title = "Defines the mathematical projection index used to separate the classes. 1 = LDA, 2 = PDA, 3 = Lr (etc.).",
                shiny::selectInput("rule", "PPtree: Projection Pursuit Rule", choices = 1:8, selected = 1)
              ),
              shiny::div(
                title = "The projection index strategy.",
                shiny::selectInput("pp_method", "PPtree: PPmethod", choices = c("LDA", "PDA"), selected = "LDA")
              )
            ),
            shiny::conditionalPanel(
              condition = "input.selected_models && input.selected_models.indexOf('PPtreeExtclass') > -1",
              shiny::div(
                title = "Controls when the tree stops splitting. Higher values result in smaller, simpler trees.",
                shiny::numericInput("stop", "PPtreeExt: Stopping Rule", value = 4, min = 1)
              )
            ),
            shiny::conditionalPanel(
              condition = "input.selected_models && input.selected_models.indexOf('rpart') > -1",
              shiny::hr(),
              shiny::div(
                title = "Complexity parameter. Lower values (e.g. 0.001) produce larger, more complex decision trees that risk overfitting.",
                shiny::numericInput("rpart_cp", "rpart: Complexity Parameter", value = 0.01, min = 0, step = 0.01)
              )
            ),
            shiny::conditionalPanel(
              condition = "input.selected_models && input.selected_models.indexOf('randomForest') > -1",
              shiny::hr(),
              shiny::div(
                title = "Total number of decision trees to grow. Higher numbers increase accuracy but take longer to compute.",
                shiny::numericInput("rf_ntree", "Random Forest: Number of Trees", value = 500, min = 10, step = 50)
              ),
              shiny::div(
                title = "Number of variables randomly sampled as candidates at each split. Leave blank for default (sqrt of total features).",
                shiny::numericInput("rf_mtry", "Random Forest: mtry", value = NA, min = 1, step = 1)
              )
            ),
            shiny::conditionalPanel(
              condition = "input.selected_models && input.selected_models.indexOf('ppforest2') > -1",
              shiny::hr(),
              shiny::div(
                title = "Total number of projection pursuit trees to grow. Higher numbers increase accuracy but take longer to compute.",
                shiny::numericInput("pprf_size", "ppforest2: Number of Trees", value = 100, min = 10, step = 10)
              ),
              shiny::div(
                title = "Penalty parameter (lambda) for the Projection Pursuit PDA index.",
                shiny::sliderInput("pprf_lambda", "ppforest2: PDA Lambda", min = 0, max = 1, value = 0.5, step = 0.1)
              )
            ),
            lapply(names(custom_uis), function(m_name) {
              shiny::conditionalPanel(
                condition = sprintf("input.selected_models && input.selected_models.indexOf('%s') > -1", m_name),
                shiny::hr(),
                custom_uis[[m_name]]
              )
            })
          )
        ),
        shiny::wellPanel(
          shiny::tags$details(
            shiny::tags$summary("Visual Settings", style = "display: list-item; font-size: 18px; font-weight: 500; cursor: pointer; margin-bottom: 10px;"),
            shiny::uiOutput("prob_surface_ui"),
            shiny::selectInput("color_palette", "Color Palette", choices = c("classbound Default", "Dark2 (Colorblind)", "Set1"))
          )
        ),
        shiny::wellPanel(
          shiny::h4("Export Results", style = "margin-top: 0;"),
          shiny::actionButton("open_export_wizard", "Export Results...", class = "btn-success", icon = shiny::icon("download", lib = "font-awesome")),
          shiny::helpText(shiny::em("Customize export formats (PNG, PDF, etc.)"))
        ),
        shiny::wellPanel(
          style = "padding-bottom: 5px;",
          shiny::tags$details(
            shiny::tags$summary(
              "\u2753 Help & Guide",
              style = "display: list-item; font-size: 18px; font-weight: 500; cursor: pointer;"
            ),
            shiny::tags$div(
              style = "margin-top: 12px;",
              shiny::tags$input(
                id = "help_search",
                type = "text",
                placeholder = "Search help...",
                style = "width: 100%; padding: 5px 8px; border: 1px solid #ccc; border-radius: 4px; margin-bottom: 10px; font-size: 13px;"
              ),
              shiny::tags$script(shiny::HTML("
                document.addEventListener('input', function(e) {
                  if (e.target && e.target.id === 'help_search') {
                    var q = e.target.value.toLowerCase();
                    document.querySelectorAll('.help-section').forEach(function(sec) {
                      var text = sec.innerText.toLowerCase();
                      sec.style.display = (!q || text.includes(q)) ? '' : 'none';
                    });
                  }
                });
              ")),
              shiny::tags$div(
                class = "help-section",
                shiny::tags$details(
                  shiny::tags$summary("Getting Started", style = "display: list-item; cursor: pointer; font-weight: bold; padding: 4px 0;"),
                  shiny::tags$p("Classbound lets you fit a classifier, compute its decision boundary, and visualize how it partitions the feature space. Use the sidebar to choose a data source, select models, and adjust settings. Plots update automatically.",
                    style = "font-size: 12px; margin: 6px 0;"
                  ),
                  shiny::actionLink("help_getting_started", "Full documentation \u2192", style = "font-size: 12px;")
                )
              ),
              shiny::tags$div(
                class = "help-section",
                shiny::tags$details(
                  shiny::tags$summary("Data Modes", style = "display: list-item; cursor: pointer; font-weight: bold; padding: 4px 0;"),
                  shiny::tags$p(shiny::tags$b("Import Data"), ": uses the dataset you passed to explorapp(). Read-only; cannot be edited in the app. Clone to Draw Canvas to make edits.",
                    style = "font-size: 12px; margin: 6px 0;"
                  ),
                  shiny::tags$p(shiny::tags$b("Simulate Data"), ": generates synthetic data using Multivariate Normal (MVN) or MixSim engines. A separate independent test set is also generated - it is not a split of the training data.",
                    style = "font-size: 12px; margin: 6px 0;"
                  ),
                  shiny::tags$p(shiny::tags$b("Draw Data"), ": click or brush directly on the plot to add observations. Use 'Undo Last' or 'Clear Canvas' to manage drawn points.",
                    style = "font-size: 12px; margin: 6px 0;"
                  ),
                  shiny::actionLink("help_explorapp_guide", "Explorapp guide \u2192", style = "font-size: 12px;")
                )
              ),
              shiny::tags$div(
                class = "help-section",
                shiny::tags$details(
                  shiny::tags$summary("Interaction Modes", style = "display: list-item; cursor: pointer; font-weight: bold; padding: 4px 0;"),
                  shiny::tags$p(shiny::tags$b("Navigate"), ": zoom with Ctrl+wheel, brush to select region, double-click to reset.",
                    style = "font-size: 12px; margin: 6px 0;"
                  ),
                  shiny::tags$p(shiny::tags$b("Draw Point"), ": each click adds one observation at that location.",
                    style = "font-size: 12px; margin: 6px 0;"
                  ),
                  shiny::tags$p(shiny::tags$b("Draw Cluster"), ": drag to paint clusters of observations. Point Density and Brush Size control how many points are added.",
                    style = "font-size: 12px; margin: 6px 0;"
                  )
                )
              ),
              shiny::tags$div(
                class = "help-section",
                shiny::tags$details(
                  shiny::tags$summary("2D Slice vs Projection", style = "display: list-item; cursor: pointer; font-weight: bold; padding: 4px 0;"),
                  shiny::tags$p(shiny::tags$b("2D Slice"), ": when the model uses more than two features, you choose two for the axes. All other features are fixed at their training-set median.",
                    style = "font-size: 12px; margin: 6px 0;"
                  ),
                  shiny::tags$p(shiny::tags$b("Projection"), ": Projects the data onto the first two principal components (PC1 and PC2) using PCA. Training points are overlaid with depth fading, where more transparent points are farther from the projection plane.",
                    style = "font-size: 12px; margin: 6px 0;"
                  ),
                  shiny::actionLink("help_high_dimensional", "High-dimensional guide \u2192", style = "font-size: 12px;")
                )
              ),
              shiny::tags$div(
                class = "help-section",
                shiny::tags$details(
                  shiny::tags$summary("Probability Surface", style = "display: list-item; cursor: pointer; font-weight: bold; padding: 4px 0;"),
                  shiny::tags$p("The probability surface shades decision regions by model confidence: deep colors = high confidence, faded colors = uncertain near boundaries. Only available for classifiers that provide class probabilities (e.g., rpart, randomForest). SVMs and PPtree models always show flat regions.",
                    style = "font-size: 12px; margin: 6px 0;"
                  )
                )
              ),
              shiny::tags$div(
                class = "help-section",
                shiny::tags$details(
                  shiny::tags$summary("Outlier Injection", style = "display: list-item; cursor: pointer; font-weight: bold; padding: 4px 0;"),
                  shiny::tags$p("Injected outliers become part of the training data and affect the fitted boundary. Outlier Magnitude controls how far they are placed from the class distribution. Use Clear to remove outliers without clearing the rest of the dataset.",
                    style = "font-size: 12px; margin: 6px 0;"
                  )
                )
              ),
              shiny::tags$div(
                class = "help-section",
                shiny::tags$details(
                  shiny::tags$summary("Performance Metrics", style = "display: list-item; cursor: pointer; font-weight: bold; padding: 4px 0;"),
                  shiny::tags$p("Training metrics (accuracy, Kappa, error) are computed on the training data. Test Error uses an independently generated test set (not a split of the training data). Test Error is only available in Simulate Data mode.",
                    style = "font-size: 12px; margin: 6px 0;"
                  )
                )
              ),
              shiny::tags$div(
                class = "help-section",
                shiny::tags$details(
                  shiny::tags$summary("Export", style = "display: list-item; cursor: pointer; font-weight: bold; padding: 4px 0;"),
                  shiny::tags$p("The Export Wizard lets you download data (CSV), fitted models (RDS), plots (PNG/PDF), and a reproduce script that regenerates the plots from exported files.",
                    style = "font-size: 12px; margin: 6px 0;"
                  ),
                  shiny::actionLink("help_export", "Export details \u2192", style = "font-size: 12px;")
                )
              )
            )
          )
        )
      ),
      shiny::mainPanel(
        shiny::uiOutput("data_stats_ui"),
        shiny::uiOutput("plot_grid"),
        shiny::conditionalPanel(
          condition = "input.selected_models && input.selected_models.length > 0",
          shiny::hr(),
          shiny::h4("Training Performance Metrics"),
          DT::dataTableOutput("metrics_table"),
          shiny::hr(),
          shiny::h4("Visualization Info"),
          DT::dataTableOutput("vis_info_ui")
        )
      )
    )
  )

  # Server
  server <- function(input, output, session) {
    # --- External Help Links ---
    shiny::observeEvent(input$help_getting_started, {
      utils::browseURL("https://natydasilva.github.io/classbound/articles/getting-started.html")
    })
    shiny::observeEvent(input$help_explorapp_guide, {
      utils::browseURL("https://natydasilva.github.io/classbound/articles/explorapp-guide.html")
    })
    shiny::observeEvent(input$help_high_dimensional, {
      utils::browseURL("https://natydasilva.github.io/classbound/articles/high-dimensional.html")
    })
    shiny::observeEvent(input$help_export, {
      utils::browseURL("https://natydasilva.github.io/classbound/articles/explorapp-guide.html#export")
    })

    # Initialize current_data based on passed 'data'
    initial_imp_data <- data.frame(Sim = character(), X1 = numeric(), X2 = numeric())
    initial_imp_classes <- c("Class 1", "Class 2", "Class 3")

    if (!is.null(data)) {
      initial_imp_data <- stats::na.omit(data)
      # Standardize target column to "Sim" for internal app logic
      colnames(initial_imp_data)[colnames(initial_imp_data) == target_col] <- "Sim"

      # Keep only numeric features and the target
      numeric_cols <- sapply(initial_imp_data, is.numeric)
      numeric_cols[which(colnames(initial_imp_data) == "Sim")] <- TRUE
      initial_imp_data <- initial_imp_data[, numeric_cols, drop = FALSE]

      if (is.factor(initial_imp_data$Sim)) {
        initial_imp_classes <- levels(initial_imp_data$Sim)
      } else {
        natural_sort <- function(x) {
          nums <- suppressWarnings(as.numeric(gsub(".*?(\\d+)$", "\\1", x)))
          if (all(!is.na(nums))) {
            prefix <- gsub("\\d+$", "", x)
            return(x[order(prefix, nums)])
          }
          sort(x)
        }
        initial_imp_classes <- natural_sort(unique(as.character(initial_imp_data$Sim)))
      }
      if (length(initial_imp_classes) == 0) {
        initial_imp_classes <- c("Class 1", "Class 2", "Class 3")
      }
    }

    initial_drawn_data <- data.frame(Sim = character(), X1 = numeric(), X2 = numeric())
    initial_drawn_classes <- c("Class 1", "Class 2", "Class 3")

    # Start with empty simulated data; the observer below will auto-fetch UI defaults.
    initial_sim_data <- data.frame(Sim = character(), X1 = numeric(), X2 = numeric())
    initial_sim_classes <- character()

    # Auto-fetch simulation data on first view of the Simulate Data tab
    auto_fetch_obs <- shiny::observe(suspended = TRUE, {
      shiny::req(input$data_mode)
      # Do not run simulation math if the user is importing data or drawing data
      if (input$data_mode != "Simulate Data") {
        return()
      }

      shiny::req(input$sim_n_classes)
      shiny::req(input$sim_mean_1, input$sim_sd_1, input$sim_cor_1, input$sim_n_1)

      n <- input$sim_n_classes
      if (is.na(n) || n < 1) {
        return()
      }

      means <- list()
      covs <- list()
      ns <- c()

      for (i in 1:n) {
        m_str <- input[[paste0("sim_mean_", i)]]
        s_str <- input[[paste0("sim_sd_", i)]]
        if (is.null(m_str) || is.null(s_str)) {
          return()
        }

        m <- suppressWarnings(as.numeric(unlist(strsplit(m_str, ","))))
        sds <- suppressWarnings(as.numeric(unlist(strsplit(s_str, ","))))
        rho <- input[[paste0("sim_cor_", i)]]
        ns <- c(ns, input[[paste0("sim_n_", i)]])

        means[[i]] <- m
        dim_m <- length(m)
        if (length(sds) != dim_m) {
          return()
        }

        cov_mat <- matrix(rho, nrow = dim_m, ncol = dim_m)
        diag(cov_mat) <- 1
        sd_diag <- diag(sds, dim_m)
        cov_mat <- sd_diag %*% cov_mat %*% sd_diag
        covs[[i]] <- cov_mat
      }

      lengths <- sapply(means, length)
      if (length(unique(lengths)) > 1) {
        return()
      }
      if (any(is.na(unlist(means))) || any(is.na(unlist(covs)))) {
        return()
      }

      seed_val <- if (is.numeric(input$sim_seed) && !is.na(input$sim_seed)) as.integer(input$sim_seed) else NULL
      noise_val <- if (is.numeric(input$sim_noise)) input$sim_noise / 100 else 0
      new_data <- tryCatch(
        {
          simu_n(means = means, covs = covs, ns = ns, seed = seed_val, noise_ratio = noise_val, test_ratio = 0.3)
        },
        error = function(e) {
          NULL
        }
      )

      if (!is.null(new_data)) {
        train_dat <- if (is.list(new_data) && !is.data.frame(new_data)) new_data$train else new_data
        new_classes <- unique(as.character(train_dat$Sim))

        conf <- list(
          engine = "mvn",
          seed = seed_val,
          noise_ratio = noise_val,
          n_classes = n
        )
        for (i in 1:n) {
          conf[[paste0("class_", i)]] <- list(
            mean = input[[paste0("sim_mean_", i)]],
            sd = input[[paste0("sim_sd_", i)]],
            cor = input[[paste0("sim_cor_", i)]],
            n = input[[paste0("sim_n_", i)]]
          )
        }

        # Store in the background state
        if (is.list(new_data) && !is.data.frame(new_data)) {
          mode_states[["Simulate Data"]]$data <- new_data$train
          mode_states[["Simulate Data"]]$test_data <- new_data$test
        } else {
          mode_states[["Simulate Data"]]$data <- new_data
          mode_states[["Simulate Data"]]$test_data <- NULL
        }
        mode_states[["Simulate Data"]]$classes <- new_classes
        mode_states[["Simulate Data"]]$sim_config <- conf

        # Only overwrite the active canvas if the user is actually on the Simulate tab
        if (is.null(input$data_mode) || input$data_mode == "Simulate Data") {
          current_data(mode_states[["Simulate Data"]]$data)
          current_test_data(mode_states[["Simulate Data"]]$test_data)
          injected_outliers(data.frame())
          class_choices(new_classes)
          applied_sim_config(conf)
        }
      }

      # Destroy the observer permanently so it never runs again
      auto_fetch_obs$destroy()
    })
    auto_fetch_obs$resume()

    injected_outliers <- shiny::reactiveVal(data.frame())
    outlier_last_action <- shiny::reactiveVal(list(action = "None", coords = ""))
    undo_history <- shiny::reactiveVal(list())

    # Freehand drawing state
    active_stroke_data <- shiny::reactiveVal(NULL)
    draw_stroke_plot_id <- shiny::reactiveVal(NULL)
    last_sample_pos <- shiny::reactiveVal(NULL)

    combined_training_data <- shiny::reactive({
      dat <- current_data()
      inj <- injected_outliers()
      if (!is.null(inj) && nrow(inj) > 0) {
        dat <- rbind(dat, inj)
      }

      # Filter out inactive classes (e.g. if user reduced Base Classes)
      active <- class_choices()
      if (length(active) > 0 && "Sim" %in% colnames(dat)) {
        dat <- dat[dat$Sim %in% active, , drop = FALSE]
      }

      dat
    })


    mode_states <- shiny::reactiveValues(
      "Import Data" = list(data = if (!is.null(data)) initial_imp_data else NULL, classes = if (!is.null(data)) initial_imp_classes else NULL, outliers = data.frame(), sim_config = NULL),
      "Draw Data" = list(data = initial_drawn_data, classes = initial_drawn_classes, outliers = data.frame(), sim_config = NULL),
      "Simulate Data" = list(data = initial_sim_data, classes = initial_sim_classes, outliers = data.frame(), sim_config = NULL)
    )

    init_mode <- if (!is.null(data)) "Import Data" else "Simulate Data"
    previous_data_mode <- shiny::reactiveVal(init_mode)

    init_data <- if (!is.null(data)) initial_imp_data else initial_sim_data
    init_classes <- if (!is.null(data)) initial_imp_classes else initial_sim_classes

    current_data <- shiny::reactiveVal(init_data)
    current_test_data <- shiny::reactiveVal(NULL)
    applied_sim_config <- shiny::reactiveVal(NULL)
    class_choices <- shiny::reactiveVal(init_classes)
    zoom_xlim <- shiny::reactiveVal(NULL)
    zoom_ylim <- shiny::reactiveVal(NULL)
    ws_update_trigger <- shiny::reactiveVal(0)

    # Strip ANSI escape codes to prevent UI rendering artifacts.
    clean_err_msg <- function(msg) {
      gsub("\033\\[[0-9;]*m", "", msg)
    }

    # Map class names to colors.
    color_palette <- shiny::reactive({
      levs <- class_choices()
      if (length(levs) == 0) {
        return(stats::setNames(character(0), character(0)))
      }

      pal_choice <- if (is.null(input$color_palette)) "classbound Default" else input$color_palette

      if (pal_choice == "Dark2 (Colorblind)") {
        if (length(levs) <= 8) {
          cols <- RColorBrewer::brewer.pal(max(3, length(levs)), "Dark2")[seq_along(levs)]
          return(stats::setNames(cols, levs))
        }
      } else if (pal_choice == "Set1") {
        if (length(levs) <= 9) {
          cols <- RColorBrewer::brewer.pal(max(3, length(levs)), "Set1")[seq_along(levs)]
          return(stats::setNames(cols, levs))
        }
      }

      classbound_palette(levs)
    })

    shiny::observe({
      levs <- class_choices()
      pal_choice <- input$color_palette
      if (is.null(pal_choice)) {
        return()
      }

      if (pal_choice == "Dark2 (Colorblind)" && length(levs) > 8) {
        shiny::updateSelectInput(shiny::getDefaultReactiveDomain(), "color_palette", selected = "classbound Default")
        shiny::showNotification("Dark2 only supports 8 classes. Reverted to classbound Default.", type = "message")
      } else if (pal_choice == "Set1" && length(levs) > 9) {
        shiny::updateSelectInput(shiny::getDefaultReactiveDomain(), "color_palette", selected = "classbound Default")
        shiny::showNotification("Set1 only supports 9 classes. Reverted to classbound Default.", type = "message")
      }
    })

    shiny::observe({
      choices <- class_choices()
      if (length(choices) > 0) {
        choices <- c("Random", choices)
      }
      selected <- shiny::isolate(input$outlier_class)
      if (is.null(selected) || !(selected %in% choices)) selected <- "Random"
      shiny::updateSelectInput(shiny::getDefaultReactiveDomain(), "outlier_class", choices = choices, selected = selected)
    })

    shiny::observeEvent(input$inject_outlier_btn, {
      cd <- current_data()
      if (nrow(cd) == 0) {
        shiny::showNotification("Please generate or load data first before injecting outliers.", type = "warning")
        return()
      }

      inj <- injected_outliers()
      start_idx <- nrow(inj) + 1

      count <- if (is.numeric(input$outlier_count) && input$outlier_count > 0) floor(input$outlier_count) else 1
      target_class_orig <- input$outlier_class

      new_pts <- do.call(rbind, lapply(1:count, function(i) {
        target_class <- target_class_orig
        if (target_class == "Random") {
          avail_classes <- class_choices()
          if (length(avail_classes) > 0) {
            target_class <- sample(avail_classes, 1)
          }
        }

        # Pass the strict statistical magnitude without jittering it.
        # The corner logic inside the function uses `total_index` to separate points spatially
        # (e.g. by picking different corners) without altering their statistical distance.
        total_index <- start_idx + i - 1
        mag <- input$outlier_magnitude
        generate_outlier(cd, target_class, mag, target_col = "Sim", index = total_index)
      }))

      if (nrow(inj) == 0) {
        injected_outliers(new_pts)
      } else {
        injected_outliers(rbind(inj, new_pts))
      }

      feat_cols <- setdiff(colnames(cd), "Sim")
      num_cols <- names(which(sapply(cd[feat_cols], is.numeric)))
      vis_cols <- if (length(num_cols) >= 2) num_cols[1:2] else num_cols

      coord_text <- ""
      if (length(vis_cols) > 0) {
        if (count == 1) {
          coord_str <- paste(sapply(vis_cols, function(col) sprintf("%s = %.2f", col, new_pts[[col]][1])), collapse = ", ")
          coord_str <- paste0(coord_str, sprintf(" [%s]", new_pts[["Sim"]][1]))
          coord_text <- sprintf("<strong>Coordinates:</strong> %s", coord_str)
        } else {
          range_str <- paste(sapply(vis_cols, function(col) {
            rng <- range(new_pts[[col]], na.rm = TRUE)
            sprintf("%s = [%.2f, %.2f]", col, rng[1], rng[2])
          }), collapse = ", ")

          pts_html <- paste(sapply(1:nrow(new_pts), function(r) {
            row_vals <- paste(sapply(vis_cols, function(c) sprintf("%.2f", new_pts[[c]][r])), collapse = ", ")
            sprintf("<li>Pt %d: (%s) [%s]</li>", r, row_vals, new_pts[["Sim"]][r])
          }), collapse = "")

          coord_text <- sprintf(
            "<strong>Coordinate Range:</strong> %s<br/>
            <details style='margin-top: 5px;'>
              <summary style='cursor: pointer; outline: none; display: list-item;'>Show Exact Coordinates</summary>
              <ul style='margin-top: 5px; padding-left: 20px; max-height: 120px; overflow-y: auto; margin-bottom: 0;'>
                %s
              </ul>
            </details>",
            range_str, pts_html
          )
        }
      }

      outlier_last_action(list(
        action = sprintf("Added %d %s outlier(s)", count, target_class_orig),
        coords = coord_text
      ))
    })

    shiny::observeEvent(input$clear_outliers_btn, {
      injected_outliers(data.frame())
      outlier_last_action(list(action = "Cleared outliers", coords = ""))
    })

    output$outlier_status_ui <- shiny::renderUI({
      inj <- injected_outliers()
      last <- outlier_last_action()
      n <- nrow(inj)

      coords_html <- if (is.list(last) && !is.null(last$coords) && last$coords != "") {
        paste0("<br/>", last$coords)
      } else {
        ""
      }
      action_text <- if (is.list(last)) last$action else last

      shiny::HTML(sprintf(
        "<div style='margin-top: 5px; font-size: 0.9em; color: #555;'>
          <strong>Outliers injected:</strong> %d<br/>
          <strong>Last action:</strong> %s%s
        </div>",
        n, action_text, coords_html
      ))
    })

    output$data_stats_ui <- shiny::renderUI({
      dat <- combined_training_data()
      if (is.null(dat) || nrow(dat) == 0) {
        return(NULL)
      }

      n_obs <- nrow(dat)
      dims <- ncol(dat) - 1
      # Filter to active classes.
      # Preserve factor ordering.
      pal <- color_palette()
      present_classes <- unique(as.character(dat$Sim))
      classes <- intersect(class_choices(), present_classes)

      legend_items <- lapply(classes, function(cls) {
        col <- if (cls %in% names(pal)) pal[[cls]] else "#999999"
        cls_count <- sum(as.character(dat$Sim) == cls)
        shiny::span(
          style = "display: inline-flex; align-items: center; margin-left: 10px;",
          shiny::tags$span(style = sprintf("display: inline-block; width: 12px; height: 12px; border-radius: 50%%; margin-right: 5px; background-color: %s;", col)),
          sprintf("%s (%d)", cls, cls_count)
        )
      })

      shiny::div(
        style = "background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; padding: 15px; margin-bottom: 20px; display: flex; justify-content: space-around; align-items: center; flex-wrap: wrap; gap: 15px;",
        shiny::div(shiny::tags$b("N Observations: "), n_obs),
        shiny::div(shiny::tags$b("Dimensions: "), dims),
        shiny::div(
          style = "display: flex; align-items: center; flex-wrap: wrap;",
          shiny::tags$b("Target Classes: "),
          shiny::div(style = "display: inline-flex; flex-wrap: wrap; gap: 5px;", legend_items)
        )
      )
    })

    output$is_high_dim <- shiny::reactive({
      ncol(combined_training_data()) > 3
    })
    shiny::outputOptions(output, "is_high_dim", suspendWhenHidden = FALSE)

    output$prob_surface_ui <- shiny::renderUI({
      req_models <- input$selected_models
      if (length(req_models) == 0) {
        return(NULL)
      }

      # Check if ANY selected model supports probability
      any_support_prob <- any(vapply(req_models, function(m) {
        if (m %in% names(app_methods)) {
          isTRUE(app_methods[[m]]$supports_prob)
        } else if (m %in% names(custom_models)) {
          isTRUE(custom_models[[m]]$supports_prob)
        } else {
          FALSE
        }
      }, logical(1)))

      current_val <- shiny::isolate(isTRUE(input$show_probs))

      shiny::div(
        style = "margin-top: 15px; margin-bottom: 10px;",
        if (any_support_prob) {
          shiny::checkboxInput("show_probs", "Show Probability Surface", value = current_val)
        } else {
          shiny::div(
            title = "Probability surface unavailable: the selected model does not provide class probabilities.",
            style = "opacity: 0.6;",
            shiny::tags$fieldset(
              disabled = NA,
              shiny::checkboxInput("show_probs_disabled", "Show Probability Surface", value = FALSE)
            )
          )
        }
      )
    })

    shiny::observeEvent(input$data_mode, priority = 1, {
      zoom_xlim(NULL)
      zoom_ylim(NULL)

      old_mode <- previous_data_mode()
      new_mode <- input$data_mode

      # Save state to old mode
      if (!is.null(old_mode)) {
        mode_states[[old_mode]]$data <- current_data()
        mode_states[[old_mode]]$test_data <- current_test_data()
        mode_states[[old_mode]]$classes <- class_choices()
        mode_states[[old_mode]]$basis <- current_basis()
        mode_states[[old_mode]]$projection <- current_projection()
        mode_states[[old_mode]]$projection_info <- current_projection_info()
        mode_states[[old_mode]]$outliers <- injected_outliers()
        mode_states[[old_mode]]$tour_var <- input$tour_var
        mode_states[[old_mode]]$tour_angle <- input$tour_angle
        mode_states[[old_mode]]$tour_path <- current_path()
        mode_states[[old_mode]]$sim_config <- applied_sim_config()
      }

      # Load state from new mode
      if (!is.null(mode_states[[new_mode]]$data)) {
        current_data(mode_states[[new_mode]]$data)
        current_test_data(mode_states[[new_mode]]$test_data)

        saved_outliers <- mode_states[[new_mode]]$outliers
        if (is.null(saved_outliers)) saved_outliers <- data.frame()
        injected_outliers(saved_outliers)

        class_choices(mode_states[[new_mode]]$classes)

        current_basis(mode_states[[new_mode]]$basis)
        current_projection(mode_states[[new_mode]]$projection)
        current_projection_info(mode_states[[new_mode]]$projection_info)
        current_path(mode_states[[new_mode]]$tour_path)

        # Reset UI element for drawing
        shiny::updateSelectInput(shiny::getDefaultReactiveDomain(), "draw_class", choices = mode_states[[new_mode]]$classes, selected = mode_states[[new_mode]]$classes[1])
        applied_sim_config(mode_states[[new_mode]]$sim_config)
      }

      previous_data_mode(new_mode)

      if (new_mode != "Draw Data") {
        shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "interaction_mode", selected = "Navigate")
      }
    })

    do_clone_to_draw <- function() {
      cd <- current_data()
      cd_test <- current_test_data()
      if (is.null(cd) || nrow(cd) == 0) {
        shiny::showNotification("No data to clone.", type = "warning")
        return()
      }
      if (ncol(cd) < 3) {
        shiny::showNotification("Classbound requires at least 2 feature columns. Cannot clone 1D data to the Draw Canvas.", type = "warning")
        return()
      }
      if (ncol(cd) > 3) {
        shiny::showNotification("High-dimensional data cannot be cloned to the Draw Canvas. The underlying dataset must have exactly two features.", type = "warning")
        return()
      }

      # Clear stale stroke state before mode switch so the Draw canvas starts fresh.
      active_stroke_data(NULL)
      draw_stroke_plot_id(NULL)
      last_sample_pos(NULL)

      # Clone current state to Draw Canvas cache.
      cls <- class_choices()
      mode_states[["Draw Data"]]$data <- cd
      mode_states[["Draw Data"]]$test_data <- cd_test
      mode_states[["Draw Data"]]$classes <- cls

      inj <- injected_outliers()
      if (is.null(inj)) inj <- data.frame()
      mode_states[["Draw Data"]]$outliers <- inj

      # Trigger mode switch.
      shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "data_mode", selected = "Draw Data")
      shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "interaction_mode", selected = "Navigate")
      shiny::showNotification("Data and outliers copied to Draw Canvas!", type = "message")
    }

    shiny::observeEvent(input$clone_to_draw_sim_mvn, {
      do_clone_to_draw()
    })
    shiny::observeEvent(input$clone_to_draw_sim_mix, {
      do_clone_to_draw()
    })
    shiny::observeEvent(input$clone_to_draw_imp, {
      do_clone_to_draw()
    })

    sim_clear_handler <- function() {
      shiny::showModal(shiny::modalDialog(
        title = "Confirm Clear Canvas",
        "Are you sure you want to completely clear the simulated dataset? This action cannot be undone.",
        footer = shiny::tagList(
          shiny::modalButton("Cancel"),
          shiny::actionButton("confirm_sim_reset", "Clear Canvas", class = "btn-danger")
        ),
        size = "s"
      ))
    }

    shiny::observeEvent(input$sim_reset_mvn,
      {
        sim_clear_handler()
      },
      ignoreInit = TRUE
    )
    shiny::observeEvent(input$sim_reset_mix,
      {
        sim_clear_handler()
      },
      ignoreInit = TRUE
    )

    shiny::observeEvent(input$confirm_sim_reset, {
      shiny::removeModal()
      current_data(data.frame(Sim = character(), X1 = numeric(), X2 = numeric()))
      current_test_data(NULL)
      applied_sim_config(NULL)
      zoom_xlim(NULL)
      zoom_ylim(NULL)
      model_cache(list())
    })

    shiny::observeEvent(input$add_class, {
      new_class <- trimws(input$new_class_name)
      if (new_class != "" && !(new_class %in% class_choices())) {
        updated_choices <- c(class_choices(), new_class)
        class_choices(updated_choices)
        shiny::updateSelectInput(shiny::getDefaultReactiveDomain(), "draw_class",
          choices = updated_choices,
          selected = new_class
        )
        shiny::updateTextInput(shiny::getDefaultReactiveDomain(), "new_class_name", value = "")
      }
    })

    # Dynamic Simulation UI (Single View with DOM Preservation)
    output$sim_params_ui <- shiny::renderUI({
      shiny::req(input$sim_n_classes)
      n <- input$sim_n_classes
      if (is.na(n) || n < 1) {
        return(NULL)
      }

      active_i <- shiny::isolate(input$sim_edit_class)
      if (is.null(active_i) || suppressWarnings(as.numeric(active_i)) > n || is.na(as.numeric(active_i))) {
        active_i <- "1"
      }

      shiny::div(
        shiny::selectInput("sim_edit_class", "Active Class to Edit:",
          choices = stats::setNames(as.character(1:n), paste("Class", 1:n)),
          selected = active_i
        ),
        shiny::hr(style = "margin-top: 10px; margin-bottom: 15px;"),
        lapply(1:n, function(i) {
          # Read existing inputs if they exist to prevent wiping user data when `n` changes
          curr_mean <- shiny::isolate(input[[paste0("sim_mean_", i)]])
          curr_sd <- shiny::isolate(input[[paste0("sim_sd_", i)]])
          curr_rho <- shiny::isolate(input[[paste0("sim_cor_", i)]])
          curr_n <- shiny::isolate(input[[paste0("sim_n_", i)]])

          if (is.null(curr_mean)) curr_mean <- if (i == 1) "-1,0" else if (i == 2) "1,0" else "0,1"
          if (is.null(curr_sd)) curr_sd <- "1,1"
          if (is.null(curr_rho)) curr_rho <- 0
          if (is.null(curr_n)) curr_n <- 100

          shiny::conditionalPanel(
            condition = paste0("input.sim_edit_class == '", i, "'"),
            shiny::wellPanel(
              style = "margin-bottom: 0px;",
              shiny::textInput(paste0("sim_mean_", i), "Mean (comma separated)", value = curr_mean),
              shiny::textInput(paste0("sim_sd_", i), "Standard Deviations (comma separated)", value = curr_sd),
              shiny::numericInput(paste0("sim_cor_", i), "Correlation (rho)", value = curr_rho, min = -0.99, max = 0.99, step = 0.1),
              shiny::numericInput(paste0("sim_n_", i), "Sample Size", value = curr_n, min = 10)
            )
          )
        })
      )
    })

    shiny::observeEvent(input$sim_do, {
      shiny::req(input$sim_do > 0)
      shiny::req(input$sim_n_classes)
      n <- input$sim_n_classes
      if (is.na(n) || n < 2) {
        shiny::showNotification("At least 2 classes are required to generate classification data.", type = "error")
        return()
      }

      means <- list()
      covs <- list()
      ns <- c()

      for (i in 1:n) {
        m_str <- input[[paste0("sim_mean_", i)]]
        s_str <- input[[paste0("sim_sd_", i)]]

        if (is.null(m_str) || is.null(s_str)) {
          return()
        }

        m <- suppressWarnings(as.numeric(unlist(strsplit(m_str, ","))))
        sds <- suppressWarnings(as.numeric(unlist(strsplit(s_str, ","))))
        rho <- input[[paste0("sim_cor_", i)]]
        ns <- c(ns, input[[paste0("sim_n_", i)]])

        means[[i]] <- m

        # Construct covariance matrix using SDs and exchangeable correlation
        dim_m <- length(m)
        if (length(sds) != dim_m) {
          shiny::showNotification(sprintf("Error in Class %d: Mean has %d dimension(s), but SD has %d dimension(s).", i, dim_m, length(sds)), type = "error")
          return()
        }

        cov_mat <- matrix(rho, nrow = dim_m, ncol = dim_m)
        diag(cov_mat) <- 1

        # Scale correlation matrix to covariance matrix: Cov_ij = rho * sd_i * sd_j
        sd_diag <- diag(sds, dim_m)
        cov_mat <- sd_diag %*% cov_mat %*% sd_diag

        covs[[i]] <- cov_mat
      }

      lengths <- sapply(means, length)
      if (length(unique(lengths)) > 1) {
        # Find the first class that doesn't match Class 1 to give a highly specific error
        target_dim <- lengths[1]
        for (i in 2:length(lengths)) {
          if (lengths[i] != target_dim) {
            shiny::showNotification(sprintf("Error: Class 1 has %d dimension(s), but Class %d has %d dimension(s). All classes must match.", target_dim, i, lengths[i]), type = "error")
            return()
          }
        }
      }

      for (i in 1:length(means)) {
        if (any(is.na(means[[i]])) || any(is.na(covs[[i]]))) {
          shiny::showNotification(sprintf("Error in Class %d: Contains invalid text or missing numbers in Mean or SD.", i), type = "error")
          return()
        }
      }

      new_data <- tryCatch(
        {
          seed_val <- if (is.numeric(input$sim_seed) && !is.na(input$sim_seed)) as.integer(input$sim_seed) else NULL
          noise_val <- if (is.numeric(input$sim_noise)) input$sim_noise / 100 else 0
          simu_n(means = means, covs = covs, ns = ns, seed = seed_val, noise_ratio = noise_val, test_ratio = 0.3)
        },
        error = function(e) {
          clean_msg <- clean_err_msg(e$message)
          shiny::showNotification(paste("Data generation failed:", clean_msg), type = "error")
          NULL
        }
      )

      if (is.null(new_data)) {
        return()
      }

      if (is.list(new_data) && !is.data.frame(new_data)) {
        current_data(new_data$train)
        current_test_data(new_data$test)
      } else {
        current_data(new_data)
        current_test_data(NULL)
      }

      conf <- list(
        engine = "mvn",
        seed = seed_val,
        noise_ratio = noise_val,
        n_classes = n
      )
      for (i in 1:n) {
        conf[[paste0("class_", i)]] <- list(
          mean = input[[paste0("sim_mean_", i)]],
          sd = input[[paste0("sim_sd_", i)]],
          cor = input[[paste0("sim_cor_", i)]],
          n = input[[paste0("sim_n_", i)]]
        )
      }
      applied_sim_config(conf)

      injected_outliers(data.frame())
      train_dat <- if (is.list(new_data) && !is.data.frame(new_data)) new_data$train else new_data
      new_classes <- unique(as.character(train_dat$Sim))
      class_choices(new_classes)
      shiny::updateSelectInput(shiny::getDefaultReactiveDomain(), "draw_class", choices = new_classes, selected = new_classes[1])
      zoom_xlim(NULL)
      zoom_ylim(NULL)
    })

    shiny::observeEvent(input$sim_mixsim_do, {
      shiny::req(input$sim_mixsim_do > 0)
      shiny::req(input$sim_mixsim_k, input$sim_mixsim_p, input$sim_mixsim_omega, input$sim_mixsim_n)

      if (is.na(input$sim_mixsim_k) || input$sim_mixsim_k < 2) {
        shiny::showNotification("At least 2 classes are required to generate classification data.", type = "error")
        return()
      }

      new_data <- tryCatch(
        {
          seed_val <- if (is.numeric(input$sim_seed) && !is.na(input$sim_seed)) as.integer(input$sim_seed) else NULL
          noise_val <- if (is.numeric(input$sim_noise)) input$sim_noise / 100 else 0
          simulate_mixsim(
            n = input$sim_mixsim_n,
            K = input$sim_mixsim_k,
            p = input$sim_mixsim_p,
            MaxOmega = input$sim_mixsim_omega,
            seed = seed_val,
            noise_ratio = noise_val,
            test_ratio = 0.3
          )
        },
        error = function(e) {
          clean_msg <- clean_err_msg(e$message)
          shiny::showNotification(paste("MixSim data generation failed:", clean_msg), type = "error")
          NULL
        }
      )

      if (is.null(new_data)) {
        return()
      }

      if (is.list(new_data) && !is.data.frame(new_data)) {
        current_data(new_data$train)
        current_test_data(new_data$test)
      } else {
        current_data(new_data)
        current_test_data(NULL)
      }

      applied_sim_config(list(
        engine = "mixsim",
        seed = seed_val,
        noise_ratio = noise_val,
        n = input$sim_mixsim_n,
        k = input$sim_mixsim_k,
        p = input$sim_mixsim_p,
        omega = input$sim_mixsim_omega
      ))

      injected_outliers(data.frame())
      train_dat <- if (is.list(new_data) && !is.data.frame(new_data)) new_data$train else new_data
      new_classes <- unique(as.character(train_dat$Sim))
      class_choices(new_classes)
      shiny::updateSelectInput(shiny::getDefaultReactiveDomain(), "draw_class", choices = new_classes, selected = new_classes[1])
      zoom_xlim(NULL)
      zoom_ylim(NULL)
    })

    shiny::observeEvent(input$undo_draw, {
      hist <- undo_history()
      if (length(hist) > 0) {
        last_state <- hist[[length(hist)]]
        current_data(last_state)
        hist[[length(hist)]] <- NULL
        undo_history(hist)
      }
    })

    shiny::observeEvent(input$clear,
      {
        shiny::showModal(shiny::modalDialog(
          title = "Confirm Clear Canvas",
          "Are you sure you want to completely clear the drawn dataset? This action cannot be undone.",
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton("confirm_draw_reset", "Clear Canvas", class = "btn-danger")
          ),
          size = "s"
        ))
      },
      ignoreInit = TRUE
    )

    shiny::observeEvent(input$confirm_draw_reset,
      {
        shiny::removeModal()
        current_data(data.frame(Sim = character(), X1 = numeric(), X2 = numeric()))
        current_test_data(NULL)
        injected_outliers(data.frame())

        new_classes <- c("Class 1", "Class 2", "Class 3")
        mode_states[["Draw Data"]]$classes <- new_classes
        class_choices(new_classes)

        shiny::updateSelectInput(shiny::getDefaultReactiveDomain(), "draw_class", choices = new_classes, selected = "Class 1")
        zoom_xlim(NULL)
        zoom_ylim(NULL)
      },
      ignoreInit = TRUE
    )

    last_valid_draw_classes <- shiny::reactiveVal(3)
    pending_draw_classes <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$draw_total_classes,
      {
        shiny::req(input$draw_total_classes > 0)
        n <- as.integer(input$draw_total_classes)
        old_n <- last_valid_draw_classes()

        if (n == old_n) {
          return()
        }

        new_classes <- paste0("Class ", 1:n)
        old_classes <- class_choices()
        custom_classes <- old_classes[!grepl("^Class \\d+$", old_classes)]
        final_classes <- unique(c(new_classes, custom_classes))

        dat <- current_data()
        inj <- injected_outliers()

        removed_classes <- setdiff(old_classes, final_classes)
        has_removed_data <- FALSE
        if (length(removed_classes) > 0) {
          if ("Sim" %in% colnames(dat) && any(dat$Sim %in% removed_classes)) has_removed_data <- TRUE
          if ("Sim" %in% colnames(inj) && any(inj$Sim %in% removed_classes)) has_removed_data <- TRUE
        }

        if (has_removed_data) {
          pending_draw_classes(n)
          shiny::updateNumericInput(shiny::getDefaultReactiveDomain(), "draw_total_classes", value = old_n)

          shiny::showModal(shiny::modalDialog(
            title = "Warning: Data Deletion",
            "Reducing the number of base classes will permanently delete the points drawn for the removed classes. Are you sure you want to proceed?",
            footer = shiny::tagList(
              shiny::actionButton("cancel_remove_classes", "Cancel"),
              shiny::actionButton("confirm_remove_classes", "Yes, Delete Data", class = "btn-danger", style = "color: white; background-color: #d9534f; border-color: #d43f3a;")
            )
          ))
        } else {
          last_valid_draw_classes(n)
          mode_states[["Draw Data"]]$classes <- final_classes
          class_choices(final_classes)
          shiny::updateSelectInput(shiny::getDefaultReactiveDomain(), "draw_class",
            choices = final_classes, selected = final_classes[1]
          )
        }
      },
      ignoreInit = TRUE
    )

    shiny::observeEvent(input$cancel_remove_classes, {
      shiny::removeModal()
      old_n <- last_valid_draw_classes()
      shiny::updateNumericInput(shiny::getDefaultReactiveDomain(), "draw_total_classes", value = old_n)
    })

    shiny::observeEvent(input$confirm_remove_classes, {
      shiny::removeModal()
      n <- pending_draw_classes()
      shiny::req(n)

      last_valid_draw_classes(n)
      shiny::updateNumericInput(shiny::getDefaultReactiveDomain(), "draw_total_classes", value = n)

      new_classes <- paste0("Class ", 1:n)
      old_classes <- class_choices()
      custom_classes <- old_classes[!grepl("^Class \\d+$", old_classes)]
      final_classes <- unique(c(new_classes, custom_classes))

      # Hard delete data
      dat <- current_data()
      inj <- injected_outliers()
      if ("Sim" %in% colnames(dat)) dat <- dat[dat$Sim %in% final_classes, , drop = FALSE]
      if ("Sim" %in% colnames(inj)) inj <- inj[inj$Sim %in% final_classes, , drop = FALSE]

      current_data(dat)
      injected_outliers(inj)

      mode_states[["Draw Data"]]$classes <- final_classes
      class_choices(final_classes)
      shiny::updateSelectInput(shiny::getDefaultReactiveDomain(), "draw_class",
        choices = final_classes, selected = final_classes[1]
      )
    })

    shiny::observeEvent(input$plot_dblclick, {
      zoom_xlim(NULL)
      zoom_ylim(NULL)
    })

    shiny::observeEvent(input$plot_wheel, {
      w <- input$plot_wheel
      shiny::req(w)
      h <- input[[paste0("hover_", w$plot_id)]]
      shiny::req(h)

      curr_x <- zoom_xlim()
      if (is.null(curr_x)) curr_x <- c(h$domain$left, h$domain$right)
      curr_y <- zoom_ylim()
      if (is.null(curr_y)) curr_y <- c(h$domain$bottom, h$domain$top)

      factor <- if (w$direction > 0) 1.25^w$direction else 0.8^abs(w$direction)

      new_width <- abs(curr_x[2] - curr_x[1]) * factor
      new_height <- abs(curr_y[2] - curr_y[1]) * factor

      prop_x <- (h$x - min(curr_x)) / abs(curr_x[2] - curr_x[1])
      prop_y <- (h$y - min(curr_y)) / abs(curr_y[2] - curr_y[1])

      new_xmin <- h$x - (new_width * prop_x)
      new_xmax <- h$x + (new_width * (1 - prop_x))

      new_ymin <- h$y - (new_height * prop_y)
      new_ymax <- h$y + (new_height * (1 - prop_y))

      if (curr_x[1] > curr_x[2]) {
        zoom_xlim(c(max(new_xmin, new_xmax), min(new_xmin, new_xmax)))
      } else {
        zoom_xlim(c(min(new_xmin, new_xmax), max(new_xmin, new_xmax)))
      }

      if (curr_y[1] > curr_y[2]) {
        zoom_ylim(c(max(new_ymin, new_ymax), min(new_ymin, new_ymax)))
      } else {
        zoom_ylim(c(min(new_ymin, new_ymax), max(new_ymin, new_ymax)))
      }
    })

    shiny::observeEvent(input$plot_brush, {
      if (input$data_mode != "Draw Data" || input$interaction_mode == "Navigate") {
        b <- input$plot_brush
        if (!is.null(b)) {
          # Because the physical plot panel is forced to be square via aspect.ratio=1,
          # we must expand the brushed region to match the data's aspect ratio to avoid distortion.
          orig_w <- abs(b$domain$right - b$domain$left)
          orig_h <- abs(b$domain$top - b$domain$bottom)

          if (orig_w > 0 && orig_h > 0) {
            b_w <- abs(b$xmax - b$xmin)
            b_h <- abs(b$ymax - b$ymin)

            ratio <- orig_w / orig_h
            if (b_w / b_h > ratio) {
              # Brush is too wide relative to domain, expand height proportionally
              new_w <- b_w
              new_h <- b_w / ratio
            } else {
              # Brush is too tall relative to domain, expand width proportionally
              new_h <- b_h
              new_w <- b_h * ratio
            }

            mid_x <- (b$xmin + b$xmax) / 2
            mid_y <- (b$ymin + b$ymax) / 2

            zoom_xlim(c(mid_x - new_w / 2, mid_x + new_w / 2))
            zoom_ylim(c(mid_y - new_h / 2, mid_y + new_h / 2))
          } else {
            zoom_xlim(c(b$xmin, b$xmax))
            zoom_ylim(c(b$ymin, b$ymax))
          }
        }
      }
    })

    # Single Point Draw Observer (uses native Shiny click)
    shiny::observeEvent(input$plot_click, {
      if (input$data_mode != "Draw Data" || input$interaction_mode != "Draw Point") {
        return()
      }
      if (ncol(current_data()) > 3) {
        return()
      }

      click <- input$plot_click
      if (is.null(click) || is.null(click$x) || is.null(click$y)) {
        return()
      }

      cd <- current_data()
      feat_cols <- setdiff(colnames(cd), "Sim")

      pts <- data.frame(Sim = shiny::isolate(input$draw_class))
      pts[[feat_cols[1]]] <- click$x
      pts[[feat_cols[2]]] <- click$y

      # Lock limits to prevent jumping only if we already have a stable coordinate space
      if (is.null(shiny::isolate(zoom_xlim())) && nrow(cd) > 1) {
        if (!is.null(click$domain)) {
          zoom_xlim(c(click$domain$left, click$domain$right))
          zoom_ylim(c(click$domain$bottom, click$domain$top))
        }
      }

      hist <- undo_history()
      hist[[length(hist) + 1]] <- cd
      undo_history(hist)
      current_data(rbind(cd, pts))
    })

    # Freehand Draw Observers (Cluster only)
    shiny::observeEvent(input$draw_stroke_start, {
      if (input$data_mode != "Draw Data" || input$interaction_mode != "Draw Cluster") {
        return()
      }
      if (ncol(current_data()) > 3) {
        return()
      }

      plot_id <- input$draw_stroke_start$plot_id
      draw_stroke_plot_id(plot_id)

      if (!is.null(input$draw_stroke_start$mapped)) {
        active_stroke_data(data.frame(x = input$draw_stroke_start$mapped$x, y = input$draw_stroke_start$mapped$y))
      } else {
        active_stroke_data(data.frame())
      }
      last_sample_pos(NULL)
    })

    shiny::observeEvent(input$draw_stroke_end, {
      if (input$data_mode != "Draw Data" || input$interaction_mode != "Draw Cluster") {
        return()
      }

      stk <- active_stroke_data()
      if (is.null(stk)) {
        return()
      }

      # Handle Click without movement: if stroke is empty, generate exactly one cluster/point at last known hover
      if (nrow(stk) == 0) {
        pid <- draw_stroke_plot_id()
        if (!is.null(pid)) {
          h <- input[[paste0("hover_", pid)]]
          if (!is.null(h)) {
            cd <- current_data()
            feat_cols <- setdiff(colnames(cd), "Sim")

            is_cluster <- shiny::isolate(input$interaction_mode) == "Draw Cluster"
            n <- if (is_cluster) shiny::isolate(input$brush_size) else 1
            spread <- if (is_cluster) shiny::isolate(input$brush_spread) else 0.05

            # Screen-space Brush Conversion (used in both JS and R)
            # Formula: R_pixels = max(5, spread * 400)
            r_pixels <- max(5, spread * 400)

            pixel_width <- if (!is.null(h$range)) abs(h$range$right - h$range$left) else 500
            pixel_height <- if (!is.null(h$range)) abs(h$range$bottom - h$range$top) else 500
            data_width <- if (!is.null(h$domain)) abs(h$domain$right - h$domain$left) else 1
            data_height <- if (!is.null(h$domain)) abs(h$domain$top - h$domain$bottom) else 1

            data_x_per_pixel <- data_width / pixel_width
            data_y_per_pixel <- data_height / pixel_height

            sd_x <- r_pixels * data_x_per_pixel
            sd_y <- r_pixels * data_y_per_pixel
            if (!is_cluster) {
              sd_x <- 0
              sd_y <- 0
            }

            pts <- data.frame(Sim = rep(shiny::isolate(input$draw_class), n))
            pts[[feat_cols[1]]] <- stats::rnorm(n, mean = h$x, sd = sd_x)
            pts[[feat_cols[2]]] <- stats::rnorm(n, mean = h$y, sd = sd_y)
            stk <- pts
          }
        }
      }

      # If the plot is auto-scaling (limits are NULL), lock the limits to exactly what the
      # Lock axis limits so the newly drawn stroke isn't distorted by auto-scaling.
      if (is.null(shiny::isolate(zoom_xlim()))) {
        pid <- draw_stroke_plot_id()
        h_end <- shiny::isolate(input[[paste0("hover_", pid)]])
        if (!is.null(h_end$domain)) {
          zoom_xlim(c(h_end$domain$left, h_end$domain$right))
          zoom_ylim(c(h_end$domain$bottom, h_end$domain$top))
        }
      }

      if (nrow(stk) > 0) {
        cd <- current_data()
        hist <- undo_history()
        hist[[length(hist) + 1]] <- cd
        undo_history(hist)

        current_data(rbind(cd, stk))
      }

      active_stroke_data(NULL)
      draw_stroke_plot_id(NULL)
      last_sample_pos(NULL)
    })

    shiny::observe({
      pid <- draw_stroke_plot_id()
      shiny::req(pid)
      if (input$data_mode != "Draw Data" || input$interaction_mode != "Draw Cluster") {
        return()
      }

      h <- input[[paste0("hover_", pid)]]
      shiny::req(h)

      last_pos <- last_sample_pos()

      is_cluster <- shiny::isolate(input$interaction_mode) == "Draw Cluster"
      spread <- if (is_cluster) shiny::isolate(input$brush_spread) else 0.05

      # Screen-space Brush Conversion (used in both JS and R)
      # Formula: R_pixels = max(5, spread * 400)
      r_pixels <- max(5, spread * 400)

      pixel_width <- if (!is.null(h$range)) abs(h$range$right - h$range$left) else 500
      pixel_height <- if (!is.null(h$range)) abs(h$range$bottom - h$range$top) else 500
      data_width <- if (!is.null(h$domain)) abs(h$domain$right - h$domain$left) else 1
      data_height <- if (!is.null(h$domain)) abs(h$domain$top - h$domain$bottom) else 1

      data_x_per_pixel <- data_width / pixel_width
      data_y_per_pixel <- data_height / pixel_height

      if (!is.null(last_pos)) {
        # Screen-space interpolation distance
        pixel_dx <- (h$x - last_pos$x) / data_x_per_pixel
        pixel_dy <- (h$y - last_pos$y) / data_y_per_pixel
        pixel_dist <- sqrt(pixel_dx^2 + pixel_dy^2)

        # Fixed threshold of ~1% of plot's relevant pixel dimension
        threshold_pixels <- max(pixel_width, pixel_height) * 0.01
        if (pixel_dist < threshold_pixels) {
          return()
        }

        num_steps <- floor(pixel_dist / threshold_pixels)
        if (num_steps > 100) num_steps <- 100

        centers <- data.frame(
          x = seq(last_pos$x, h$x, length.out = num_steps + 1)[-1],
          y = seq(last_pos$y, h$y, length.out = num_steps + 1)[-1]
        )
      } else {
        centers <- data.frame(x = h$x, y = h$y)
      }

      last_sample_pos(list(x = h$x, y = h$y))

      n <- if (is_cluster) shiny::isolate(input$brush_size) else 1

      sd_x <- r_pixels * data_x_per_pixel
      sd_y <- r_pixels * data_y_per_pixel
      if (!is_cluster) {
        sd_x <- 0
        sd_y <- 0
      }

      cd <- shiny::isolate(current_data())
      feat_cols <- setdiff(colnames(cd), "Sim")

      pts_list <- lapply(seq_len(nrow(centers)), function(i) {
        p <- data.frame(Sim = rep(shiny::isolate(input$draw_class), n))
        p[[feat_cols[1]]] <- stats::rnorm(n, mean = centers$x[i], sd = sd_x)
        p[[feat_cols[2]]] <- stats::rnorm(n, mean = centers$y[i], sd = sd_y)
        p
      })

      pts <- do.call(rbind, pts_list)

      curr_stk <- shiny::isolate(active_stroke_data())
      if (is.null(curr_stk)) curr_stk <- data.frame()
      active_stroke_data(rbind(curr_stk, pts))
    })

    # State reset on canvas clear, mode switch, dataset switch
    shiny::observeEvent(c(input$clear, input$data_mode, input$interaction_mode), {
      active_stroke_data(NULL)
      draw_stroke_plot_id(NULL)
      last_sample_pos(NULL)
    })

    # Helper to standardize formatting across all data preview tables
    create_preview_table <- function(dat, empty_msg) {
      if (is.null(dat) || nrow(dat) == 0) {
        return(DT::datatable(data.frame(Message = empty_msg), options = list(dom = "t"), rownames = FALSE))
      }
      num_cols <- sapply(dat, is.numeric)
      dat[num_cols] <- lapply(dat[num_cols], round, 2)
      DT::datatable(
        dat,
        filter = "top",
        options = list(pageLength = 5, lengthMenu = c(5, 10, 20), scrollX = TRUE),
        rownames = FALSE
      )
    }

    output$drawn_points_table <- DT::renderDataTable(
      {
        create_preview_table(combined_training_data(), "No points drawn yet")
      },
      server = FALSE
    )

    output$import_data_table <- DT::renderDataTable(
      {
        create_preview_table(combined_training_data(), "No imported data available")
      },
      server = FALSE
    )

    output$sim_data_table_train <- DT::renderDataTable(
      {
        create_preview_table(combined_training_data(), "No simulated training data available")
      },
      server = FALSE
    )

    output$sim_data_table_test <- DT::renderDataTable(
      {
        create_preview_table(current_test_data(), "No simulated testing data available")
      },
      server = FALSE
    )

    output$workspace_import_ui <- shiny::renderUI({
      input$refresh_ws # Take dependency on refresh button
      models <- find_workspace_models(.GlobalEnv)

      shiny::tagList(
        if (length(models) > 0) {
          shiny::selectInput("ws_model", "Select Workspace Object", choices = models)
        } else {
          shiny::p("No workflow or model_fit objects found in global environment.")
        },
        shiny::div(
          style = "display: flex; gap: 10px; flex-wrap: wrap;",
          shiny::actionButton("refresh_ws", "Refresh", icon = shiny::icon("sync")),
          if (length(models) > 0) shiny::actionButton("add_ws_model", "Add to Comparison", class = "btn-primary")
        )
      )
    })

    shiny::observeEvent(input$add_ws_model, {
      shiny::req(input$ws_model)
      obj_name <- input$ws_model

      tryCatch(
        {
          obj <- get(obj_name, envir = .GlobalEnv)
          friendly_name <- paste0("Workspace: ", obj_name)

          is_new <- !(friendly_name %in% names(app_methods))

          # Add to the closure variables
          app_methods[[friendly_name]] <<- list(
            fn = obj,
            args = list(),
            supports_prob = inherits(obj, "workflow") # Generally workflows support probability
          )
          predict_args[[friendly_name]] <<- function(...) list()

          # Update UI choices
          shiny::updateCheckboxGroupInput(
            shiny::getDefaultReactiveDomain(),
            "selected_models",
            choices = names(app_methods),
            selected = unique(c(input$selected_models, friendly_name))
          )

          if (is_new) {
            shiny::showNotification(paste("Added", obj_name, "to models."), type = "message")
          } else {
            shiny::showNotification(paste("Updated existing model", obj_name, "from workspace."), type = "message")
          }

          # Invalidate comparison_state to force re-fit.
          ws_update_trigger(ws_update_trigger() + 1)
        },
        error = function(e) {
          clean_msg <- clean_err_msg(e$message)
          shiny::showNotification(paste("Failed to import model:", clean_msg), type = "error")
        }
      )
    })

    # Debounced model selection to coalesce rapid checkbox changes
    selected_models_d <- shiny::debounce(shiny::reactive({
      input$selected_models
    }), 400)

    output$plot_grid <- shiny::renderUI({
      models <- selected_models_d()
      if (is.null(models) || length(models) == 0) {
        return(shiny::p("Please select at least one model to compare."))
      }

      plot_outputs <- lapply(models, function(m) {
        safe_id <- gsub("[^a-zA-Z0-9_\\-]", "_", paste0("plot_", m))
        shiny::div(
          # Aspect ratio wrapper to keep container roughly square to match ggplot2 aspect.ratio=1
          style = "width: 100%; max-width: 450px; aspect-ratio: 1.15 / 1; margin: 0 auto;",
          shiny::plotOutput(
            outputId = safe_id,
            height = "100%",
            click = "plot_click",
            dblclick = "plot_dblclick",
            brush = shiny::brushOpts(id = "plot_brush", resetOnNew = TRUE),
            hover = shiny::hoverOpts(id = paste0("hover_", safe_id), delay = 50, delayType = "throttle")
          )
        )
      })

      shiny::div(
        # Responsive grid: centers the plots, prevents extreme horizontal stretching
        style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(350px, 1fr)); gap: 16px; align-items: center; justify-items: center;",
        plot_outputs
      )
    })

    # Manual Tour Steering (High-Dimensional Mode)

    current_basis <- shiny::reactiveVal(NULL)
    current_path <- shiny::reactiveVal(NULL)
    current_projection <- shiny::reactiveVal(NULL)
    current_projection_info <- shiny::reactiveVal(NULL)

    # Reactive for previous slice selections (collision avoidance)
    prev_slice <- shiny::reactiveValues(x = NULL, y = NULL)

    output$tour_panel <- shiny::renderUI({
      dat <- current_data()
      if (ncol(dat) > 3) {
        mode <- shiny::isolate(input$data_mode)
        saved_var <- mode_states[[mode]]$tour_var
        saved_angle <- mode_states[[mode]]$tour_angle

        feat_choices <- setdiff(colnames(dat), "Sim")
        num_choices <- names(which(sapply(dat[feat_choices], is.numeric)))
        if (is.null(saved_var) || !(saved_var %in% feat_choices)) {
          saved_var <- feat_choices[1]
        }
        if (is.null(saved_angle)) saved_angle <- 0

        tourr_available <- requireNamespace("tourr", quietly = TRUE)
        hd_modes <- if (tourr_available) c("Projection" = "projection", "2D Slice" = "slice") else c("2D Slice" = "slice")
        default_hd <- if (tourr_available) "projection" else "slice"

        # Default slice X/Y to first two numeric predictors
        default_sx <- if (length(num_choices) >= 1) num_choices[1] else feat_choices[1]
        default_sy <- if (length(num_choices) >= 2) num_choices[2] else if (length(feat_choices) >= 2) feat_choices[2] else feat_choices[1]

        shiny::wellPanel(
          shiny::tags$details(
            open = "open",
            shiny::tags$summary("High-Dimensional Visualization", style = "display: list-item; font-size: 18px; font-weight: 500; cursor: pointer; margin-bottom: 10px;"),
            shiny::p("High-dimensional data detected."),
            shiny::radioButtons("hd_view_mode", "Visualization Mode", choices = hd_modes, selected = default_hd, inline = TRUE),
            if (!tourr_available) shiny::p(shiny::tags$em("Note: Install the 'tourr' package to enable Projection mode."), style = "color: #888; font-size: 12px;"),
            shiny::conditionalPanel(
              condition = "input.hd_view_mode == 'projection'",
              shiny::selectInput("tour_var", "Manipulation Variable", choices = feat_choices, selected = saved_var),
              shiny::sliderInput("tour_angle", "Rotation Angle", min = 0, max = 1, value = saved_angle, step = 0.01),
              shiny::uiOutput("proj_note_ui")
            ),
            shiny::conditionalPanel(
              condition = "input.hd_view_mode == 'slice'",
              shiny::selectInput("slice_x", "X-axis Feature", choices = num_choices, selected = default_sx),
              shiny::selectInput("slice_y", "Y-axis Feature", choices = num_choices, selected = default_sy),
              shiny::uiOutput("slice_note_ui")
            )
          )
        )
      }
    })

    # Slice feature collision observer
    shiny::observe({
      sx <- input$slice_x
      sy <- input$slice_y
      if (is.null(sx) || is.null(sy)) {
        return()
      }

      dat <- current_data()
      if (is.null(dat) || ncol(dat) <= 3) {
        return()
      }
      feat_choices <- setdiff(colnames(dat), "Sim")
      num_choices <- names(which(sapply(dat[feat_choices], is.numeric)))

      if (sx == sy && length(num_choices) > 1) {
        # Determine which changed
        if (!is.null(prev_slice$x) && sx != prev_slice$x) {
          # X was changed to match Y, reassign Y
          alt <- setdiff(num_choices, sx)[1]
          shiny::updateSelectInput(shiny::getDefaultReactiveDomain(), "slice_y", selected = alt)
        } else {
          # Y was changed to match X, reassign X
          alt <- setdiff(num_choices, sy)[1]
          shiny::updateSelectInput(shiny::getDefaultReactiveDomain(), "slice_x", selected = alt)
        }
      }
      prev_slice$x <- sx
      prev_slice$y <- sy
    })

    # Dynamic explanatory note for 2D Slice
    output$slice_note_ui <- shiny::renderUI({
      sx <- input$slice_x
      sy <- input$slice_y
      if (is.null(sx) || is.null(sy)) {
        return(NULL)
      }
      shiny::p(
        shiny::tags$em(paste0(
          "2D Slice: Shows ", sx, " vs ", sy,
          ". Other numeric features are fixed at their median and categorical features at their mode. ",
          "This is a 2D slice of the full multivariate boundary, so it may not represent the complete decision boundary accurately."
        )),
        style = "color: #666; font-size: 12px; margin-top: 5px;"
      )
    })

    # Dynamic explanatory note for Projection
    output$proj_note_ui <- shiny::renderUI({
      shiny::p(
        shiny::tags$em("Projection: Shows a 2D linear projection of the high-dimensional space. All features contribute simultaneously via a weighted linear combination."),
        style = "color: #666; font-size: 12px; margin-top: 5px;"
      )
    })

    # Initialize PCA basis when high-dim data loads
    shiny::observe({
      dat <- current_data()
      if (ncol(dat) > 3 && nrow(dat) >= 2 && is.null(current_basis())) {
        if (!requireNamespace("tourr", quietly = TRUE)) {
          # tourr unavailable; 2D Slice mode will be used instead
          shiny::showNotification("The 'tourr' package is not installed. Falling back to 2D Slice mode for high-dimensional data.", type = "error", duration = 10)
          return()
        }

        num_dat <- dat[, setdiff(colnames(dat), "Sim"), drop = FALSE]

        # Check variance to prevent scaling errors
        vars <- apply(num_dat, 2, stats::var)
        scale_flag <- if (any(vars == 0, na.rm = TRUE)) FALSE else TRUE

        if (!scale_flag) {
          shiny::showNotification("Notice: One or more features have zero variance. PCA scaling has been disabled.", type = "warning")
        }

        pca <- stats::prcomp(num_dat, scale. = scale_flag)
        basis <- pca$rotation[, 1:2]
        basis <- tourr::orthonormalise(basis)
        current_basis(basis)
        current_projection(basis)
        scale_val <- if (is.numeric(pca$scale)) pca$scale else NULL
        current_projection_info(list(center = pca$center, scale = scale_val))
      } else if (ncol(dat) <= 3) {
        current_basis(NULL)
        current_path(NULL)
        current_projection(NULL)
        current_projection_info(NULL)
      }
    })

    # Generate Geodesic Path when variable changes
    shiny::observeEvent(input$tour_var, {
      shiny::req(input$tour_var)
      shiny::req(current_basis())
      if (!requireNamespace("tourr", quietly = TRUE)) {
        shiny::showNotification("The 'tourr' package is required for projection steering. Please install it.", type = "error")
        return()
      }

      dat <- current_data()
      num_cols <- setdiff(colnames(dat), "Sim")
      var_idx <- which(num_cols == input$tour_var)

      # Orthonormalize target basis.
      start_basis <- current_basis()
      target_basis <- start_basis
      target_basis[var_idx, ] <- 0
      target_basis <- tourr::orthonormalise(target_basis)

      path <- tourr::geodesic_path(start_basis, target_basis)
      current_path(path)

      # Reset slider when variable changes
      shiny::updateSliderInput(shiny::getDefaultReactiveDomain(), "tour_angle", value = 0)
    })

    # Interpolate basis when slider moves
    shiny::observeEvent(input$tour_angle, {
      shiny::req(current_path())

      frac <- input$tour_angle
      path <- current_path()

      new_proj <- path$interpolate(frac)
      current_projection(new_proj)
    })

    # Model cache for hash-based reuse across checkbox changes
    model_cache <- shiny::reactiveVal(list())

    comparison_state <- shiny::reactive({
      ws_update_trigger()
      models <- selected_models_d()
      dat <- combined_training_data()
      if (is.null(models) || length(models) == 0 || nrow(dat) < 2 || length(unique(dat$Sim)) < 2) {
        return(NULL)
      }

      state <- list(train_data = dat, models = list())
      cache <- shiny::isolate(model_cache())

      for (m in models) {
        config <- app_methods[[m]]
        if (!is.null(config)) {
          plot_opts <- list()
          if (!is.null(config$fit_args_fn) && is.function(config$fit_args_fn)) {
            custom_opts <- config$fit_args_fn(input)
            plot_opts <- c(plot_opts, custom_opts)
          }

          fit_args <- utils::modifyList(config$args, plot_opts)
          pred_args <- predict_args[[m]](as.numeric(input$rule))

          # Deterministic cache key
          cache_key <- rlang::hash(list(data = dat, model = m, fit_args = fit_args, pred_args = pred_args))

          if (!is.null(cache[[cache_key]])) {
            state$models[[m]] <- cache[[cache_key]]
          } else {
            tryCatch(
              {
                cb_mod <- fit_model(dat, Sim ~ ., classifier = config$fn, fit_args = fit_args)
                preds <- predict_model(cb_mod, dat, predict_args = pred_args)

                test_dat <- shiny::isolate(current_test_data())
                test_preds <- NULL
                if (!is.null(test_dat) && nrow(test_dat) > 0) {
                  test_preds <- predict_model(cb_mod, test_dat, predict_args = pred_args)
                }

                entry <- list(
                  model = cb_mod,
                  predictions = preds,
                  predict_args = pred_args,
                  test_predictions = test_preds
                )
                state$models[[m]] <- entry
                cache[[cache_key]] <- entry
              },
              error = function(e) {
                clean_msg <- clean_err_msg(e$message)
                shiny::showNotification(paste("Model", m, "failed to fit:", clean_msg), type = "error", duration = 10, id = paste0("err_", m))
              }
            )
          }
        }
      }

      model_cache(cache)
      state
    })

    # Register renderPlot for selected models (uses debounced selection)
    shiny::observe({
      models <- selected_models_d()
      if (is.null(models)) {
        return()
      }

      lapply(models, function(m) {
        safe_id <- gsub("[^a-zA-Z0-9_\\-]", "_", paste0("plot_", m))
        output[[safe_id]] <- shiny::renderPlot({
          # Graceful abort if model was deselected
          if (!(m %in% selected_models_d())) {
            return(NULL)
          }

          dat <- combined_training_data()
          if (!is.null(dat) && nrow(dat) > 0) {
            shiny::validate(shiny::need(ncol(dat) >= 3, "Classbound requires at least 2 feature columns for 2D visualization."))
          }

          title <- switch(m,
            "rpart" = "Decision Tree (rpart)",
            "PPtreeViz" = "PPtreeViz",
            "PPtreeExtclass" = "PPtreeExtclass",
            "PPtreeExt_split" = "PPtreeExt_split",
            "randomForest" = "Random Forest",
            m
          )

          if (!is.null(dat) && nrow(dat) > 0) {
            dat$Sim <- factor(dat$Sim, levels = class_choices())
          }

          if (is.null(dat) || nrow(dat) < 2 || length(unique(dat$Sim)) < 2) {
            feat_cols <- setdiff(colnames(dat), "Sim")
            x_name <- if (length(feat_cols) >= 1) feat_cols[1] else "X1"
            y_name <- if (length(feat_cols) >= 2) feat_cols[2] else "X2"

            p <- ggplot2::ggplot(dat) +
              ggplot2::labs(title = paste0(title, " (Waiting for data)"), x = x_name, y = y_name) +
              ggplot2::theme_minimal() +
              ggplot2::theme(aspect.ratio = 1, legend.position = "none")

            if (nrow(dat) == 0) {
              p <- p + ggplot2::xlim(-4, 4) + ggplot2::ylim(-4, 4)
            }

            if (nrow(dat) > 0) {
              p <- p + ggplot2::geom_point(ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]], color = .data$Sim), size = 3) +
                ggplot2::scale_color_manual(values = color_palette(), drop = FALSE)
            }

            if (!is.null(zoom_xlim()) && !is.null(zoom_ylim())) {
              p <- p + ggplot2::coord_cartesian(xlim = zoom_xlim(), ylim = zoom_ylim(), expand = FALSE)
            }

            return(p)
          }

          state <- comparison_state()
          if (is.null(state) || is.null(state$models[[m]])) {
            return(
              ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0, y = 0, label = "Model failed to fit") +
                ggplot2::labs(title = title) +
                ggplot2::theme_minimal() +
                ggplot2::theme(aspect.ratio = 1)
            )
          }

          cb_mod <- state$models[[m]]$model

          # Determine slice_x/slice_y if in 2D Slice mode
          use_slice_x <- NULL
          use_slice_y <- NULL
          use_proj <- current_projection()
          use_proj_info <- current_projection_info()
          if (ncol(dat) > 3 && isTRUE(input$hd_view_mode == "slice")) {
            use_proj <- NULL
            use_proj_info <- NULL
            use_slice_x <- input$slice_x
            use_slice_y <- input$slice_y
          }

          tryCatch(
            {
              p <- withCallingHandlers(
                create_boundary_plot(
                  cb_mod = cb_mod,
                  data = dat,
                  title = title,
                  class_levels = class_choices(),
                  class_colors = color_palette(),
                  proj_matrix = use_proj,
                  proj_info = use_proj_info,
                  zoom_x = zoom_xlim(),
                  zoom_y = zoom_ylim(),
                  show_probs = isTRUE(input$show_probs),
                  resolution = if (!is.null(input$grid_resolution)) input$grid_resolution else 100,
                  predict_args = state$models[[m]]$predict_args,
                  n_outliers = nrow(injected_outliers()),
                  highlight_outliers = isTRUE(input$highlight_outliers),
                  slice_x = use_slice_x,
                  slice_y = use_slice_y
                ),
                warning = function(w) {
                  msg <- conditionMessage(w)
                  if (grepl("automatically imputed|projects points flat", msg, ignore.case = TRUE)) {
                    invokeRestart("muffleWarning")
                  }
                }
              )

              stk <- shiny::isolate(active_stroke_data())
              if (!is.null(stk) && nrow(stk) > 0 && identical(draw_stroke_plot_id(), safe_id)) {
                feat_cols <- setdiff(colnames(dat), "Sim")
                x_name <- if (length(feat_cols) >= 1) feat_cols[1] else "X1"
                y_name <- if (length(feat_cols) >= 2) feat_cols[2] else "X2"

                # Overlay active stroke preview using selected class color
                preview_color <- tryCatch(color_palette()[[shiny::isolate(input$draw_class)]], error = function(e) "grey50")
                if (is.null(preview_color) || is.na(preview_color)) preview_color <- "grey50"
                p <- p + ggplot2::geom_point(data = stk, ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]]), color = preview_color, size = 1.5, alpha = 0.5)
              }

              p
            },
            error = function(e) {
              # Handle degenerate data gracefully.
              err_msg <- clean_err_msg(e$message)
              err_text <- paste(strwrap(paste("Error:", err_msg), width = 40), collapse = "\n")

              ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0, y = 0, label = err_text) +
                ggplot2::labs(title = title) +
                ggplot2::theme_minimal() +
                ggplot2::theme(aspect.ratio = 1)
            }
          )
        })
      })
    })

    metrics_df <- shiny::reactive({
      state <- comparison_state()
      if (is.null(state) || length(state$models) == 0) {
        return(NULL)
      }

      dat <- state$train_data
      true_labels <- dat$Sim

      metrics_list <- lapply(names(state$models), function(m) {
        preds <- state$models[[m]]$predictions$class

        if (is.null(preds) || length(preds) != length(true_labels)) {
          return(data.frame(
            Model = m,
            `Training Accuracy` = NA,
            `Training Kappa` = NA,
            `Training Error` = NA,
            check.names = FALSE
          ))
        }

        acc <- sum(preds == true_labels) / length(true_labels)
        err <- 1 - acc

        levs <- class_choices()
        if (is.null(levs)) levs <- unique(c(as.character(preds), as.character(true_labels)))
        tab <- table(factor(preds, levels = levs), factor(true_labels, levels = levs))
        row_sums <- rowSums(tab)
        col_sums <- colSums(tab)
        pe <- sum((row_sums * col_sums) / sum(tab)^2)
        kappa <- if (is.nan(pe) || pe == 1) 1 else (acc - pe) / (1 - pe)

        df_row <- data.frame(
          Model = m,
          `Training Accuracy` = round(acc, 4),
          `Training Kappa` = round(kappa, 4),
          `Training Error` = round(err, 4),
          check.names = FALSE
        )

        test_preds <- state$models[[m]]$test_predictions$class
        if (!is.null(test_preds)) {
          test_true <- shiny::isolate(current_test_data()$Sim)
          if (!is.null(test_true) && length(test_preds) == length(test_true)) {
            test_acc <- sum(test_preds == test_true) / length(test_true)
            test_err <- 1 - test_acc
            df_row$`Test Error` <- round(test_err, 4)
          }
        }

        df_row
      })

      do.call(rbind, metrics_list)
    })

    output$metrics_table <- DT::renderDataTable({
      res <- metrics_df()
      if (is.null(res)) {
        return(NULL)
      }

      DT::datatable(
        res,
        options = list(
          dom = "t",
          paging = FALSE,
          searching = FALSE,
          ordering = TRUE
        ),
        rownames = FALSE,
        class = "cell-border stripe hover"
      )
    })

    output$vis_info_ui <- DT::renderDataTable({
      res_val <- if (!is.null(input$grid_resolution)) input$grid_resolution else 100
      tot_points <- res_val * res_val
      render_mode <- "raster"
      proj_active <- !is.null(current_projection())

      dat <- combined_training_data()
      range_str <- if (!is.null(dat) && nrow(dat) > 0) {
        feat_cols <- setdiff(colnames(dat), "Sim")
        if (length(feat_cols) > 0) {
          paste(sapply(feat_cols[seq_len(min(2, length(feat_cols)))], function(f) {
            r <- range(dat[[f]], na.rm = TRUE)
            sprintf("%s: [%.1f, %.1f]", f, r[1], r[2])
          }), collapse = " | ")
        } else {
          "N/A"
        }
      } else {
        "N/A"
      }

      df <- data.frame(
        Metric = c("Grid Resolution", "Total Grid Points", "Feature Ranges", "Rendering Engine", "Space"),
        Value = c(
          sprintf("%d x %d", res_val, res_val),
          format(tot_points, big.mark = ","),
          range_str,
          render_mode,
          if (proj_active) "High-Dimensional Projection" else "Native 2D Feature Space"
        ),
        stringsAsFactors = FALSE
      )

      DT::datatable(
        df,
        rownames = FALSE,
        colnames = rep("", ncol(df)),
        options = list(
          dom = "t",
          bSort = FALSE,
          paging = FALSE,
          language = list(emptyTable = "Waiting for data...")
        ),
        selection = "none"
      )
    })

    shiny::observeEvent(input$open_export_wizard, {
      shiny::showModal(shiny::modalDialog(
        title = "Export Wizard",
        shiny::checkboxGroupInput(
          "export_includes",
          "Include in Export:",
          choices = c("Data", "Fitted Models", "Plots", "Grid Predictions", "Performance Metrics", "Reproduce Script", "Configuration", "Session Info"),
          selected = c("Data", "Fitted Models", "Plots", "Performance Metrics", "Reproduce Script")
        ),
        shiny::conditionalPanel(
          condition = "input.export_includes.indexOf('Plots') > -1",
          shiny::checkboxGroupInput(
            "export_format",
            "Plot Format:",
            choices = c("PNG", "PDF"),
            selected = c("PNG")
          ),
          shiny::conditionalPanel(
            condition = "input.export_format.indexOf('PNG') > -1",
            shiny::selectInput("export_dpi", "PNG Resolution (DPI):", choices = c("150", "300", "600"), selected = "300")
          )
        ),
        footer = shiny::tagList(
          shiny::modalButton("Cancel"),
          shiny::downloadButton("export_download", "Download ZIP", class = "btn-success")
        ),
        size = "s"
      ))
    })

    shiny::observeEvent(input$export_includes,
      {
        shiny::req(input$export_includes)
        includes <- input$export_includes
        if ("Reproduce Script" %in% includes) {
          required <- c("Data", "Fitted Models")
          missing <- setdiff(required, includes)
          if (length(missing) > 0) {
            shiny::updateCheckboxGroupInput(
              session = session,
              inputId = "export_includes",
              selected = unique(c(includes, required))
            )
          }
        }
      },
      ignoreInit = TRUE
    )

    output$export_download <- shiny::downloadHandler(
      filename = function() {
        paste("classbound_export_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip", sep = "")
      },
      content = function(file) {
        shiny::removeModal()

        temp_dir <- tempdir()
        export_dir <- file.path(temp_dir, paste0("export_", as.integer(Sys.time())))
        dir.create(export_dir)

        dat <- combined_training_data()
        state <- comparison_state()
        metrics <- metrics_df()
        includes <- input$export_includes

        if ("Data" %in% includes && !is.null(dat) && nrow(dat) > 0) {
          export_data_csv(dat, file.path(export_dir, "data.csv"))
          test_dat <- shiny::isolate(current_test_data())
          if (!is.null(test_dat) && nrow(test_dat) > 0) {
            export_data_csv(test_dat, file.path(export_dir, "test_data.csv"))
          }
        }

        if (!is.null(state) && length(state$models) > 0) {
          if ("Fitted Models" %in% includes) {
            mod_list <- lapply(state$models, function(x) x$model)
            export_models_rds(mod_list, file.path(export_dir, "models.rds"))
          }

          if ("Plots" %in% includes) {
            plots <- lapply(names(state$models), function(m) {
              title <- switch(m,
                "rpart" = "Decision Tree (rpart)",
                "PPtreeViz" = "PPtreeViz",
                "PPtreeExtclass" = "PPtreeExtclass",
                "PPtreeExt_split" = "PPtreeExt_split",
                "randomForest" = "Random Forest",
                m
              )
              tryCatch(
                {
                  # Determine slice params for export
                  exp_proj <- current_projection()
                  exp_proj_info <- current_projection_info()
                  exp_sx <- NULL
                  exp_sy <- NULL
                  if (ncol(dat) > 3 && isTRUE(input$hd_view_mode == "slice")) {
                    exp_proj <- NULL
                    exp_proj_info <- NULL
                    exp_sx <- input$slice_x
                    exp_sy <- input$slice_y
                  }
                  p_tmp <- withCallingHandlers(
                    create_boundary_plot(
                      cb_mod = state$models[[m]]$model,
                      data = dat,
                      title = title,
                      class_levels = class_choices(),
                      class_colors = color_palette(),
                      proj_matrix = exp_proj,
                      proj_info = exp_proj_info,
                      zoom_x = zoom_xlim(),
                      zoom_y = zoom_ylim(),
                      show_probs = isTRUE(input$show_probs),
                      resolution = if (!is.null(input$grid_resolution)) input$grid_resolution else 100,
                      predict_args = state$models[[m]]$predict_args,
                      n_outliers = nrow(injected_outliers()),
                      highlight_outliers = isTRUE(input$highlight_outliers),
                      slice_x = exp_sx,
                      slice_y = exp_sy
                    ),
                    warning = function(w) {
                      msg <- conditionMessage(w)
                      if (grepl("automatically imputed|projects points flat", msg, ignore.case = TRUE)) {
                        invokeRestart("muffleWarning")
                      }
                    }
                  )

                  stk <- active_stroke_data()
                  if (!is.null(stk) && nrow(stk) > 0) {
                    feat_cols <- setdiff(colnames(dat), "Sim")
                    x_name <- if (length(feat_cols) >= 1) feat_cols[1] else "X1"
                    y_name <- if (length(feat_cols) >= 2) feat_cols[2] else "X2"
                    p_tmp <- p_tmp + ggplot2::geom_point(data = stk, ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]]), color = "black", size = 1.5, alpha = 0.6)
                  }
                  p_tmp
                },
                error = function(e) NULL
              )
            })

            names(plots) <- names(state$models)
            plots <- plots[!sapply(plots, is.null)]

            if (length(plots) > 0) {
              if ("PDF" %in% input$export_format) {
                export_plots_pdf(plots, file.path(export_dir, "plots.pdf"))
              }
              if ("PNG" %in% input$export_format) {
                png_dir <- file.path(export_dir, "plots")
                export_plots_png(plots, png_dir, dpi = input$export_dpi)
              }
            }
          }
        }

        if ("Performance Metrics" %in% includes && !is.null(metrics)) {
          export_metrics_csv(metrics, file.path(export_dir, "metrics.csv"))
        }

        if ("Session Info" %in% includes) {
          utils::capture.output(utils::sessionInfo(), file = file.path(export_dir, "session_info.txt"))
        }

        # Export grid predictions for each model as CSV.
        if ("Grid Predictions" %in% includes && !is.null(state) && length(state$models) > 0) {
          grid_dir <- file.path(export_dir, "grid_predictions")
          dir.create(grid_dir)
          res_val <- if (!is.null(input$grid_resolution)) input$grid_resolution else 100
          for (m in names(state$models)) {
            safe_name <- gsub("[^A-Za-z0-9_.-]", "_", m)
            export_grid_csv(
              cb_mod = state$models[[m]]$model,
              data = dat,
              resolution = res_val,
              proj_matrix = if (!is.null(dat) && ncol(dat) > 3 && isTRUE(input$hd_view_mode == "slice")) NULL else current_projection(),
              proj_info = if (!is.null(dat) && ncol(dat) > 3 && isTRUE(input$hd_view_mode == "slice")) NULL else current_projection_info(),
              predict_args = state$models[[m]]$predict_args,
              slice_x = if (!is.null(dat) && ncol(dat) > 3 && isTRUE(input$hd_view_mode == "slice")) input$slice_x else NULL,
              slice_y = if (!is.null(dat) && ncol(dat) > 3 && isTRUE(input$hd_view_mode == "slice")) input$slice_y else NULL,
              file = file.path(grid_dir, paste0(safe_name, ".csv"))
            )
          }
        }

        # Export current UI configuration as JSON for reproducibility.
        if ("Configuration" %in% includes) {
          all_inputs <- shiny::isolate(shiny::reactiveValuesToList(input))

          base_keys <- c("data_mode", "show_probs", "grid_resolution", "color_palette", "highlight_outliers", "outlier_class", "outlier_count", "outlier_magnitude", "draw_class")
          active_mode <- input$data_mode
          mode_keys <- character(0)

          if (isTRUE(active_mode == "Simulate Data")) {
            mode_keys <- c(
              "sim_engine", "sim_seed", "sim_noise", "sim_n_classes",
              paste0("sim_mean_", 1:10), paste0("sim_sd_", 1:10),
              paste0("sim_cor_", 1:10), paste0("sim_n_", 1:10),
              "sim_mixsim_k", "sim_mixsim_p", "sim_mixsim_omega", "sim_mixsim_n"
            )
          } else if (isTRUE(active_mode == "Draw Data")) {
            mode_keys <- c(
              "interaction_mode", "draw_total_classes",
              "brush_size", "brush_spread"
            )
          }

          hd_keys <- character(0)
          if (!is.null(dat) && ncol(dat) > 3) {
            hd_keys <- c("hd_view_mode")
            if (isTRUE(input$hd_view_mode == "projection")) {
              hd_keys <- c(hd_keys, "tour_var", "tour_angle")
            } else if (isTRUE(input$hd_view_mode == "slice")) {
              hd_keys <- c(hd_keys, "slice_x", "slice_y")
            }
          }

          model_keys <- c()
          selected_models <- if (!is.null(state)) names(state$models) else c()

          if ("rpart" %in% selected_models) model_keys <- c(model_keys, "rpart_cp")
          if ("randomForest" %in% selected_models) model_keys <- c(model_keys, "rf_ntree", "rf_mtry")
          if ("PPtreeViz" %in% selected_models) model_keys <- c(model_keys, "rule", "pp_method")
          if ("PPtreeExtclass" %in% selected_models || "PPtreeExt_split" %in% selected_models) model_keys <- c(model_keys, "stop")
          if ("ppforest2" %in% selected_models) model_keys <- c(model_keys, "pprf_size", "pprf_lambda")

          all_keys <- unique(c(base_keys, mode_keys, hd_keys, model_keys))
          safe_inputs <- all_inputs[intersect(names(all_inputs), all_keys)]
          safe_inputs <- Filter(function(x) is.numeric(x) || is.character(x) || is.logical(x), safe_inputs)

          config_vals <- list(
            ui_state = safe_inputs,
            timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
            r_version = paste(R.version$major, R.version$minor, sep = "."),
            classbound_version = as.character(utils::packageVersion("classbound"))
          )

          sim_conf <- applied_sim_config()
          if (!is.null(sim_conf)) {
            if (!is.null(sim_conf$seed) && !is.na(sim_conf$seed)) {
              sim_conf$reproducible_from_settings <- TRUE
            } else {
              sim_conf$reproducible_from_settings <- FALSE
              sim_conf$reproducibility_note <- "No random seed was provided. Regenerating the simulation from these settings will produce a new random realization. To reproduce the generated result exactly, export the Reproduce Script."
            }
            config_vals$simulation_settings <- sim_conf
          }

          if (!is.null(state) && length(state$models) > 0) {
            config_vals$fitted_models <- names(state$models)
          }
          export_config_json(config_vals, file.path(export_dir, "config.json"))
        }

        # Export a self-contained R script to reproduce the results.
        if ("Reproduce Script" %in% includes && !is.null(state) && length(state$models) > 0) {
          proj <- if (!is.null(dat) && ncol(dat) > 3 && isTRUE(input$hd_view_mode == "slice")) NULL else current_projection()
          proj_info <- if (!is.null(dat) && ncol(dat) > 3 && isTRUE(input$hd_view_mode == "slice")) NULL else current_projection_info()
          if (!is.null(proj)) {
            saveRDS(
              list(basis = proj, center = proj_info$center, scale = proj_info$scale),
              file.path(export_dir, "projection.rds")
            )
          }
          export_reproduce_script(
            model_names = names(state$models),
            has_projection = !is.null(proj),
            slice_x = if (!is.null(dat) && ncol(dat) > 3 && isTRUE(input$hd_view_mode == "slice")) input$slice_x else NULL,
            slice_y = if (!is.null(dat) && ncol(dat) > 3 && isTRUE(input$hd_view_mode == "slice")) input$slice_y else NULL,
            resolution = if (!is.null(input$grid_resolution)) input$grid_resolution else 100,
            show_probs = isTRUE(input$show_probs),
            zoom_x = zoom_xlim(),
            zoom_y = zoom_ylim(),
            file = file.path(export_dir, "reproduce.R")
          )
        }

        owd <- setwd(export_dir)
        on.exit(setwd(owd))
        utils::zip(zipfile = file, files = list.files(export_dir, recursive = TRUE), extras = "-q")
      }
    )
  }

  shiny::shinyApp(ui = ui, server = server)
}

#' Render a decision boundary plot for the Shiny app
#'
#' @param cb_mod A classbound_model object.
#' @param data The training data to overlay.
#' @param title The title of the plot.
#' @param class_levels The factor levels for the classes.
#' @param proj_matrix Optional projection matrix.
#' @param proj_info Optional list with center and scale for projection.
#' @param zoom_x Optional x-axis limits.
#' @param zoom_y Optional y-axis limits.
#' @param show_probs Whether to show probability gradients.
#' @param resolution The grid resolution.
#' @param slice_x Optional. Feature name to use for the X axis in 2D Slice mode.
#' @param slice_y Optional. Feature name to use for the Y axis in 2D Slice mode.
#' @return A ggplot object.
#' @keywords internal
create_boundary_plot <- function(cb_mod, data, title, class_levels = NULL, class_colors = NULL, proj_matrix = NULL, proj_info = NULL, zoom_x = NULL, zoom_y = NULL, show_probs = FALSE, resolution = 100, predict_args = list(), n_outliers = 0, highlight_outliers = TRUE, slice_x = NULL, slice_y = NULL) {
  if (!is.null(proj_matrix)) {
    proj_list <- list(basis = proj_matrix, center = proj_info$center, scale = proj_info$scale)
    x_mat <- as.matrix(data[, rownames(proj_matrix)])
    if (!is.null(proj_info$center)) x_mat <- sweep(x_mat, 2, proj_info$center, "-")
    if (!is.null(proj_info$scale)) x_mat <- sweep(x_mat, 2, proj_info$scale, "/")
    z_mat <- x_mat %*% proj_matrix

    r1 <- range(z_mat[, 1])
    r2 <- range(z_mat[, 2])
    pad1 <- max(diff(r1) * 0.06, 0.5)
    pad2 <- max(diff(r2) * 0.06, 0.5)

    range_list <- list()
    if (!is.null(zoom_x) && !is.null(zoom_y)) {
      range_list[["PC1"]] <- c(min(r1[1] - pad1, min(zoom_x)), max(r1[2] + pad1, max(zoom_x)))
      range_list[["PC2"]] <- c(min(r2[1] - pad2, min(zoom_y)), max(r2[2] + pad2, max(zoom_y)))
    } else {
      range_list[["PC1"]] <- r1 + c(-pad1, pad1)
      range_list[["PC2"]] <- r2 + c(-pad2, pad2)
    }

    cb_bound <- boundary_compute(cb_mod, feature_range = range_list, resolution = resolution, projection = proj_list, predict_args = predict_args)
    x_col_label <- "PC1"
    y_col_label <- "PC2"
  } else {
    feat_cols <- setdiff(colnames(data), "Sim")
    # Use slice_x/slice_y if provided; otherwise default to first two features
    x_name <- if (!is.null(slice_x) && slice_x %in% feat_cols) slice_x else feat_cols[1]
    y_name <- if (!is.null(slice_y) && slice_y %in% feat_cols) slice_y else feat_cols[2]

    r1 <- range(data[[x_name]])
    r2 <- range(data[[y_name]])
    pad1 <- max(diff(r1) * 0.06, 0.5)
    pad2 <- max(diff(r2) * 0.06, 0.5)

    range_list <- list()
    if (!is.null(zoom_x) && !is.null(zoom_y)) {
      range_list[[x_name]] <- c(min(r1[1] - pad1, min(zoom_x)), max(r1[2] + pad1, max(zoom_x)))
      range_list[[y_name]] <- c(min(r2[1] - pad2, min(zoom_y)), max(r2[2] + pad2, max(zoom_y)))
    } else {
      range_list[[x_name]] <- r1 + c(-pad1, pad1)
      range_list[[y_name]] <- r2 + c(-pad2, pad2)
    }

    cb_bound <- boundary_compute(cb_mod, feature_range = range_list, resolution = resolution, predict_args = predict_args)
    x_col_label <- x_name
    y_col_label <- y_name
  }

  # Lock factor levels if provided to prevent color shifting
  if (!is.null(class_levels)) {
    data$Sim <- factor(data$Sim, levels = class_levels)
    cb_bound$boundary_data$prediction <- factor(cb_bound$boundary_data$prediction, levels = class_levels)
  }

  # Tag outlier rows for visual differentiation (only at plot time, not during model fitting)
  if (n_outliers > 0) {
    data$is_outlier <- FALSE
    data$is_outlier[(nrow(data) - n_outliers + 1):nrow(data)] <- TRUE
  }

  # Render plot with explicit color palette to stay in sync with the UI legend
  p <- plot_boundary(cb_bound, obs_data = data, x_col = x_col_label, y_col = y_col_label, true_label = "Sim", show_gradient = show_probs, colors = class_colors, highlight_outliers = highlight_outliers, xlim = zoom_x, ylim = zoom_y) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(aspect.ratio = 1, legend.position = "none")

  p
}

#' Generate a deterministic extreme outlier based on dataset bounds
#'
#' @param data A data frame containing the features and target column.
#' @param class_label The class label to assign to the outlier.
#' @param magnitude A numeric value indicating how far outside the bounding box to place the outlier.
#' @param target_col The name of the target column in the data.
#' @param index An integer index used to cycle through bounding box corners.
#' @return A one-row data frame with the generated outlier.
#' @keywords internal
generate_outlier <- function(data, class_label, magnitude, target_col = "Sim", index = 1) {
  outlier <- data[1, , drop = FALSE]
  outlier[[target_col]] <- class_label

  feat_cols <- setdiff(colnames(data), target_col)
  num_cols <- names(which(sapply(data[feat_cols], is.numeric)))

  class_offset <- suppressWarnings(sum(utf8ToInt(as.character(class_label))))
  corner <- (index + class_offset) %% 4

  target_data <- data[data[[target_col]] == class_label, , drop = FALSE]
  if (nrow(target_data) == 0) {
    target_data <- data
  }

  visual_cols <- if (length(num_cols) >= 2) num_cols[1:2] else num_cols

  apply_tukey <- function(col_name, is_x) {
    valid_vals <- target_data[[col_name]][!is.na(target_data[[col_name]])]

    # Degenerate 1: Target class is completely missing values for this feature
    if (length(valid_vals) == 0) {
      global_vals <- data[[col_name]][!is.na(data[[col_name]])]
      if (length(global_vals) == 0) {
        return(NA_real_)
      }
      base_val <- stats::median(global_vals)
      scale <- max(global_vals) - min(global_vals)
      if (scale == 0) {
        return(unname(base_val))
      }

      use_max <- if (is_x) (corner %in% c(0, 3)) else (corner %in% c(0, 1))
      if (use_max) {
        return(unname(base_val + magnitude * scale))
      }
      return(unname(base_val - magnitude * scale))
    }

    # Calculate target class Tukey stats
    q1 <- stats::quantile(valid_vals, 0.25)
    q3 <- stats::quantile(valid_vals, 0.75)
    iqr <- q3 - q1

    # Degenerate 2: Zero IQR in target class
    if (iqr == 0) {
      iqr <- max(valid_vals) - min(valid_vals)
    }

    # Degenerate 3: Zero variance (collinear) in target class or insufficient samples
    if (iqr == 0) {
      # Use global scale but originate from the class's constant value
      global_vals <- data[[col_name]][!is.na(data[[col_name]])]
      if (length(global_vals) > 0) {
        global_scale <- max(global_vals) - min(global_vals)
        if (global_scale > 0) {
          use_max <- if (is_x) (corner %in% c(0, 3)) else (corner %in% c(0, 1))
          base_val <- valid_vals[1]
          if (use_max) {
            return(unname(base_val + magnitude * global_scale))
          }
          return(unname(base_val - magnitude * global_scale))
        }
      }
      return(unname(valid_vals[1]))
    }

    # Valid Tukey
    use_max <- if (is_x) (corner %in% c(0, 3)) else (corner %in% c(0, 1))
    if (use_max) {
      return(unname(q3 + magnitude * iqr))
    }
    return(unname(q1 - magnitude * iqr))
  }

  if (length(visual_cols) == 2) {
    mat <- as.matrix(target_data[, visual_cols])
    mat_clean <- mat[stats::complete.cases(mat), , drop = FALSE]

    # Require genuinely valid, finite, invertible 2D covariance matrix (n > 2 for full rank)
    if (nrow(mat_clean) > 2) {
      mu <- colMeans(mat_clean)
      Sigma <- tryCatch(stats::cov(mat_clean), error = function(e) NULL)

      # Check if singular, zero variance, or invalid
      if (is.null(Sigma) || any(is.na(Sigma)) || any(diag(Sigma) == 0) || abs(det(Sigma)) < 1e-10) {
        outlier[[visual_cols[1]]] <- apply_tukey(visual_cols[1], TRUE)
        outlier[[visual_cols[2]]] <- apply_tukey(visual_cols[2], FALSE)
      } else {
        eig <- eigen(Sigma)
        val <- eig$values
        vec <- eig$vectors

        # Directions: 0 = +major, 1 = -major, 2 = +minor, 3 = -minor
        if (corner %in% c(0, 1)) {
          dir <- vec[, 1]
          scale <- sqrt(val[1])
          sign_mult <- if (corner == 0) 1 else -1
        } else {
          dir <- vec[, 2]
          scale <- sqrt(val[2])
          sign_mult <- if (corner == 2) 1 else -1
        }

        pt <- mu + sign_mult * magnitude * scale * dir
        outlier[[visual_cols[1]]] <- unname(pt[1])
        outlier[[visual_cols[2]]] <- unname(pt[2])
      }
    } else {
      # Insufficient data for 2D covariance -> fallback to Tukey
      outlier[[visual_cols[1]]] <- apply_tukey(visual_cols[1], TRUE)
      outlier[[visual_cols[2]]] <- apply_tukey(visual_cols[2], FALSE)
    }
  } else if (length(visual_cols) == 1) {
    outlier[[visual_cols[1]]] <- apply_tukey(visual_cols[1], TRUE)
  }

  # Non-visualized numeric columns
  for (col in setdiff(num_cols, visual_cols)) {
    outlier[[col]] <- stats::median(target_data[[col]], na.rm = TRUE)
  }

  # Categorical columns
  for (col in setdiff(feat_cols, num_cols)) {
    freqs <- table(target_data[[col]])
    if (length(freqs) == 0 || max(freqs) == 0) {
      freqs <- table(data[[col]])
    }
    if (length(freqs) > 0) {
      mode_val <- names(freqs)[which.max(freqs)]
      if (is.factor(data[[col]])) {
        outlier[[col]] <- factor(mode_val, levels = levels(data[[col]]))
      } else {
        outlier[[col]] <- mode_val
      }
    } else {
      outlier[[col]] <- NA
    }
  }

  return(outlier)
}
