#' Shiny app to compare PPtree, PPtreeExt and rpart boundaries in 2D with different simulation scenarios
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
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' must be installed to use explorapp().", call. = FALSE)
  }

  if (!requireNamespace("MixSim", quietly = TRUE)) {
    stop("Package 'MixSim' must be installed to use explorapp().", call. = FALSE)
  }
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("Package 'DT' must be installed to use explorapp().", call. = FALSE)
  }

  if (!is.null(data)) {
    if (is.null(target_col) || !(target_col %in% colnames(data))) {
      stop("Please provide a valid 'target_col' that exists in 'data'.", call. = FALSE)
    }
  }

  # UI to Package API Mapping
  app_methods <- list(
    "PPtreeViz"       = list(fn = PPtreeViz::PPTreeclass, args = list(PPmethod = "LDA"), supports_prob = FALSE),
    "Rpart"           = list(fn = rpart::rpart, args = list(), supports_prob = TRUE),
    "PPtreeExt_split" = list(fn = PPtreeExt::PPtreeExt_split, args = list(PPmethod = "LDA"), supports_prob = FALSE),
    "PPtreeExtclass"  = list(fn = PPtreeExt::PPtreeExtclass, args = list(PPmethod = "LDA"), supports_prob = FALSE),
    "RandomForest"    = list(fn = randomForest::randomForest, args = list(), supports_prob = TRUE),
    "PPforest"        = list(fn = function(formula, data, ...) {
      y_name <- all.vars(formula[[2]])
      PPforest::PPforest(data = as.data.frame(data), y = y_name, ...)
    }, args = list(PPmethod = "LDA", size.tr = 1, size.p = 1), supports_prob = FALSE)
  )

  predict_args <- list(
    "PPtreeViz"       = function(ru) list(Rule = ru),
    "Rpart"           = function(...) list(),
    "PPtreeExt_split" = function(ru) list(Rule = ru),
    "PPtreeExtclass"  = function(...) list(),
    "RandomForest"    = function(...) list(),
    "PPforest"        = function(...) list()
  )

  # Dynamically add Tidymodels presets
  if (requireNamespace("parsnip", quietly = TRUE)) {
    tm_friendly <- c(
      "rpart" = "Decision Tree (rpart)",
      "randomForest" = "Random Forest",
      "kernlab" = "SVM (kernlab)",
      "nnet" = "Neural Net (nnet)"
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
    shiny::tags$head(shiny::tags$style(shiny::HTML(".col-sm-4 { max-height: 90vh; overflow-y: auto; }"))),
    shiny::titlePanel("Classbound Exploration & Comparison"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::wellPanel(
          shiny::radioButtons("data_mode", "Data Mode", choices = if(!is.null(data)) c("Import Data", "Simulate Data", "Draw Data") else c("Simulate Data", "Draw Data"), selected = if (!is.null(data)) "Import Data" else "Simulate Data"),

          shiny::conditionalPanel(
            condition = "input.data_mode == 'Simulate Data'",
            shiny::radioButtons("sim_engine", "Simulation Engine", choices = c("Multivariate Normal (MVN)" = "mvn", "MixSim" = "mixsim"))
          ),

          shiny::hr(),
          shiny::radioButtons("interaction_mode", "Interaction Mode", choices = c("Navigate", "Draw Point", "Draw Cluster")),
          shiny::helpText("Tip: Brush on any plot to zoom in, and double-click to reset the view."),

          shiny::conditionalPanel(
            condition = "input.data_mode == 'Simulate Data'",
            shiny::hr(),
            shiny::conditionalPanel(
              condition = "input.sim_engine == 'mvn'",
              shiny::numericInput("sim_n_classes", "Number of Classes", value = 3, min = 2, max = 10),
              shiny::helpText("Tip: Class parameters appear dynamically based on the number of classes chosen above."),
              shiny::uiOutput("sim_params_ui"),
              shiny::div(style = "display: flex; gap: 10px; flex-wrap: wrap;",
                shiny::actionButton("sim_do", "Generate Data", class = "btn-primary"),
                shiny::actionButton("clone_to_draw_sim_mvn", "Clone to Draw Canvas", icon = shiny::icon("copy"), class = "btn-warning")
              )
            ),

            shiny::conditionalPanel(
              condition = "input.sim_engine == 'mixsim'",
              shiny::numericInput("sim_mixsim_k", "Number of Classes (K)", value = 3, min = 2, max = 10),
              shiny::numericInput("sim_mixsim_p", "Dimensions (p)", value = 2, min = 2, max = 10),
              shiny::numericInput("sim_mixsim_omega", "Max Overlap (MaxOmega)", value = 0.05, min = 0.01, max = 0.5, step = 0.01),
              shiny::numericInput("sim_mixsim_n", "Sample Size", value = 300, min = 10, step = 100),
              shiny::div(style = "display: flex; gap: 10px; flex-wrap: wrap;",
                shiny::actionButton("sim_mixsim_do", "Generate Data", class = "btn-primary"),
                shiny::actionButton("clone_to_draw_sim_mix", "Clone to Draw Canvas", icon = shiny::icon("copy"), class = "btn-warning")
              )
            )
          ),

          shiny::conditionalPanel(
            condition = "input.data_mode == 'Draw Data'",
            shiny::hr(),
            shiny::p("Click or brush on any plot to add data."),
            shiny::numericInput("draw_total_classes", "Total Classes", value = 3, min = 1, max = 26, step = 1),
            shiny::radioButtons("draw_class", "Active Class", choices = c("Class 1", "Class 2", "Class 3")),
            shiny::div(
              style = "display: flex; gap: 10px; align-items: baseline; margin-bottom: 15px;",
              shiny::textInput("new_class_name", label = NULL, placeholder = "New class name...", width = "100%"),
              shiny::actionButton("add_class", "Add Class", class = "btn-info")
            ),
            shiny::numericInput("brush_size", "Brush Density (points/brush)", value = 20, min = 1),
            shiny::div(
              style = "display: flex; gap: 10px; margin-bottom: 15px; flex-wrap: wrap;",
              shiny::actionButton("undo_draw", "Undo Last", icon = shiny::icon("undo"), class = "btn-default"),
              shiny::actionButton("clear", "Clear Canvas", class = "btn-danger")
            ),
            shiny::h5("Drawn Points", style = "margin-top: 15px; font-weight: bold;"),
            DT::dataTableOutput("drawn_points_table")
          ),
          shiny::conditionalPanel(
            condition = "input.data_mode == 'Import Data'",
            shiny::p("Using dataset provided via console."),
            shiny::hr(),
            shiny::actionButton("clone_to_draw_imp", "Clone to Draw Canvas", icon = shiny::icon("copy"), class = "btn-warning")
          )
        ),
        shiny::wellPanel(
          shiny::tags$details(
            shiny::tags$summary("Model Configuration", style = "display: list-item; font-size: 18px; font-weight: 500; cursor: pointer; margin-bottom: 10px;"),
              shiny::checkboxGroupInput(
                "selected_models",
                "Models to Compare",
                choices = names(app_methods),
                selected = c("RandomForest", "Rpart")
              ),
            shiny::uiOutput("prob_surface_ui"),
            shiny::sliderInput("grid_resolution", "Grid Resolution", min = 50, max = 300, value = 100, step = 25),
            shiny::conditionalPanel(
              condition = "input.selected_models && (input.selected_models.indexOf('PPtreeViz') > -1 || input.selected_models.indexOf('PPtreeExtclass') > -1 || input.selected_models.indexOf('PPtreeExt_split') > -1)",
              shiny::hr(),
              shiny::selectInput("rule", "PPtree: Projection Pursuit Rule", choices = 1:8, selected = 1)
            ),
            shiny::conditionalPanel(
              condition = "input.selected_models && input.selected_models.indexOf('PPtreeExtclass') > -1",
              shiny::numericInput("stop", "PPtreeExt: Stopping Rule", value = 4, min = 1)
            ),
            shiny::conditionalPanel(
              condition = "input.selected_models && input.selected_models.indexOf('Rpart') > -1",
              shiny::hr(),
              shiny::numericInput("rpart_cp", "Rpart: Complexity Parameter", value = 0.01, min = 0, step = 0.01)
            ),
            shiny::conditionalPanel(
              condition = "input.selected_models && input.selected_models.indexOf('RandomForest') > -1",
              shiny::hr(),
              shiny::numericInput("rf_ntree", "Random Forest: Number of Trees", value = 500, min = 10, step = 50)
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
            shiny::tags$summary("Outlier Injection", style = "display: list-item; font-size: 18px; font-weight: 500; cursor: pointer; margin-bottom: 10px;"),
            shiny::helpText("Add extreme outliers to see how different models react to them."),
            shiny::selectInput("outlier_class", "Outlier Class", choices = c("Random", "Class 1", "Class 2", "Class 3")),
            shiny::sliderInput("outlier_magnitude", "Outlier Magnitude", min = 2, max = 15, value = 5, step = 0.5),
            shiny::numericInput("outlier_count", "Number of Outliers", value = 1, min = 1, max = 20, step = 1),
            shiny::div(
              style = "display: flex; gap: 10px;",
              shiny::actionButton("inject_outlier_btn", "Inject Outliers", class = "btn-danger"),
              shiny::actionButton("clear_outliers_btn", "Clear", class = "btn-default")
            )
          )
        ),
        shiny::wellPanel(
          shiny::tags$details(
            shiny::tags$summary("Import Workspace Models", style = "display: list-item; font-size: 18px; font-weight: 500; cursor: pointer; margin-bottom: 10px;"),
            shiny::helpText("This feature looks for fitted `workflow`, `model_fit`, or `model_spec` objects in your R Global Environment."),
            shiny::uiOutput("workspace_import_ui")
          )
        ),
        shiny::wellPanel(
          shiny::h4("Export Results", style = "margin-top: 0;"),
          shiny::actionButton("open_export_wizard", "Export Results...", class = "btn-success", icon = shiny::icon("download", lib = "font-awesome")),
          shiny::helpText(shiny::em("Customize export formats (PNG, PDF, etc.)"))
        ),
        shiny::uiOutput("tour_panel")
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
  server <- function(input, output) {
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

      initial_imp_classes <- as.character(unique(initial_imp_data$Sim))
      if (length(initial_imp_classes) == 0) {
        initial_imp_classes <- c("Class 1", "Class 2", "Class 3")
      }
    }

    initial_drawn_data <- data.frame(Sim = character(), X1 = numeric(), X2 = numeric())
    initial_drawn_classes <- c("Class 1", "Class 2", "Class 3")

    # Generate initial simulated dataset using simu_n helper
    initial_sim_data <- tryCatch({
      simu_n(
        means = list(c(-1, 0), c(1, 0)),
        covs = list(diag(c(1, 1)), diag(c(1, 1))),
        ns = c(50, 50),
        class_names = c("Class 1", "Class 2")
      )
    }, error = function(e) {
      data.frame(Sim = character(), X1 = numeric(), X2 = numeric())
    })
    initial_sim_classes <- c("Class 1", "Class 2")

    injected_outliers <- shiny::reactiveVal(data.frame())
    undo_history <- shiny::reactiveVal(list())

    combined_training_data <- shiny::reactive({
      dat <- current_data()
      inj <- injected_outliers()
      if (!is.null(inj) && nrow(inj) > 0) {
        dat <- rbind(dat, inj)
      }
      dat
    })

    mode_states <- shiny::reactiveValues(
      "Import Data" = list(data = if (!is.null(data)) initial_imp_data else NULL, classes = if (!is.null(data)) initial_imp_classes else NULL, outliers = data.frame()),
      "Draw Data" = list(data = initial_drawn_data, classes = initial_drawn_classes, outliers = data.frame()),
      "Simulate Data" = list(data = initial_sim_data, classes = initial_sim_classes, outliers = data.frame())
    )

    init_mode <- if (!is.null(data)) "Import Data" else "Simulate Data"
    previous_data_mode <- shiny::reactiveVal(init_mode)

    init_data <- if (!is.null(data)) initial_imp_data else initial_sim_data
    init_classes <- if (!is.null(data)) initial_imp_classes else initial_sim_classes

    current_data <- shiny::reactiveVal(init_data)
    class_choices <- shiny::reactiveVal(init_classes)
    zoom_xlim <- shiny::reactiveVal(NULL)
    zoom_ylim <- shiny::reactiveVal(NULL)

    shiny::observe({
      choices <- class_choices()
      if (length(choices) > 0) {
        choices <- c("Random", choices)
      }
      shiny::updateSelectInput(shiny::getDefaultReactiveDomain(), "outlier_class", choices = choices)
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

      new_pts <- do.call(rbind, lapply(1:count, function(i) {
        target_class <- input$outlier_class
        if (target_class == "Random") {
          avail_classes <- class_choices()
          if (length(avail_classes) > 0) {
            target_class <- sample(avail_classes, 1)
          }
        }

        # Add slight variation to magnitude relative to the chosen magnitude
        # This prevents outliers at the same corner from perfectly overlapping
        # We base the cycle on total_index so manual clicks also receive the variance
        total_index <- start_idx + i - 1
        cycle <- floor((total_index - 1) / 4)

        # Increase step_size to 10% for a more noticeable gap between cycles
        step_size <- input$outlier_magnitude * 0.1
        mag <- input$outlier_magnitude + (cycle * step_size)
        generate_extreme_outlier(cd, target_class, mag, target_col = "Sim", index = total_index)
      }))

      if (nrow(inj) == 0) {
        injected_outliers(new_pts)
      } else {
        injected_outliers(rbind(inj, new_pts))
      }
    })

    shiny::observeEvent(input$clear_outliers_btn, {
      injected_outliers(data.frame())
    })

    output$data_stats_ui <- shiny::renderUI({
      dat <- combined_training_data()
      if (is.null(dat) || nrow(dat) == 0) {
        return(NULL)
      }

      n_obs <- nrow(dat)
      dims <- ncol(dat) - 1
      classes <- class_choices()

      colors <- scales::hue_pal()(length(classes))
      legend_items <- lapply(seq_along(classes), function(i) {
        shiny::span(
          style = "display: inline-flex; align-items: center; margin-left: 10px;",
          shiny::tags$span(style = sprintf("display: inline-block; width: 12px; height: 12px; border-radius: 50%%; margin-right: 5px; background-color: %s;", colors[i])),
          classes[i]
        )
      })

      shiny::div(
        style = "background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; padding: 15px; margin-bottom: 20px; display: flex; justify-content: space-around; align-items: center;",
        shiny::div(shiny::tags$b("N Observations: "), n_obs),
        shiny::div(shiny::tags$b("Dimensions: "), dims),
        shiny::div(
          shiny::tags$b("Target Classes: "),
          shiny::div(style = "display: inline-flex;", legend_items)
        )
      )
    })

    output$is_high_dim <- shiny::reactive({
      ncol(combined_training_data()) > 3
    })
    shiny::outputOptions(output, "is_high_dim", suspendWhenHidden = FALSE)

    output$prob_surface_ui <- shiny::renderUI({
      req_models <- input$selected_models
      if (length(req_models) == 0) return(NULL)

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
            style = "opacity: 0.6; pointer-events: none;",
            title = "None of the selected classifiers support probability outputs. Falling back to hard decision boundaries.",
            shiny::checkboxInput("show_probs_disabled", "Show Probability Surface", value = current_val)
          )
        }
      )
    })

    shiny::observeEvent(input$data_mode, {
      zoom_xlim(NULL)
      zoom_ylim(NULL)

      old_mode <- previous_data_mode()
      new_mode <- input$data_mode

      # Save state to old mode
      if (!is.null(old_mode)) {
        mode_states[[old_mode]]$data <- current_data()
        mode_states[[old_mode]]$classes <- class_choices()
        mode_states[[old_mode]]$basis <- current_basis()
        mode_states[[old_mode]]$projection <- current_projection()
        mode_states[[old_mode]]$projection_info <- current_projection_info()
        mode_states[[old_mode]]$outliers <- injected_outliers()
      }

      # Load state from new mode
      if (!is.null(mode_states[[new_mode]]$data)) {
        current_data(mode_states[[new_mode]]$data)

        saved_outliers <- mode_states[[new_mode]]$outliers
        if (is.null(saved_outliers)) saved_outliers <- data.frame()
        injected_outliers(saved_outliers)

        class_choices(mode_states[[new_mode]]$classes)

        current_basis(mode_states[[new_mode]]$basis)
        current_projection(mode_states[[new_mode]]$projection)
        current_projection_info(mode_states[[new_mode]]$projection_info)

        # Reset UI element for drawing
        shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "draw_class", choices = mode_states[[new_mode]]$classes, selected = mode_states[[new_mode]]$classes[1])
      }

      previous_data_mode(new_mode)

      # Dynamically update interaction mode choices
      if (new_mode == "Draw Data") {
        shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "interaction_mode", choices = c("Navigate", "Draw Point", "Draw Cluster"))
      } else {
        shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "interaction_mode", choices = c("Navigate"))
      }
    })

    # Clone to Draw Canvas logic
    do_clone_to_draw <- function() {
      cd <- current_data()
      if (is.null(cd) || nrow(cd) == 0) {
        shiny::showNotification("No data to clone.", type = "warning")
        return()
      }

      # Copy data and outliers to Draw Data cache directly
      mode_states[["Draw Data"]]$data <- cd
      mode_states[["Draw Data"]]$classes <- class_choices()

      inj <- injected_outliers()
      if (is.null(inj)) inj <- data.frame()
      mode_states[["Draw Data"]]$outliers <- inj

      # Switch mode (the data_mode observer will handle loading the updated cache)
      shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "data_mode", selected = "Draw Data")
      shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "interaction_mode", selected = "Navigate")
      shiny::showNotification("Data and outliers copied to Draw Canvas!", type = "message")
    }

    shiny::observeEvent(input$clone_to_draw_sim_mvn, { do_clone_to_draw() })
    shiny::observeEvent(input$clone_to_draw_sim_mix, { do_clone_to_draw() })
    shiny::observeEvent(input$clone_to_draw_imp, { do_clone_to_draw() })

    shiny::observeEvent(input$add_class, {
      new_class <- trimws(input$new_class_name)
      if (new_class != "" && !(new_class %in% class_choices())) {
        updated_choices <- c(class_choices(), new_class)
        class_choices(updated_choices)
        shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "draw_class",
          choices = updated_choices,
          selected = new_class
        )
        shiny::updateTextInput(shiny::getDefaultReactiveDomain(), "new_class_name", value = "")
      }
    })

    # Dynamic Simulation UI
    output$sim_params_ui <- shiny::renderUI({
      shiny::req(input$sim_n_classes)
      n <- input$sim_n_classes
      if (is.na(n) || n < 1) {
        return(NULL)
      }

      lapply(1:n, function(i) {
        shiny::wellPanel(
          shiny::h5(paste("Class", i)),
          shiny::textInput(paste0("sim_mean_", i), "Mean (comma separated)", value = if (i == 1) "-1,0" else if (i == 2) "1,0" else "0,1"),
          shiny::textInput(paste0("sim_sd_", i), "Standard Deviations (comma separated)", value = "1,1"),
          shiny::numericInput(paste0("sim_cor_", i), "Correlation (rho)", value = 0, min = -0.99, max = 0.99, step = 0.1),
          shiny::numericInput(paste0("sim_n_", i), "Sample Size", value = 100, min = 10)
        )
      })
    })

    shiny::observeEvent(input$sim_do, {
      shiny::req(input$sim_do > 0)
      shiny::req(input$sim_n_classes)
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

        # Guard against NULL inputs if UI hasn't fully rendered
        if (is.null(m_str) || is.null(s_str)) return()

        m <- suppressWarnings(as.numeric(unlist(strsplit(m_str, ","))))
        sds <- suppressWarnings(as.numeric(unlist(strsplit(s_str, ","))))
        rho <- input[[paste0("sim_cor_", i)]]
        ns <- c(ns, input[[paste0("sim_n_", i)]])

        means[[i]] <- m

        # Construct covariance matrix using SDs and exchangeable correlation
        dim_m <- length(m)
        if (length(sds) != dim_m) {
          shiny::showNotification(sprintf("Class %d: Length of SDs must match length of means.", i), type = "error")
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
        shiny::showNotification("Error: All mean vectors must have the same number of dimensions.", type = "error")
        return()
      }
      if (any(is.na(unlist(means))) || any(is.na(unlist(covs)))) {
        shiny::showNotification("Error: Invalid numeric input in mean or SDs.", type = "error")
        return()
      }

      new_data <- tryCatch({
        simu_n(means = means, covs = covs, ns = ns)
      }, error = function(e) {
        shiny::showNotification(paste("Data generation failed:", e$message), type = "error")
        NULL
      })

      if (is.null(new_data)) return()

      current_data(new_data)
      injected_outliers(data.frame())
      new_classes <- unique(as.character(new_data$Sim))
      class_choices(new_classes)
      shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "draw_class", choices = new_classes, selected = new_classes[1])
      zoom_xlim(NULL)
      zoom_ylim(NULL)
    })

    shiny::observeEvent(input$sim_mixsim_do, {
      shiny::req(input$sim_mixsim_do > 0)
      shiny::req(input$sim_mixsim_k, input$sim_mixsim_p, input$sim_mixsim_omega, input$sim_mixsim_n)

      new_data <- tryCatch({
        simulate_mixsim(
          n = input$sim_mixsim_n,
          K = input$sim_mixsim_k,
          p = input$sim_mixsim_p,
          MaxOmega = input$sim_mixsim_omega
        )
      }, error = function(e) {
        shiny::showNotification(paste("MixSim data generation failed:", e$message), type = "error")
        NULL
      })

      if (is.null(new_data)) return()

      current_data(new_data)
      injected_outliers(data.frame())
      new_classes <- unique(as.character(new_data$Sim))
      class_choices(new_classes)
      shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "draw_class", choices = new_classes, selected = new_classes[1])
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

    shiny::observeEvent(input$clear, {
      current_data(data.frame(Sim = character(), X1 = numeric(), X2 = numeric()))
      injected_outliers(data.frame())

      new_classes <- c("Class 1", "Class 2", "Class 3")
      mode_states[["Draw Data"]]$classes <- new_classes
      class_choices(new_classes)

      shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "draw_class", choices = new_classes, selected = "Class 1")
      zoom_xlim(NULL)
      zoom_ylim(NULL)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$draw_total_classes, {
      shiny::req(input$draw_total_classes > 0)
      n <- as.integer(input$draw_total_classes)

      new_classes <- paste0("Class ", 1:n)

      mode_states[["Draw Data"]]$classes <- new_classes
      class_choices(new_classes)

      # Update UI to select first class
      shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "draw_class",
        choices = new_classes,
        selected = new_classes[1]
      )
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$plot_click, {
      if (input$data_mode == "Draw Data" && input$interaction_mode == "Draw Point") {
        if (ncol(current_data()) > 3) {
          return()
        } # Prevent drawing on high-dim data
        cd <- current_data()
        feat_cols <- setdiff(colnames(cd), "Sim")
        pt <- data.frame(Sim = input$draw_class)
        pt[[feat_cols[1]]] <- input$plot_click$x
        pt[[feat_cols[2]]] <- input$plot_click$y

        # Save history
        hist <- undo_history()
        hist[[length(hist) + 1]] <- cd
        undo_history(hist)

        current_data(rbind(cd, pt))
      }
    })

    shiny::observeEvent(input$plot_dblclick, {
      if (input$data_mode != "Draw Data" || input$interaction_mode == "Navigate") {
        zoom_xlim(NULL)
        zoom_ylim(NULL)
      }
    })

    shiny::observeEvent(input$plot_brush, {
      if (input$data_mode != "Draw Data" || input$interaction_mode == "Navigate") {
        b <- input$plot_brush
        cd <- current_data()

        # Allow zooming regardless of data size
        if (!is.null(b)) {
          zoom_xlim(c(b$xmin, b$xmax))
          zoom_ylim(c(b$ymin, b$ymax))
        }
      } else if (input$data_mode == "Draw Data" && input$interaction_mode == "Draw Cluster") {
        if (ncol(current_data()) > 3) {
          return()
        } # Prevent brushing on high-dim data
        b <- input$plot_brush
        n <- shiny::isolate(input$brush_size)

        # Ensure we don't try to generate NA points
        if (!is.null(b) && is.numeric(n) && n > 0) {
          cd <- current_data()
          feat_cols <- setdiff(colnames(cd), "Sim")
          pts <- data.frame(Sim = rep(shiny::isolate(input$draw_class), n))

          # Use mvrnorm to generate a Gaussian cluster in the brushed area
          mu_x <- (b$xmin + b$xmax) / 2
          mu_y <- (b$ymin + b$ymax) / 2
          sd_x <- (b$xmax - b$xmin) / 4
          sd_y <- (b$ymax - b$ymin) / 4

          # Fallback to runif if SD is very small (brush too small)
          if (sd_x > 0 && sd_y > 0) {
            if (!requireNamespace("MASS", quietly = TRUE)) {
              shiny::showNotification("The 'MASS' package is required for this feature. Please install it.", type = "error")
              return()
            }
            cluster <- MASS::mvrnorm(n, mu = c(mu_x, mu_y), Sigma = diag(c(sd_x^2, sd_y^2)))
            pts[[feat_cols[1]]] <- cluster[, 1]
            pts[[feat_cols[2]]] <- cluster[, 2]
          } else {
            pts[[feat_cols[1]]] <- stats::runif(n, b$xmin, b$xmax)
            pts[[feat_cols[2]]] <- stats::runif(n, b$ymin, b$ymax)
          }

          # Save history
          hist <- undo_history()
          hist[[length(hist) + 1]] <- cd
          undo_history(hist)

          current_data(rbind(cd, pts))
        }
      }
    })

    output$drawn_points_table <- DT::renderDataTable({
      dat <- current_data()
      if (is.null(dat) || nrow(dat) == 0) {
        return(DT::datatable(data.frame(Message = "No points drawn yet"), options = list(dom = 't'), rownames = FALSE))
      }
      # Show only the drawn points (which we'll assume is the current_data)
      # Format numeric columns
      num_cols <- sapply(dat, is.numeric)
      dat[num_cols] <- lapply(dat[num_cols], round, 2)

      DT::datatable(
        dat,
        options = list(pageLength = 5, lengthMenu = c(5, 10, 20), scrollX = TRUE),
        rownames = FALSE
      )
    })

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

      tryCatch({
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
      }, error = function(e) {
        shiny::showNotification(paste("Failed to import model:", e$message), type = "error")
      })
    })

    output$plot_grid <- shiny::renderUI({
      models <- input$selected_models
      if (is.null(models) || length(models) == 0) {
        return(shiny::p("Please select at least one model to compare."))
      }

      plot_outputs <- lapply(models, function(m) {
        col_width <- max(4, floor(12 / length(models)))
        safe_id <- gsub("[^a-zA-Z0-9_\\-]", "_", paste0("plot_", m))

        shiny::column(
          col_width,
          shiny::plotOutput(
            outputId = safe_id,
            click = "plot_click",
            dblclick = "plot_dblclick",
            brush = shiny::brushOpts(id = "plot_brush", resetOnNew = TRUE)
          )
        )
      })

      shiny::fluidRow(plot_outputs)
    })

    # Manual Tour Steering (High-Dimensional Mode)

    current_basis <- shiny::reactiveVal(NULL)
    current_path <- shiny::reactiveVal(NULL)
    current_projection <- shiny::reactiveVal(NULL)
    current_projection_info <- shiny::reactiveVal(NULL)

    output$tour_panel <- shiny::renderUI({
      dat <- current_data()
      if (ncol(dat) > 3) {
        shiny::wellPanel(
          shiny::h4("Manual Tour Steering"),
          shiny::p("High-dimensional data detected. Use the controls below to steer the 2D projection plane."),
          shiny::selectInput("tour_var", "Manipulation Variable", choices = setdiff(colnames(dat), "Sim")),
          shiny::sliderInput("tour_angle", "Rotation Angle", min = 0, max = 1, value = 0, step = 0.01)
        )
      }
    })

    # Initialize PCA basis when high-dim data loads
    shiny::observe({
      dat <- current_data()
      if (ncol(dat) > 3 && nrow(dat) >= 2 && is.null(current_basis())) {
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

      dat <- current_data()
      num_cols <- setdiff(colnames(dat), "Sim")
      var_idx <- which(num_cols == input$tour_var)

      # Following tourr::radial_tour logic to find the target basis:
      # Zero out the variable's contribution and orthonormalize
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

    comparison_state <- shiny::reactive({
      models <- input$selected_models
      dat <- combined_training_data()
      if (is.null(models) || length(models) == 0 || nrow(dat) < 2 || length(unique(dat$Sim)) < 2) {
        return(NULL)
      }

      state <- list(train_data = dat, models = list())

      for (m in models) {
        config <- app_methods[[m]]
        if (!is.null(config)) {
          plot_opts <- list()
          if (m == "PPtreeExtclass") plot_opts <- list(stop = input$stop)
          if (m == "Rpart") plot_opts <- list(control = rpart::rpart.control(cp = input$rpart_cp))
          if (m == "RandomForest") plot_opts <- list(ntree = input$rf_ntree)

          if (!is.null(config$fit_args_fn) && is.function(config$fit_args_fn)) {
            custom_opts <- config$fit_args_fn(input)
            plot_opts <- c(plot_opts, custom_opts)
          }

          fit_args <- c(config$args, plot_opts)

          tryCatch({
            cb_mod <- fit_model(dat, Sim ~ ., classifier = config$fn, fit_args = fit_args)
            pred_args <- predict_args[[m]](as.numeric(input$rule))
            preds <- predict_model(cb_mod, dat, predict_args = pred_args)

            state$models[[m]] <- list(
              model = cb_mod,
              predictions = preds,
              predict_args = pred_args
            )
          }, error = function(e) {
            warning(paste("Model", m, "failed to fit:", e$message))
            shiny::showNotification(paste("Model", m, "failed to fit:", e$message), type = "error", duration = 10)
          })
        }
      }

      state
    })

    # We must observe changes to selected models to register renderPlot for newly selected ones
    shiny::observe({
      models <- input$selected_models
      if (is.null(models)) {
        return()
      }

      lapply(models, function(m) {
        safe_id <- gsub("[^a-zA-Z0-9_\\-]", "_", paste0("plot_", m))
        output[[safe_id]] <- shiny::renderPlot({
          dat <- combined_training_data()

          title <- switch(m,
            "Rpart" = "Decision Tree (Rpart)",
            "PPtreeViz" = "PPtreeViz",
            "PPtreeExtclass" = "PPtreeExtclass",
            "PPtreeExt_split" = "PPtreeExt_split",
            "RandomForest" = "Random Forest",
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
                ggplot2::scale_color_discrete(drop = FALSE)
            }

            if (!is.null(zoom_xlim()) && !is.null(zoom_ylim())) {
              p <- p + ggplot2::coord_cartesian(xlim = zoom_xlim(), ylim = zoom_ylim())
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

          tryCatch(
            {
              create_boundary_plot(
                cb_mod = cb_mod,
                data = dat,
                title = title,
                class_levels = class_choices(),
                proj_matrix = current_projection(),
                proj_info = current_projection_info(),
                zoom_x = zoom_xlim(),
                zoom_y = zoom_ylim(),
                show_probs = isTRUE(input$show_probs),
                resolution = if (!is.null(input$grid_resolution)) input$grid_resolution else 100,
                predict_args = state$models[[m]]$predict_args
              )
            },
            error = function(e) {
              # Graceful fallback if model fails (e.g. data is degenerate)
              err_msg <- gsub("\033\\[[0-9;]*m", "", e$message)
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
             `Training Error` = NA,
             `Training Kappa` = NA,
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

        data.frame(
          Model = m,
          `Training Accuracy` = round(acc, 4),
          `Training Error` = round(err, 4),
          `Training Kappa` = round(kappa, 4),
          check.names = FALSE
        )
      })

      do.call(rbind, metrics_list)
    })

    output$metrics_table <- DT::renderDataTable({
      res <- metrics_df()
      if (is.null(res)) return(NULL)

      DT::datatable(
        res,
        options = list(
          dom = 't',
          paging = FALSE,
          searching = FALSE,
          ordering = TRUE
        ),
        rownames = FALSE,
        class = 'cell-border stripe hover'
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
          if(proj_active) "High-Dimensional Projection" else "Native 2D Feature Space"
        ),
        stringsAsFactors = FALSE
      )

      DT::datatable(
        df,
        rownames = FALSE,
        colnames = rep("", ncol(df)),
        options = list(
          dom = 't',
          bSort = FALSE,
          paging = FALSE,
          language = list(emptyTable = "Waiting for data...")
        ),
        selection = 'none'
      )
    })

    shiny::observeEvent(input$open_export_wizard, {
      shiny::showModal(shiny::modalDialog(
        title = "Export Wizard",
        shiny::checkboxGroupInput(
          "export_includes",
          "Include in Export:",
          choices = c("Plots", "Data", "Performance Metrics", "Fitted Models", "Session Info"),
          selected = c("Plots", "Data", "Performance Metrics", "Fitted Models", "Session Info")
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

    output$export_download <- shiny::downloadHandler(
      filename = function() {
        paste("classbound_export_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip", sep = "")
      },
      content = function(file) {
        # Close the modal upon clicking download
        shiny::removeModal()

        temp_dir <- tempdir()
        export_dir <- file.path(temp_dir, paste0("export_", as.integer(Sys.time())))
        dir.create(export_dir)

        dat <- combined_training_data()
        state <- comparison_state()
        metrics <- metrics_df()
        includes <- input$export_includes

        # 1. Export Data
        if ("Data" %in% includes && !is.null(dat) && nrow(dat) > 0) {
          export_data_csv(dat, file.path(export_dir, "data.csv"))
        }

        if (!is.null(state) && length(state$models) > 0) {
          # 2. Export Models
          if ("Fitted Models" %in% includes) {
            mod_list <- lapply(state$models, function(x) x$model)
            export_models_rds(mod_list, file.path(export_dir, "models.rds"))
          }

          # 3. Export Plots
          if ("Plots" %in% includes) {
            plots <- lapply(names(state$models), function(m) {
              title <- switch(m,
                "Rpart" = "Decision Tree (Rpart)",
                "PPtreeViz" = "PPtreeViz",
                "PPtreeExtclass" = "PPtreeExtclass",
                "PPtreeExt_split" = "PPtreeExt_split",
                "RandomForest" = "Random Forest",
                m
              )
              tryCatch({
                create_boundary_plot(
                  cb_mod = state$models[[m]]$model,
                  data = dat,
                  title = title,
                  class_levels = class_choices(),
                  proj_matrix = current_projection(),
                  proj_info = current_projection_info(),
                  zoom_x = zoom_xlim(),
                  zoom_y = zoom_ylim(),
                  show_probs = isTRUE(input$show_probs),
                  resolution = if (!is.null(input$grid_resolution)) input$grid_resolution else 100,
                  predict_args = state$models[[m]]$predict_args
                )
              }, error = function(e) NULL)
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

        # 4. Export Metrics
        if ("Performance Metrics" %in% includes && !is.null(metrics)) {
          export_metrics_csv(metrics, file.path(export_dir, "metrics.csv"))
        }

        # 5. Export Session Info
        if ("Session Info" %in% includes) {
          utils::capture.output(utils::sessionInfo(), file = file.path(export_dir, "session_info.txt"))
        }

        # Zip files
        owd <- setwd(export_dir)
        on.exit(setwd(owd))
        utils::zip(zipfile = file, files = list.files(export_dir), extras = "-q")
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
#' @param predict_args Arguments passed to predict
#' @return A ggplot object.
#' @keywords internal
create_boundary_plot <- function(cb_mod, data, title, class_levels = NULL, proj_matrix = NULL, proj_info = NULL, zoom_x = NULL, zoom_y = NULL, show_probs = FALSE, resolution = 100, predict_args = list()) {
  if (!is.null(proj_matrix)) {
    proj_list <- list(basis = proj_matrix, center = proj_info$center, scale = proj_info$scale)
    x_mat <- as.matrix(data[, rownames(proj_matrix)])
    if (!is.null(proj_info$center)) x_mat <- sweep(x_mat, 2, proj_info$center, "-")
    if (!is.null(proj_info$scale)) x_mat <- sweep(x_mat, 2, proj_info$scale, "/")
    z_mat <- x_mat %*% proj_matrix

    range_list <- list(
      PC1 = range(z_mat[, 1]) + c(-0.5, 0.5),
      PC2 = range(z_mat[, 2]) + c(-0.5, 0.5)
    )
    cb_bound <- boundary_compute(cb_mod, range = range_list, resolution = resolution, projection = proj_list, predict_args = predict_args)
    x_col_label <- "PC1"
    y_col_label <- "PC2"
  } else {
    feat_cols <- setdiff(colnames(data), "Sim")
    x_name <- feat_cols[1]
    y_name <- feat_cols[2]
    range_list <- list()
    range_list[[x_name]] <- range(data[[x_name]]) + c(-0.5, 0.5)
    range_list[[y_name]] <- range(data[[y_name]]) + c(-0.5, 0.5)

    cb_bound <- boundary_compute(cb_mod, range = range_list, resolution = resolution, predict_args = predict_args)
    x_col_label <- x_name
    y_col_label <- y_name
  }

  # Lock factor levels if provided to prevent color shifting
  if (!is.null(class_levels)) {
    data$Sim <- factor(data$Sim, levels = class_levels)
    cb_bound$boundary_data$prediction <- factor(cb_bound$boundary_data$prediction, levels = class_levels)
  }

  # Render plot
  p <- plot_boundary(cb_bound, obs_data = data, x_col = x_col_label, y_col = y_col_label, true_label = "Sim", show_gradient = show_probs) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(aspect.ratio = 1, legend.position = "none") +
    ggplot2::scale_color_discrete(drop = FALSE)

  if (!is.null(zoom_x) && !is.null(zoom_y)) {
    p <- suppressMessages(p + ggplot2::coord_cartesian(xlim = zoom_x, ylim = zoom_y, expand = FALSE))
  }

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
generate_extreme_outlier <- function(data, class_label, magnitude, target_col = "Sim", index = 1) {
  outlier <- data[1, , drop = FALSE]
  outlier[[target_col]] <- class_label

  feat_cols <- setdiff(colnames(data), target_col)
  num_cols <- names(which(sapply(data[feat_cols], is.numeric)))

  # Determine which corner to use (0 to 3) based on the index and class
  # Offset the index by a class hash so different classes start at different corners
  class_offset <- suppressWarnings(sum(utf8ToInt(as.character(class_label))))
  corner <- (index + class_offset) %% 4

  # Determine class-based magnitude offset to prevent different classes from overlapping at the same corner
  class_levels <- if (is.factor(data[[target_col]])) levels(data[[target_col]]) else unique(data[[target_col]])
  class_idx <- match(class_label, class_levels)
  if (is.na(class_idx)) class_idx <- length(class_levels) + 1

  # Apply a larger fraction of the magnitude (e.g. 3% per class index) to visibly separate classes
  class_sub_step <- (class_idx - 1) * 0.03 * magnitude

  for (col in feat_cols) {
    if (is.numeric(data[[col]])) {
      if (length(num_cols) >= 2 && col %in% num_cols[1:2]) {
        min_val <- min(data[[col]], na.rm = TRUE)
        max_val <- max(data[[col]], na.rm = TRUE)
        range_val <- diff(c(min_val, max_val))
        if (range_val == 0) range_val <- 1

        expansion <- (magnitude + class_sub_step) * range_val * 0.1

        # Corner mapping:
        # 0: max X1, max X2
        # 1: min X1, max X2
        # 2: min X1, min X2
        # 3: max X1, min X2
        use_max <- if (col == num_cols[1]) (corner %in% c(0, 3)) else (corner %in% c(0, 1))

        if (use_max) {
          outlier[[col]] <- max_val + expansion
        } else {
          outlier[[col]] <- min_val - expansion
        }
      } else {
        outlier[[col]] <- stats::median(data[[col]], na.rm = TRUE)
      }
    } else {
      freqs <- table(data[[col]])
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
  }

  return(outlier)
}
