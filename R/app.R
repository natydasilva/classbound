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
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' must be installed to use explorapp().", call. = FALSE)
  }
  if (!requireNamespace("MixSim", quietly = TRUE)) {
    stop("Package 'MixSim' must be installed to use explorapp().", call. = FALSE)
  }

  # UI to Package API Mapping
  app_methods <- list(
    "PPtreeViz"       = list(fn = PPtreeViz::PPTreeclass, args = list(PPmethod = "LDA")),
    "Rpart"           = list(fn = rpart::rpart, args = list()),
    "PPtreeExt_split" = list(fn = PPtreeExt::PPtreeExt_split, args = list(PPmethod = "LDA")),
    "PPtreeExtclass"  = list(fn = PPtreeExt::PPtreeExtclass, args = list(PPmethod = "LDA")),
    "RandomForest"    = list(fn = randomForest::randomForest, args = list())
  )

  predict_args <- list(
    "PPtreeViz"       = function(ru) list(Rule = ru),
    "Rpart"           = function(...) list(),
    "PPtreeExt_split" = function(ru) list(Rule = ru),
    "PPtreeExtclass"  = function(...) list(),
    "RandomForest"    = function(...) list()
  )

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
          shiny::radioButtons("data_mode", "Data Mode", choices = c("Draw Data", "Simulate Data", "Imported Data"), selected = if (!is.null(data)) "Imported Data" else "Draw Data"),
          shiny::helpText("Tip: Brush on any plot to zoom in, and double-click to reset the view."),
          shiny::conditionalPanel(
            condition = "input.data_mode == 'Draw Data'",
            shiny::hr(),
            shiny::radioButtons("interaction_mode", "Interaction Mode", choices = c("Draw Point", "Draw Cluster", "Navigate"), inline = TRUE),
            shiny::p("Click or brush on any plot to add data."),
            shiny::radioButtons("draw_class", "Active Class", choices = c("Class A", "Class B", "Class C"), inline = TRUE),
            shiny::div(
              style = "display: flex; gap: 10px; align-items: baseline; margin-bottom: 15px;",
              shiny::textInput("new_class_name", label = NULL, placeholder = "New class name...", width = "100%"),
              shiny::actionButton("add_class", "Add Class")
            ),
            shiny::numericInput("brush_size", "Brush Density (points/brush)", value = 20, min = 1),
            shiny::actionButton("clear", "Clear Canvas")
          ),
          shiny::conditionalPanel(
            condition = "input.data_mode == 'Simulate Data'",
            shiny::numericInput("sim_n_classes", "Number of Classes", value = 3, min = 2, max = 10),
            shiny::uiOutput("sim_params_ui"),
            shiny::actionButton("sim_do", "Generate Data")
          ),
          shiny::conditionalPanel(
            condition = "input.data_mode == 'Imported Data'",
            shiny::p("Using dataset provided via console.")
          )
        ),
        shiny::wellPanel(
          shiny::tags$details(
            shiny::tags$summary("Model Configuration", style = "display: list-item; font-size: 18px; font-weight: 500; cursor: pointer; margin-bottom: 10px;"),
            shiny::checkboxGroupInput(
              "selected_models",
              "Models to Compare",
              choices = names(app_methods),
              selected = c("Rpart", "PPtreeViz", "PPtreeExtclass")
            ),
            shiny::conditionalPanel(
              condition = "input.selected_models && (input.selected_models.indexOf('PPtreeViz') > -1 || input.selected_models.indexOf('PPtreeExtclass') > -1 || input.selected_models.indexOf('PPtreeExt_split') > -1)",
              shiny::hr(),
              shiny::selectInput("rule", "Projection Pursuit Rule", choices = 1:8, selected = 1)
            ),
            shiny::conditionalPanel(
              condition = "input.selected_models && input.selected_models.indexOf('PPtreeExtclass') > -1",
              shiny::numericInput("stop", "Stopping Rule", value = 4, min = 1)
            ),
            shiny::conditionalPanel(
              condition = "input.selected_models && input.selected_models.indexOf('Rpart') > -1",
              shiny::hr(),
              shiny::numericInput("rpart_cp", "Complexity Parameter", value = 0.01, min = 0, step = 0.01)
            ),
            shiny::conditionalPanel(
              condition = "input.selected_models && input.selected_models.indexOf('RandomForest') > -1",
              shiny::hr(),
              shiny::numericInput("rf_ntree", "Number of Trees", value = 500, min = 10, step = 50)
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
        shiny::uiOutput("tour_panel")
      ),
      shiny::mainPanel(
        shiny::uiOutput("data_stats_ui"),
        shiny::uiOutput("plot_grid")
      )
    )
  )


  create_boundary_plot <- function(data, test, meth, title, ru = 1, fit_opts = list(), class_levels = NULL, proj_matrix = NULL, proj_info = NULL, zoom_x = NULL, zoom_y = NULL) {
    config <- app_methods[[meth]]
    if (is.null(config)) stop("Unknown method: ", meth)
    fit_args <- c(config$args, fit_opts)
    predict_args <- predict_args[[meth]](ru)

    # Fit the model via the core pipeline
    cb_mod <- fit_model(data, Sim ~ ., classifier = config$fn, fit_args = fit_args)

    if (!is.null(proj_matrix)) {
      proj_list <- list(basis = proj_matrix, center = proj_info$center, scale = proj_info$scale)
      x_mat <- as.matrix(data[, rownames(proj_matrix)])
      if (!is.null(proj_info$scale)) x_mat <- sweep(x_mat, 2, proj_info$scale, "/")
      if (!is.null(proj_info$center)) x_mat <- sweep(x_mat, 2, proj_info$center, "-")
      z_mat <- x_mat %*% proj_matrix

      range_list <- list(
        PC1 = range(z_mat[, 1]) + c(-0.5, 0.5),
        PC2 = range(z_mat[, 2]) + c(-0.5, 0.5)
      )
      cb_bound <- boundary_compute(cb_mod, range = range_list, resolution = 100, projection = proj_list)
      x_col_label <- "PC1"
      y_col_label <- "PC2"
    } else {
      feat_cols <- setdiff(colnames(data), "Sim")
      x_name <- feat_cols[1]
      y_name <- feat_cols[2]
      range_list <- list()
      range_list[[x_name]] <- range(data[[x_name]]) + c(-0.5, 0.5)
      range_list[[y_name]] <- range(data[[y_name]]) + c(-0.5, 0.5)

      cb_bound <- boundary_compute(cb_mod, range = range_list, resolution = 100)
      x_col_label <- x_name
      y_col_label <- y_name
    }

    # Calculate test error for the UI title
    preds <- predict_model(cb_mod, test, predict_args = predict_args)
    err <- round(mean(preds$class != test$Sim) * 100, 3)

    # Lock factor levels if provided to prevent color shifting
    if (!is.null(class_levels)) {
      data$Sim <- factor(data$Sim, levels = class_levels)
      cb_bound$boundary_data$prediction <- factor(cb_bound$boundary_data$prediction, levels = class_levels)
    }

    # Render plot
    p <- plot_boundary(cb_bound, obs_data = data, x_col = x_col_label, y_col = y_col_label, true_label = "Sim") +
      ggplot2::ggtitle(paste0(title, " (Error ", err, "%)")) +
      ggplot2::theme(aspect.ratio = 1, legend.position = "none") +
      ggplot2::scale_color_discrete(drop = FALSE)

    if (!is.null(zoom_x) && !is.null(zoom_y)) {
      p <- p + ggplot2::coord_cartesian(xlim = zoom_x, ylim = zoom_y)
    }

    p
  }

  # Server
  server <- function(input, output) {
    # Initialize current_data based on passed 'data'
    initial_data <- data.frame(Sim = character(), X1 = numeric(), X2 = numeric())
    initial_classes <- c("Class A", "Class B", "Class C")

    if (!is.null(data)) {
      if (is.null(target_col) || !(target_col %in% colnames(data))) {
        stop("Please provide a valid 'target_col' that exists in 'data'.")
      }
      initial_data <- stats::na.omit(data)
      # Standardize target column to "Sim" for internal app logic
      colnames(initial_data)[colnames(initial_data) == target_col] <- "Sim"

      # Keep only numeric features and the target
      numeric_cols <- sapply(initial_data, is.numeric)
      numeric_cols[which(colnames(initial_data) == "Sim")] <- TRUE
      initial_data <- initial_data[, numeric_cols, drop = FALSE]

      initial_classes <- as.character(unique(initial_data$Sim))
      if (length(initial_classes) == 0) {
        initial_classes <- c("Class A", "Class B", "Class C")
      }
    }

    current_data <- shiny::reactiveVal(initial_data)
    imported_data <- shiny::reactiveVal(if (!is.null(data)) initial_data else NULL)
    class_choices <- shiny::reactiveVal(initial_classes)
    zoom_xlim <- shiny::reactiveVal(NULL)
    zoom_ylim <- shiny::reactiveVal(NULL)

    output$data_stats_ui <- shiny::renderUI({
      dat <- current_data()
      if (nrow(dat) == 0) {
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
      ncol(current_data()) > 3
    })
    shiny::outputOptions(output, "is_high_dim", suspendWhenHidden = FALSE)

    shiny::observeEvent(input$data_mode, {
      zoom_xlim(NULL)
      zoom_ylim(NULL)

      if (input$data_mode == "Imported Data") {
        if (!is.null(imported_data())) {
          current_data(imported_data())
          new_classes <- unique(as.character(imported_data()$Sim))
          class_choices(new_classes)
          shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "draw_class", choices = new_classes, selected = new_classes[1], inline = TRUE)
        }
      } else {
        # Switch to Draw/Simulate Data: completely reset the canvas and classes
        current_data(data.frame(Sim = character(), X1 = numeric(), X2 = numeric()))
        default_classes <- c("Class A", "Class B", "Class C")
        class_choices(default_classes)
        shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "draw_class", choices = default_classes, selected = default_classes[1], inline = TRUE)
      }
    })

    shiny::observeEvent(input$add_class, {
      new_class <- trimws(input$new_class_name)
      if (new_class != "" && !(new_class %in% class_choices())) {
        updated_choices <- c(class_choices(), new_class)
        class_choices(updated_choices)
        shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "draw_class",
          choices = updated_choices,
          selected = new_class,
          inline = TRUE
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
          shiny::textInput(paste0("sim_cov_", i), "Covariance Diag (comma separated)", value = "1,1"),
          shiny::numericInput(paste0("sim_n_", i), "Sample Size", value = 100, min = 10)
        )
      })
    })

    shiny::observeEvent(input$sim_do, {
      shiny::req(input$sim_n_classes)
      n <- input$sim_n_classes
      if (is.na(n) || n < 1) {
        return()
      }

      means <- list()
      covs <- list()
      ns <- c()

      for (i in 1:n) {
        m <- as.numeric(unlist(strsplit(input[[paste0("sim_mean_", i)]], ",")))
        cv <- as.numeric(unlist(strsplit(input[[paste0("sim_cov_", i)]], ",")))
        ns <- c(ns, input[[paste0("sim_n_", i)]])

        means[[i]] <- m
        covs[[i]] <- diag(cv, length(m))
      }

      lengths <- sapply(means, length)
      if (length(unique(lengths)) > 1) {
        shiny::showNotification("Error: All mean vectors must have the same number of dimensions.", type = "error")
        return()
      }
      cv_lengths <- sapply(covs, nrow)
      if (any(cv_lengths != lengths)) {
        shiny::showNotification("Error: Covariance diagonals must have the same length as their mean vectors.", type = "error")
        return()
      }
      if (any(is.na(unlist(means))) || any(is.na(unlist(covs)))) {
        shiny::showNotification("Error: Invalid numeric input in mean or covariance.", type = "error")
        return()
      }

      new_data <- simu_n(means = means, covs = covs, ns = ns)
      current_data(new_data)
      new_classes <- unique(as.character(new_data$Sim))
      class_choices(new_classes)
      shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "draw_class", choices = new_classes, selected = new_classes[1], inline = TRUE)
      zoom_xlim(NULL)
      zoom_ylim(NULL)
    })

    shiny::observeEvent(input$clear, {
      current_data(data.frame(Sim = character(), X1 = numeric(), X2 = numeric()))
      class_choices(c("Class A", "Class B", "Class C"))
      shiny::updateRadioButtons(shiny::getDefaultReactiveDomain(), "draw_class", choices = c("Class A", "Class B", "Class C"), selected = "Class A", inline = TRUE)
      zoom_xlim(NULL)
      zoom_ylim(NULL)
    })

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
          pts <- data.frame(Sim = shiny::isolate(input$draw_class))
          pts[[feat_cols[1]]] <- stats::runif(n, b$xmin, b$xmax)
          pts[[feat_cols[2]]] <- stats::runif(n, b$ymin, b$ymax)
          current_data(rbind(cd, pts))
        }
      }
    })

    output$plot_grid <- shiny::renderUI({
      models <- input$selected_models
      if (is.null(models) || length(models) == 0) {
        return(shiny::p("Please select at least one model to compare."))
      }

      plot_outputs <- lapply(models, function(m) {
        col_width <- max(4, floor(12 / length(models)))

        shiny::column(
          col_width,
          shiny::plotOutput(
            outputId = paste0("plot_", m),
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

        pca <- stats::prcomp(num_dat, scale. = scale_flag)
        basis <- pca$rotation[, 1:2]
        basis <- tourr::orthonormalise(basis)
        current_basis(basis)
        current_projection(basis)
        current_projection_info(list(center = pca$center, scale = pca$scale))
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

    #

    # We must observe changes to selected models to register renderPlot for newly selected ones
    shiny::observe({
      models <- input$selected_models
      if (is.null(models)) {
        return()
      }

      lapply(models, function(m) {
        output[[paste0("plot_", m)]] <- shiny::renderPlot({
          dat <- current_data()

          title <- switch(m,
            "Rpart" = "Decision Tree (Rpart)",
            "PPtreeViz" = "PPtreeViz",
            "PPtreeExtclass" = "PPtreeExtclass",
            "PPtreeExt_split" = "PPtreeExt_split",
            "RandomForest" = "Random Forest",
            m
          )

          if (nrow(dat) < 2 || length(unique(dat$Sim)) < 2) {
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
              dat$Sim <- factor(dat$Sim, levels = class_choices())
              p <- p + ggplot2::geom_point(ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]], color = .data$Sim), size = 3) +
                ggplot2::scale_color_discrete(drop = FALSE)
            }
            return(p)
          }

          # Safe fit and plot
          plot_opts <- list()
          if (m == "PPtreeExtclass") plot_opts <- list(stop = input$stop)
          if (m == "Rpart") plot_opts <- list(control = rpart::rpart.control(cp = input$rpart_cp))
          if (m == "RandomForest") plot_opts <- list(ntree = input$rf_ntree)

          c_mod <- app_methods[[m]]
          if (!is.null(c_mod$fit_args_fn) && is.function(c_mod$fit_args_fn)) {
            custom_opts <- c_mod$fit_args_fn(input)
            plot_opts <- c(plot_opts, custom_opts)
          }

          tryCatch(
            {
              create_boundary_plot(
                data = dat,
                test = dat,
                meth = m,
                title = m,
                ru = as.numeric(input$rule),
                fit_opts = plot_opts,
                class_levels = class_choices(),
                proj_matrix = current_projection(),
                proj_info = current_projection_info(),
                zoom_x = zoom_xlim(),
                zoom_y = zoom_ylim()
              )
            },
            error = function(e) {
              # Graceful fallback if model fails (e.g. data is degenerate)
              ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0, y = 0, label = paste("Error:", e$message)) +
                ggplot2::labs(title = title) +
                ggplot2::theme_minimal() +
                ggplot2::theme(aspect.ratio = 1)
            }
          )
        })
      })
    })
  }

  shiny::shinyApp(ui = ui, server = server)
}
