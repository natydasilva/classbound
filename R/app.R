#' Shiny app to compare PPtree, PPtreeExt and rpart boundaries in 2D with different simulation scenarios
#'
#' @usage explorapp()
#' @return No return value, called for side effects. Shinyapp is launched.
#' @export
#' @examples
#' if (interactive()) {
#'   explorapp()
#' }
#'
explorapp <- function() {
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

  X1 <- NULL
  X2 <- NULL
  ppred <- NULL
  Sim <- NULL
  pred <- NULL
  predict <- NULL


  # UI ----------------------------------------------------------------------
  ui <- shiny::fluidPage(shiny::mainPanel(
    shiny::tabsetPanel(
      shiny::tabPanel(
        "Basic-Sim",
        shiny::fluidRow(
          shiny::column(
            3,
            shiny::selectInput(
              inputId = "rule",
              label = "Rule",
              choices = 1:8,
              selected = 1
            )
          ),
          shiny::column(
            3,
            shiny::selectInput(
              inputId = "modi",
              label = "Modification",
              choices = c(
                "Subsetting clases" = "1",
                "Multiple splits" = "3"
              ),
              selected = 3
            )
          )
        ),
        shiny::fluidRow(
          shiny::column(
            4,
            shiny::textInput(
              inputId = "mean",
              label = "Group means ",
              value =
                "-1, 0.6, 0, -0.6, 2,-1"
            )
          ),
          shiny::column(
            4,
            shiny::textInput(
              inputId = "cor",
              label = "Correlations",
              value = "0.95, 0.5, 0.95"
            )
          ),
          shiny::column(
            4,
            shiny::textInput(
              inputId = "sample",
              label = "Group sample",
              value = "100, 100, 100"
            )
          )
        ),
        shiny::fluidRow(shiny::actionButton("do", label = "OK")),
        shiny::fluidRow(shiny::plotOutput("distPlot"))
      ),
      shiny::tabPanel(
        "SIM-Outliers",
        shiny::fluidRow(
          shiny::column(
            4,
            shiny::selectInput(
              inputId = "rule2",
              label = "Rule",
              choices = 1:8,
              selected = 1
            )
          ),
          shiny::column(
            3,
            shiny::selectInput(
              inputId = "modi2",
              label = "Modification",
              choices = c(
                "Subsetting clases" = "1",
                "Multiple splits" = "3"
              ),
              selected = 3
            )
          )
        ),
        shiny::fluidRow(
          shiny::column(
            4,
            shiny::textInput(
              inputId = "mean2",
              label = "Group means ",
              value =
                "-1, 0.6, 0, -0.6, 2,-1"
            )
          ),
          shiny::column(
            4,
            shiny::textInput(
              inputId = "cor2",
              label = "Correlations",
              value = "0.95, 0.95, 0.95"
            )
          ),
          shiny::column(
            4,
            shiny::textInput(
              inputId = "sample2",
              label = "Group sample",
              value = "100, 100, 100"
            )
          )
        ),
        shiny::fluidRow(shiny::column(
          4,
          shiny::selectInput(
            inputId = "group",
            label = "Add outliers to class",
            choices = 1:3,
            selected = 2
          )
        )),
        shiny::fluidRow(
          shiny::column(
            4,
            shiny::textInput(
              inputId = "meanout",
              label = "Out. X1, X2 means ",
              value =
                "-3, 3"
            )
          ),
          shiny::column(
            4,
            shiny::textInput(
              inputId = "sdout",
              label = "Out. X1, X2 sd ",
              value = ".2,.2"
            )
          ),
          shiny::column(
            4,
            shiny::textInput(
              inputId = "sampleout",
              label = "Out. sample size",
              value = "50"
            )
          ),
          shiny::fluidRow(shiny::actionButton("do2", label = "OK"))
        ),
        shiny::fluidRow(shiny::plotOutput("distPlot2"))
      ),
      ##
      shiny::tabPanel(
        "MixSim",
        shiny::fluidRow(
          shiny::column(
            4,
            shiny::selectInput(
              inputId = "rule3",
              label = "Rule",
              choices = 1:8,
              selected = 1
            )
          ),
          shiny::column(
            4,
            shiny::selectInput(
              inputId = "modi3",
              label = "Modification",
              choices = c(
                "Subsetting clases" = "1",
                "Multiple splits" = "3"
              ),
              selected = 3
            )
          )
        ),
        shiny::fluidRow(
          shiny::column(4, shiny::numericInput(
            "size",
            label = "Sample size", value = 500
          )),
          shiny::column(
            4,
            shiny::numericInput("BarOmega", label = "BarOmega desired average overlap", value = 0.05)
          )
        ),
        shiny::fluidRow(
          shiny::column(
            4,
            shiny::numericInput("MaxOmega", label = "MaxOmega desired maximum overlap", value = 0.15)
          ),
          shiny::column(
            4,
            shiny::numericInput("K", label = "K number of components", value = 4)
          )
        ),
        # shiny::numericInput("p", label = "number of dimensions", value = 5),
        shiny::fluidRow(shiny::actionButton("simmaitra", "OK")),
        shiny::fluidRow(shiny::plotOutput("plotsmaitra"))
        ##
      )
    )
  ))



# UI to Package API Mapping -----------------------------------------------
APP_METHODS <- list(
  "Original" = list(fn = PPtreeViz::PPTreeclass,     args = list(PPmethod = "LDA")),
  "Rpart"    = list(fn = rpart::rpart,               args = list()),
  "Modified" = list(fn = PPtreeExt::PPtreeExt_split, args = list(PPmethod = "LDA")),
  "MOD"      = list(fn = PPtreeExt::PPtreeExtclass,  args = list(PPmethod = "LDA"))
)

PREDICT_ARGS <- list(
  "Original" = function(ru) list(Rule = ru),
  "Rpart"    = function(...) list(type = "class"),
  "Modified" = function(ru) list(Rule = ru),
  "MOD"      = function(...) list()
)

create_boundary_plot <- function(data, test, meth, title, ru = 1, fit_opts = list()) {
  config <- APP_METHODS[[meth]]
  fit_args <- c(config$args, fit_opts)
  predict_args <- PREDICT_ARGS[[meth]](ru)

  # Fit the model via the core pipeline
  cb_mod <- fit_model(data, Sim ~ ., classifier = config$fn, fit_args = fit_args)

  # Compute decision boundary (standard 0.5 padding on standardized data)
  range_list <- list(
    X1 = range(data$X1) + c(-0.5, 0.5),
    X2 = range(data$X2) + c(-0.5, 0.5)
  )
  cb_bound <- boundary_compute(cb_mod, range = range_list, resolution = 100)

  # Calculate test error for the UI title
  preds <- predict_model(cb_mod, test, predict_args = predict_args)
  err <- round(mean(preds$class != test$Sim) * 100, 3)

  # Render plot
  plot_boundary(cb_bound, obs_data = data, x_col = "X1", y_col = "X2", true_label = "Sim") +
    ggplot2::ggtitle(paste0(title, " (test error ", err, "%)")) +
    ggplot2::theme(aspect.ratio = 1, legend.position = "none") +
    ggplot2::scale_fill_brewer(name = "Class", type = "qual", palette = "Dark2")
}

  # Server ------------------------------------------------------------------
  server <- function(input, output) {
    output$distPlot <- shiny::renderPlot({
      if (input$do) {
        x1 <- shiny::isolate(as.numeric(unlist(strsplit(input$mean, ","))))
        x2 <- shiny::isolate(as.numeric(unlist(strsplit(input$cor, ","))))
        x3 <- shiny::isolate(as.numeric(unlist(strsplit(input$sample, ","))))
        x4 <- shiny::isolate(as.numeric(input$stop))
        dat.pl2 <-
          shiny::isolate(simu3(
            x1[1], x1[2], x1[3], x1[4], x1[5], x1[6],
            x2[1], x2[2], x2[3], x3[1], x3[2], x3[3]
          ))
        dat.test <-
          shiny::isolate(simu3(
            x1[1], x1[2], x1[3], x1[4], x1[5], x1[6],
            x2[1], x2[2], x2[3], round(x3[1] * 0.25), round(x3[2] * 0.25), round(x3[3] * 0.25)
          ))

        if (input$modi == 1) {
          modpl <-
            create_boundary_plot(
              ru = 1,
              data = dat.pl2,
              test = dat.test,
              meth = "Modified",
              title = "PPtreeExt: Subsetting clases",
              fit_opts = list(entro = FALSE)
            )
        }
        # if (input$modi == 2) {
        #   #entropy mp groups
        #   modpl <-
        #     ppbound(
        #       ru =  1, #as.numeric(input$rule),
        #       data = dat.pl2,
        #       test = dat.test,
        #       meth = "Modified" ,
        #       entro = TRUE,
        #       title = "Modified 2 "
        #     )
        # }
        if (input$modi == 3) {
          modpl <-
            create_boundary_plot(
              data = dat.pl2,
              test = dat.test,
              meth = "MOD",
              title = "PPtreeExt: Multiple splits",
              fit_opts = list(strule = x4, tot = sum(x3))
            )
        }

        gridExtra::grid.arrange(
          create_boundary_plot(
            ru = as.numeric(input$rule),
            data = dat.pl2,
            test = dat.test,
            meth = "Rpart",
            title = "Rpart"
          ),
          create_boundary_plot(
            ru = as.numeric(input$rule),
            data = dat.pl2,
            test = dat.test,
            meth = "Original",
            title = "PPtree"
          ),
          modpl,
          ncol = 3
        )
      }
    })


    output$distPlot2 <- shiny::renderPlot({
      if (input$do2) {
        x1 <-
          shiny::isolate(as.numeric(unlist(strsplit(input$mean2, ","))))
        x2 <-
          shiny::isolate(as.numeric(unlist(strsplit(input$cor2, ","))))
        x3 <-
          shiny::isolate(as.numeric(unlist(strsplit(input$sample2, ","))))
        x4 <-
          shiny::isolate(as.numeric(unlist(strsplit(input$meanout, ","))))
        x5 <-
          shiny::isolate(as.numeric(unlist(strsplit(input$sdout, ","))))
        x6 <-
          shiny::isolate(as.numeric(unlist(strsplit(input$sampleout, ","))))
        x7 <-
          shiny::isolate(input$stop2)
        dat.pl2 <- simu3(
          x1[1], x1[2], x1[3], x1[4], x1[5], x1[6],
          x2[1], x2[2], x2[3], x3[1], x3[2], x3[3]
        )
        dat.test <-
          shiny::isolate(simu3(
            x1[1], x1[2], x1[3], x1[4], x1[5], x1[6],
            x2[1], x2[2], x2[3], round(x3[1] * 0.25), round(x3[2] * 0.25), round(x3[3] * 0.25)
          ))

        aux <-
          data.frame(
            Sim = rep(paste("sim", as.numeric(input$group), sep = ""), x6),
            X1 = stats::rnorm(
              n = x6,
              mean = x4[1],
              sd = x5[1]
            ),
            X2 = stats::rnorm(
              n = x6,
              mean = x4[2],
              sd = x5[2]
            )
          )

        aux2 <-
          data.frame(
            Sim = rep(paste("sim", as.numeric(input$group), sep = ""), round(x6 * 0.25)),
            X1 = stats::rnorm(
              n = round(x6 * 0.25),
              mean = x4[1],
              sd = x5[1]
            ),
            X2 = stats::rnorm(
              n = round(x6 * 0.25),
              mean = x4[2],
              sd = x5[2]
            )
          )
        dat.pl2 <- rbind(dat.pl2, aux)
        dat.test <- rbind(dat.test, aux2)

        if (input$modi2 == 1) {
          modpl <-
            create_boundary_plot(
              ru = as.numeric(input$rule),
              data = dat.pl2,
              test = dat.test,
              meth = "Modified",
              title = "PPtreeExt: Subsetting clases",
              fit_opts = list(entro = FALSE)
            )
        }
        if (input$modi2 == 2) {
          modpl <-
            create_boundary_plot(
              ru = as.numeric(input$rule),
              data = dat.pl2,
              test = dat.test,
              meth = "Modified",
              title = "Modified 2",
              fit_opts = list(entro = TRUE)
            )
        }
        if (input$modi2 == 3) {
          modpl <-
            create_boundary_plot(
              data = dat.pl2,
              test = dat.test,
              meth = "MOD",
              title = "PPtreeExt: Multiple splits",
              fit_opts = list(strule = x7, tot = sum(x3 + x6))
            )
        }

        gridExtra::grid.arrange(
          create_boundary_plot(
            ru = as.numeric(input$rule2),
            data = dat.pl2,
            test = dat.test,
            meth = "Rpart",
            title = "Rpart"
          ),
          create_boundary_plot(
            ru = as.numeric(input$rule2),
            data = dat.pl2,
            test = dat.test,
            meth = "Original",
            title = "PPtree"
          ),
          modpl,
          ncol = 3
        )
      }
    })


    output$plotsmaitra <- shiny::renderPlot({
      if (input$simmaitra) {
        x1 <- shiny::isolate(as.numeric(input$stop3))
        # Q <- MixSim(BarOmega = 0.01, K = 4, p = 2)

        repeat {
          Q <-
            MixSim::MixSim(
              BarOmega = shiny::isolate(as.numeric(input$BarOmega)),
              MaxOmega = shiny::isolate(as.numeric(input$MaxOmega)),
              K = shiny::isolate(as.numeric(input$K)),
              p = 2,
              # sph = FALSE,
              # ecc = 1,
              # PiLow = 1.0,
              # int = c(0.0, 1.0),
              # resN = 100,
              # eps = 1e-06,
              # lim = 1e06
            )
          if (Q$fail == 0) {
            break
          }
        }
        A <-
          MixSim::simdataset(
            n = shiny::isolate(as.numeric(input$size)),
            Pi = Q$Pi,
            Mu = Q$Mu,
            S = Q$S
          )

        Atest <-
          MixSim::simdataset(
            n = shiny::isolate(as.numeric(round(input$size * 0.25))),
            Pi = Q$Pi,
            Mu = Q$Mu,
            S = Q$S
          )
        dat.pl2 <-
          data.frame(
            Sim = paste("sim", A[[2]], sep = ""),
            X1 = scale(A[[1]][, 1]),
            X2 = scale(A[[1]][, 2])
          )
        dat.test <-
          data.frame(
            Sim = paste("sim", Atest[[2]], sep = ""),
            X1 = scale(Atest[[1]][, 1]),
            X2 = scale(Atest[[1]][, 2])
          )

        if (input$modi3 == 1) {
          modpl <-
            create_boundary_plot(
              ru = as.numeric(input$rule3),
              data = dat.pl2,
              test = dat.test,
              meth = "Modified",
              title = "PPtreeExt: Subsetting clases",
              fit_opts = list(entro = FALSE)
            )
        }
        if (input$modi3 == 2) {
          modpl <-
            create_boundary_plot(
              ru = as.numeric(input$rule3),
              data = dat.pl2,
              test = dat.test,
              meth = "Modified",
              title = "Modified 2",
              fit_opts = list(entro = TRUE)
            )
        }
        if (input$modi3 == 3) {
          modpl <-
            create_boundary_plot(
              data = dat.pl2,
              test = dat.test,
              meth = "MOD",
              title = "PPtreeExt: Multiple splits",
              fit_opts = list(strule = x1, tot = input$size)
            )
        }

        gridExtra::grid.arrange(
          create_boundary_plot(
            ru = as.numeric(input$rule3),
            data = dat.pl2,
            test = dat.test,
            meth = "Rpart",
            title = "Rpart"
          ),
          create_boundary_plot(
            ru = as.numeric(input$rule3),
            data = dat.pl2,
            test = dat.test,
            meth = "Original",
            title = "PPtree"
          ),
          modpl,
          ncol = 3
        )
      }
    })
  }

  shiny::shinyApp(ui = ui, server = server)
}
