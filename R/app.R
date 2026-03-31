#' Shiny app to compare PPtree, PPtreeExt and rpart boundaries in 2D with different simulation scenarios
 #' 
 #' @usage explorapp(ui, server) 
 #' @param ui user interface
 #' @param server server function
 #' @return No return value, called for side effects. Shinyapp is launched. 
 #' @export
 #' @examples
 #' if(interactive()){
 #' explorapp(ui,server)
 #' }
 #' 
 explorapp <-function(ui,server){
   
   X1 <- NULL
   X2 <- NULL
   ppred <- NULL
   Sim <- NULL
   pred <- NULL
   predict <- NULL
   
   # Function to simulate multivariate (bivariate) normal distributions
   #  simu3(mux1, mux2, muy1, muy2, muz1, muz2, cor1, cor2, cor3, n1 = 100, n2 = 100, n3 = 100) 
   #  mux1 mean of X1 for class 1
   #  mux2 mean of X2 for class 1
   #  muy1 mean of X1 for class 2
   #  muy2 mean of X2 for class 2
   #  muz1 mean of X1 for class 3
   #  muz2 mean of X2 for class 3
   #  cor1 correlation for class 1
   #  cor2 correlation for class 2
   #  cor3 correlation for class 3
   #  n1 number of samples for class 1
   #  n2 number of samples for class 2
   #  n3 number of samples for class 3
   #  data.frame with dimension (n1+n2+n3)x(3)
   
   simu3 <- 
     function(mux1,
              mux2,
              muy1,
              muy2,
              muz1,
              muz2,
              cor1,
              cor2,
              cor3,
              n1 = 100,
              n2 = 100,
              n3 = 100
     ) {
       bivn <- MASS::mvrnorm(n1, mu = c(mux1, mux2), Sigma = matrix(c(1, cor1, cor1, 1), 2))
       bivn2 <- MASS::mvrnorm(n2, mu = c(muy1, muy2), Sigma = matrix(c(1, cor2, cor2, 1), 2))
       bivn3 <- MASS::mvrnorm(n3, mu = c(muz1, muz2), Sigma = matrix(c(1, cor3, cor3, 1), 2))
       
       d1 <- data.frame(Sim = "sim1", bivn)
       d2 <- data.frame(Sim = "sim2", bivn2)
       d3 <- data.frame(Sim = "sim3", bivn3)
       return(rbind(d1, d2, d3))
     }
   
   
   # Description: Function to generate grid values for plotting decision boundaries for modification 1
   #  modifying the choice of split points-through class subsetting
   # ppbound(ru, data , test, meth, entro , title, simM = FALSE)
   # ru = split rule = {1, 2, 3, 4, 5, 6, 7, 8}
   # data data frame with the simulated dataset 
   # test data frame simulated test data 
   # entro logical; if TRUE use entropy in the modified PPtreeExt
   # meth character; method to use "Original" for PPtree, "Rpart" for rpart, "Modified" for PPtreeExt with subsetting classes
   # title character; title for the plot
   # simM logical; if TRUE use shapes and colors for classes, if FALSE use only colors
   
   ppbound <- function(ru, data , test, meth, entro , title, simM = FALSE) {
     grilla <-
       base::expand.grid(
         X1 = seq((min(data$X1) + sign(min(data$X1)) * .5), (max(data$X1) + sign(max(data$X1)) * .5), length.out = 100),
         X2 = seq((min(data$X2) + sign(min(data$X2)) * .5), (max(data$X2) + sign(max(data$X2)) * .5), length.out = 100)
       )
     
     data$Sim <- as.factor(data$Sim)
     
     if (meth == "Original") {
       pptree <- PPtreeViz::PPTreeclass(Sim ~ ., data = data, "LDA")
       ppred.sim <-
         PPtreeViz::PPclassify(pptree, test.data = grilla, Rule = ru)
       grilla$pred <- ppred.sim[[2]]
       err <-
         round(
           PPtreeViz::PPclassify(
             pptree,
             test.data = test[, -1],
             true.class = test[, 1],
             Rule = ru
           )[[1]] / nrow(test[, -1]),
           3
         ) * 100
     }
     if (meth == "Rpart") {
       rpart.mod <- rpart::rpart(Sim ~ ., data = data)
       grilla$pred <-
         predict(rpart.mod, newdata = grilla, type = "class")
       err <-
         round(1 - sum(diag(table(
           predict(rpart.mod, newdata = test[, -1], type = "class") , test[, 1]
         ))) / nrow(test[, -1]), 3) * 100
       
     }
     
     if (meth == "Modified") {
       #Esto esta mal?? estoy usando PPtreeViz
       pptree <- PPtreeExt::PPtreeExt_split(Sim ~ ., data = data, "LDA", entro = entro)
       #ppred.sim <- PPtreeViz::PPclassify(pptree, test.data = grilla, Rule = ru)
       ppred.sim <- predict(pptree, newdata  = grilla, Rule = ru)
       #grilla$pred <- paste("sim", ppred.sim[[2]], sep = "")
       grilla$pred <- ppred.sim
       aux_test <- predict(
         pptree,
         newdata= test[, -1],
         true.class = test[, 1],
         Rule = ru
       ) 
       err <- round(1- sum(aux_test == test[, 1])/length(test[, 1]), 3) * 100
       
       
     }
     
     if (meth == "RandomForest") {
       if (!requireNamespace("randomForest", quietly = TRUE)) {
         stop(
           "Package 'randomForest' is required for method 'RandomForest'. ",
           "Please install it with: install.packages('randomForest')"
         )
       }
       ntree_val <- getOption("ppbound.randomForest.ntree", 500)
       rf.mod <- randomForest::randomForest(Sim ~ ., data = data, ntree = ntree_val)
       grilla$pred <- predict(rf.mod, newdata = grilla, type = "response")
       err <- round(1 - sum(diag(table(
         predict(rf.mod, newdata = test[, -1], type = "response"), test[, 1]
       ))) / nrow(test[, -1]), 3) * 100
     }
     
     #ruleid <- pptree$splitCutoff.node[,ru]
     if (simM) {
       pl.pp  <-
         ggplot2::ggplot(data = grilla) + ggplot2::geom_point(ggplot2::aes(
           x = X1,
           y = X2,
           color = as.factor(pred)
         ), alpha = .20) +
         ggplot2::scale_colour_brewer(name = "Class",
                                      type = "qual",
                                      palette = "Dark2") + ggplot2::theme_bw() +
         ggplot2::geom_point(
           data = data,
           ggplot2::aes(
             x = X1 ,
             y = X2,
             group = Sim,
             color = Sim
           ),
           size = I(3)
         ) +
         ggplot2::theme(legend.position = "none", aspect.ratio = 1) +
         ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::scale_x_continuous(expand = c(0, 0)) +
         ggplot2::labs(
           x = " ",
           y = "",
           title = paste(title, "\n(test error: ", err, "%)", sep = '')
         )
     } else{
       pl.pp  <-
         ggplot2::ggplot(data = grilla) + ggplot2::geom_point(ggplot2::aes(
           x = X1,
           y = X2,
           color = as.factor(pred),
           shape = as.factor(pred)
         ),
         alpha = .20) +
         ggplot2::scale_colour_brewer(name = "Class",
                                      type = "qual",
                                      palette = "Dark2") + ggplot2::theme_bw() +
         ggplot2::scale_shape_discrete(name = 'Class') + ggplot2::geom_point(
           data = data,
           ggplot2::aes(
             x = X1 ,
             y = X2,
             group = Sim,
             shape = Sim,
             color = Sim
           ),
           size = I(3)
         ) +
         ggplot2::theme(legend.position = "none", aspect.ratio = 1) +
         ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::scale_x_continuous(expand = c(0, 0)) +
         ggplot2::labs(
           x = " ",
           y = "",
           title = paste(title, "\n(test error: ", err, "%)", sep = '')
         )
     }
     
     return(pl.pp)
   }
   
   # Description: Function to generate grid values for plotting decision boundaries for modification 2
   #  modification: multiple splits
   # ppboundMOD(data , test, meth, entro , entroindiv, title, simM = FALSE, strule, tot)
   # data data frame with the simulated dataset 
   # test data frame simulated test data 
   # entro logical; if TRUE use entropy in the modified PPtreeExt
   # entroindiv logical; if TRUE use individual entropy stopping rule in PPtreeExt
   # meth character; method to use "MOD" for PPtreeExt with multiple splits
   # title character; title for the plot
   
   ppboundMOD <-
     function(data ,
              test,
              meth = "MOD",
              entro = FALSE,
              entroindiv = TRUE,
              title,
              simM = FALSE,
              strule,
              tot) {
       
       # Generate grid values to evaluate tree
       grilla <-
         base::expand.grid(
           X1 = seq((min(data$X1) + sign(min(data$X1)) * .5), (max(data$X1) + sign(max(data$X1)) * .5), length.out = 100), 
           X2 = seq((min(data$X2) + sign(min(data$X2)) * .5), (max(data$X2) + sign(max(data$X2)) * .5), length.out = 100)
         )
       
       # Sim variable must be a factor 
       data$Sim <- as.factor(data$Sim)
       pptree <-
         PPtreeExt::PPtreeExtclass(
           Sim ~ . ,
           data = data,
           PPmethod = 'LDA',
           strule = strule,
           tot = tot
         )
       
       ppred.sim <- predict(object = pptree, newdata = grilla)
       
       grilla$ppred <- ppred.sim[[2]]
       
       err <-
         round(predict(object = pptree, newdata = test[, -1], true.class = test[, 1])[[1]] /
                 nrow(test[, -1]),
               3) * 100
       
       if (simM) {
         pl.pp <-
           ggplot2::ggplot(data = grilla) + ggplot2::geom_point(ggplot2::aes(x = X1, y = X2, color = ppred), alpha = .20) +
           ggplot2::scale_colour_brewer(name = "Class",
                                        type = "qual",
                                        palette = "Dark2") +
           ggplot2::theme_bw() +
           ggplot2::scale_shape_discrete(name = 'Class') + 
           ggplot2::geom_point(
             data = data,
             ggplot2::aes(
               x = X1 ,
               y = X2,
               group = Sim,
               color = Sim
             ),
             size = I(3)
           ) +
           ggplot2::theme(legend.position = "none", aspect.ratio = 1) +
           ggplot2::scale_y_continuous(expand = c(0, 0)) +
           ggplot2::scale_x_continuous(expand = c(0, 0)) +
           ggplot2::labs(
             x = " ",
             y = "",
             title = paste(title, "\n(test error: ", err, "%)", sep = '')
           )
         
       } else {
         pl.pp <-
           ggplot2::ggplot(data = grilla) + 
           ggplot2::geom_point(ggplot2::aes(
             x = X1,
             y = X2,
             color = ppred,
             shape = ppred
           ), 
           alpha = .20) +
           ggplot2::scale_colour_brewer(name = "Class",
                                        type = "qual",
                                        palette = "Dark2") + ggplot2::theme_bw() +
           ggplot2::scale_shape_discrete(name = 'Class') + 
           ggplot2::geom_point(
             data = data,
             ggplot2::aes(
               x = X1 ,
               y = X2,
               group = Sim,
               shape = Sim,
               color = Sim
             ),
             size = I(3)
           ) +
           ggplot2::theme(legend.position = "none", aspect.ratio = 1) +
           ggplot2::scale_y_continuous(expand = c(0, 0)) +
           ggplot2::scale_x_continuous(expand = c(0, 0)) +
           ggplot2::labs(
             x = " ",
             y = "",
             title = paste(title, "\n(test error: ", err, "%)", sep = '')
           )
         
       }
       return( pl.pp)
     }
   
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
               choices = c("Subsetting clases" = "1", 
                           
                           "Multiple splits" = "3"),
               selected = 3
             )
           )
         ),
         shiny::fluidRow(
           shiny::column(
             4,
             shiny::textInput(
               inputId = 'mean',
               label = 'Group means ',
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
           ) ,
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
               choices = c("Subsetting clases" = "1", 
                           "Multiple splits" = "3"),
               selected = 3
             )
           )
         ),
         shiny::fluidRow(
           shiny::column(
             4,
             shiny::textInput(
               inputId = 'mean2',
               label = 'Group means ',
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
           ) ,
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
               inputId = 'meanout',
               label = 'Out. X1, X2 means ',
               value =
                 "-3, 3"
             )
           ),
           shiny::column(
             4,
             shiny::textInput(
               inputId = 'sdout',
               label = 'Out. X1, X2 sd ',
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
               choices = c("Subsetting clases" = "1", 
                           "Multiple splits" = "3"),
               selected = 3
             )
           )
         ),
         shiny::fluidRow(
           shiny::column(4, shiny::numericInput(
             "size", label = "Sample size", value = 500
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
         #shiny::numericInput("p", label = "number of dimensions", value = 5),
         shiny::fluidRow(shiny::actionButton("simmaitra", "OK")),
         shiny::fluidRow(shiny::plotOutput("plotsmaitra"))
         ##
         
       ),
       shiny::tabPanel(
         "Draw Data",
         shiny::fluidRow(
           shiny::column(4, shiny::radioButtons("drawClass", "Select class:", choices = c("sim1", "sim2", "sim3"), inline = TRUE)),
           shiny::column(2, shiny::selectInput("drawRule", "Rule", choices = 1:8, selected = 1)),
           shiny::column(3, shiny::selectInput("drawModi", "Modification", choices = c("Subsetting" = "1", "Multiple splits" = "3"), selected = 3))
         ),
         shiny::fluidRow(
           shiny::column(2, shiny::actionButton("drawUndo", "Undo Last")),
           shiny::column(2, shiny::actionButton("drawClear", "Clear All")),
           shiny::column(2, shiny::actionButton("drawRun", "Run Classifiers"))
         ),
         shiny::fluidRow(
           shiny::column(12, shiny::textOutput("drawCount"))
         ),
         shiny::tags$hr(),
         shiny::fluidRow(shiny::plotOutput("drawCanvas", click = "canvas_click", height = "380px")),
         shiny::fluidRow(style = "margin-top: -8px;", shiny::uiOutput("drawBoundariesUi"))
       )
     )
   ))
   
   
   # Server ------------------------------------------------------------------
   server <- function(input, output) {
     output$distPlot <- shiny::renderPlot({
       if (input$do) {
         x1 <-shiny::isolate(as.numeric(unlist(strsplit(input$mean, ","))))
         x2 <-shiny::isolate(as.numeric(unlist(strsplit(input$cor, ","))))
         x3 <-shiny::isolate(as.numeric(unlist(strsplit(input$sample, ","))))
         x4 <- shiny::isolate(as.numeric(input$stop))
         dat.pl2 <-
           shiny::isolate(simu3(x1[1], x1[2], x1[3], x1[4], x1[5], x1[6],
                                x2[1], x2[2], x2[3], x3[1], x3[2], x3[3]))
         dat.test <-
           shiny::isolate(simu3(x1[1], x1[2], x1[3], x1[4], x1[5], x1[6],
                                x2[1], x2[2], x2[3], round(x3[1]*0.25), round(x3[2]*0.25), round(x3[3]*0.25)))
         
         if (input$modi == 1) {
           modpl <-
             ppbound(
               ru = 1,#as.numeric(input$rule),
               data = dat.pl2,
               test = dat.test,
               meth = "Modified" ,
               entro = FALSE,
               title = "PPtreeExt: Subsetting clases"
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
             ppboundMOD(
               data = dat.pl2,
               test = dat.test,
               meth = "MOD",
               entro = FALSE,
               entroindiv = TRUE,
               title = "PPtreeExt: Multiple splits",
               strule = x4,
               tot = sum(x3)
             )
         }
         
         gridExtra::grid.arrange(
           ppbound(
             ru =  as.numeric(input$rule),
             data = dat.pl2,
             test = dat.test,
             meth = "Rpart",
             entro = TRUE ,
             title = "Rpart"
           ),
           ppbound(
             ru =  as.numeric(input$rule),
             data = dat.pl2,
             test =  dat.test,
             meth = "Original" ,
             entro = FALSE,
             title = "PPtree"
           ),
           #ppbound(ru =  as.numeric(input$rule),  data = dat.pl2, meth = "Modified" , entro = TRUE),
           modpl,
           ppbound(
             ru = as.numeric(input$rule),
             data = dat.pl2,
             test = dat.test,
             meth = "RandomForest",
             entro = FALSE,
             title = "Random Forest"
           ),
           ncol = 4
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
         dat.pl2 <- simu3(x1[1], x1[2], x1[3], x1[4], x1[5], x1[6],
                          x2[1], x2[2], x2[3], x3[1], x3[2], x3[3])
         dat.test <-
           shiny::isolate(simu3(x1[1], x1[2], x1[3], x1[4], x1[5], x1[6],
                                x2[1], x2[2], x2[3], round(x3[1]*0.25), round(x3[2]*0.25), round(x3[3]*0.25)))
         
         aux <-
           data.frame(
             Sim = rep(paste("sim", as.numeric(input$group), sep = ""), x6),
             X1 = stats::rnorm(n = x6,
                               mean = x4[1],
                               sd = x5[1]),
             X2 = stats::rnorm(n = x6,
                               mean = x4[2],
                               sd = x5[2])
           )
         
         aux2 <-
           data.frame(
             Sim = rep(paste("sim", as.numeric(input$group), sep = ""), round(x6*0.25)),
             X1 = stats::rnorm(n = round(x6*0.25),
                               mean = x4[1],
                               sd = x5[1]),
             X2 = stats::rnorm(n = round(x6*0.25),
                               mean = x4[2],
                               sd = x5[2])
           )
         dat.pl2 <- rbind(dat.pl2, aux)
         dat.test <- rbind(dat.test, aux2)
         
         if (input$modi2 == 1) {
           modpl <-
             ppbound(
               ru =  as.numeric(input$rule),
               data = dat.pl2,
               test = dat.test,
               meth = "Modified" ,
               entro = FALSE,
               title = "PPtreeExt: Subsetting clases"
             )
         }
         if (input$modi2 == 2) {
           modpl <-
             ppbound(
               ru =  as.numeric(input$rule),
               data = dat.pl2,
               test = dat.test,
               meth = "Modified" ,
               entro = TRUE,
               title = "Modified 2"
             )
         }
         if (input$modi2 == 3) {
           modpl <-
             ppboundMOD(
               data = dat.pl2,
               test = dat.test,
               meth = "MOD",
               entro = FALSE,
               entroindiv = TRUE,
               title = "PPtreeExt: Multiple splits",
               strule = x7,
               tot = sum(x3 + x6)
             )
         }
         
         gridExtra::grid.arrange(
           ppbound(
             ru =  as.numeric(input$rule2),
             data = dat.pl2,
             test = dat.test,
             meth = "Rpart" ,
             entro = FALSE,
             title = "Rpart"
           ),
           ppbound(
             ru =  as.numeric(input$rule2),
             data = dat.pl2,
             test = dat.test,
             meth = "Original",
             entro = TRUE,
             title = "PPtree"
           ),
           
           # ppbound(ru =  as.numeric(input$rule2), FALSE, data = dat.pl2, meth = "Modified" , entro = TRUE),
           modpl,
           ppbound(
             ru = as.numeric(input$rule2),
             data = dat.pl2,
             test = dat.test,
             meth = "RandomForest",
             entro = FALSE,
             title = "Random Forest"
           ),
           ncol = 4
         )
         
         
       }
     })
     
     
     output$plotsmaitra <- shiny::renderPlot({
       if (input$simmaitra) {
         x1 <- shiny::isolate(as.numeric(input$stop3))
         #Q <- MixSim(BarOmega = 0.01, K = 4, p = 2)
         
         repeat {
           Q <-
             MixSim::MixSim(
               BarOmega = shiny::isolate(as.numeric(input$BarOmega)),
               MaxOmega = shiny::isolate(as.numeric(input$MaxOmega)),
               K = shiny::isolate(as.numeric(input$K)),
               p = 2,
               # sph = FALSE,
               #ecc = 1,
               # PiLow = 1.0,
               # int = c(0.0, 1.0),
               # resN = 100,
               # eps = 1e-06,
               # lim = 1e06
             )
           if (Q$fail == 0)
             break
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
             n = shiny::isolate(as.numeric(round(input$size*0.25))),
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
             ppbound(
               ru =  as.numeric(input$rule3),
               data = dat.pl2,
               test =dat.test,
               meth = "Modified" ,
               entro = FALSE,
               title = "PPtreeExt: Subsetting clases",
               simM = TRUE
             )
         }
         if (input$modi3 == 2) {
           modpl <-
             ppbound(
               ru =  as.numeric(input$rule3),
               data = dat.pl2,
               test = dat.test,
               meth = "Modified" ,
               entro = TRUE,
               title = "Modified 2",
               simM = TRUE
             )
         }
         if (input$modi3 == 3) {
           modpl <-
             ppboundMOD(
               data = dat.pl2,
               test = dat.test,
               meth = "MOD",
               entro = FALSE,
               entroindiv = TRUE,
               title = "PPtreeExt: Multiple splits",
               simM = TRUE,
               strule = x1,
               tot = input$size
             )
         }
         
         gridExtra::grid.arrange(
           ppbound(
             ru = as.numeric(input$rule3),
             data = dat.pl2,
             test = dat.test,
             meth = "Rpart",
             entro = TRUE ,
             title = "Rpart",
             simM = TRUE
           ),
           ppbound(
             ru = as.numeric(input$rule3),
             data = dat.pl2,
             test = dat.test,
             meth = "Original" ,
             entro = FALSE,
             title = "PPtree",
             simM = TRUE
           ),
           #ppbound(ru =  as.numeric(input$rule),  data = dat.pl2, meth = "Modified" , entro = TRUE),
           modpl,
           ppbound(
             ru = as.numeric(input$rule3),
             data = dat.pl2,
             test = dat.test,
             meth = "RandomForest",
             entro = FALSE,
             title = "Random Forest",
             simM = TRUE
           ),
           ncol = 4
         )
       }
     })
     
    # Keep track of points
    drawn_data <- shiny::reactiveVal(data.frame(Sim=character(), X1=numeric(), X2=numeric()))
    draw_run <- shiny::reactiveVal(0)
     
     # Catch clicks
     shiny::observeEvent(input$canvas_click, {
       click <- input$canvas_click
       if (!is.null(click)) {
         drawn_data(rbind(drawn_data(), data.frame(Sim = input$drawClass, X1 = click$x, X2 = click$y)))
       }
     })
     
     shiny::observeEvent(input$drawUndo, {
       d <- drawn_data()
       if (nrow(d) > 0) drawn_data(d[-nrow(d), , drop = FALSE])
     })
     
     shiny::observeEvent(input$drawClear, {
       drawn_data(data.frame(Sim=character(), X1=numeric(), X2=numeric()))
      draw_run(0)
     })

    shiny::observeEvent(input$drawRun, {
      draw_run(draw_run() + 1)
    })
     
     output$drawCount <- shiny::renderText({ paste("Points:", nrow(drawn_data())) })
     
     output$drawCanvas <- shiny::renderPlot({
       d <- drawn_data()
       legend_levels <- c("sim1", "sim2", "sim3")
       legend_df <- data.frame(
         X1 = rep(0, length(legend_levels)),
         X2 = rep(0, length(legend_levels)),
         Sim = factor(legend_levels, levels = legend_levels)
       )
       p <- ggplot2::ggplot() + ggplot2::theme_bw() +
            ggplot2::labs(title = "Click to place points") +
            ggplot2::geom_point(
              data = legend_df,
              ggplot2::aes(x = X1, y = X2, color = Sim, shape = Sim),
              alpha = 0
            ) +
            ggplot2::scale_colour_brewer(
              palette = "Dark2",
              drop = FALSE,
              guide = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 3))
            ) +
            ggplot2::scale_shape_discrete(drop = FALSE) +
            ggplot2::coord_fixed(ratio = 1, xlim = c(-5, 5), ylim = c(-5, 5), expand = FALSE) +
            ggplot2::theme(
              legend.position = "right",
              plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
            )

       if (nrow(d) > 0) {
         p <- p + ggplot2::geom_point(data=d, ggplot2::aes(x=X1, y=X2, color=Sim, shape=Sim), size=3)
       }
       p
     })

     output$drawBoundariesUi <- shiny::renderUI({
      if (draw_run() == 0) return(NULL)
      shiny::plotOutput("drawBoundaries", height = "460px")
     })
     
     output$drawBoundaries <- shiny::renderPlot({
       shiny::req(draw_run() > 0)

       d <- shiny::isolate(drawn_data())
       
       # Validation
       if (nrow(d) < 5) return(plot(0, main = "Need at least 5 points"))
       if (length(unique(d$Sim)) < 2) return(plot(0, main = "Need at least 2 different classes"))

       # FIX 1: Explicitly set all factor levels (sim1, sim2, sim3)
       d$Sim <- factor(d$Sim, levels = c("sim1", "sim2", "sim3"))
       
       # FIX 2: Ensure numeric columns are actually numeric
       d$X1 <- as.numeric(d$X1)
       d$X2 <- as.numeric(d$X2)
       
       # FIX 3: Expand data range slightly to avoid grid collapse
       expand_range <- function(x, pct = 0.1) {
         r <- range(x, na.rm = TRUE)
         w <- max(abs(r)) * pct
         c(r[1] - w, r[2] + w)
       }
       
       # Create test data with expanded range for grid generation
       test_data <- data.frame(
         Sim = d$Sim[1],
         X1 = seq(expand_range(d$X1)[1], expand_range(d$X1)[2], length.out = 10),
         X2 = seq(expand_range(d$X2)[1], expand_range(d$X2)[2], length.out = 10)
       )
       
       tryCatch({
         # Build plots with error handling
         rule <- as.numeric(shiny::isolate(input$drawRule))
         
         # Helper function to create error ggplot
         error_plot <- function(msg) {
           ggplot2::ggplot() + 
             ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg, size = 5) +
             ggplot2::theme_void() +
             ggplot2::labs(title = msg)
         }
         
         # RandomForest specifically fails if any factor level has 0 points
         # So we need to drop empty levels just for the data we feed to it, or it throws
         # "Can't have empty classes in y."
         
         p1 <- tryCatch(
           ppbound(rule, data = d, test = d, meth = "Rpart", entro = FALSE, title = "Rpart"),
           error = function(e) error_plot("Rpart failed")
         )
         
         p2 <- tryCatch(
           ppbound(rule, data = d, test = d, meth = "Original", entro = FALSE, title = "PPtree"),
           error = function(e) error_plot("PPtree failed")
         )
         
         if (shiny::isolate(input$drawModi) == "1") {
           p3 <- tryCatch(
             ppbound(rule, data = d, test = d, meth = "Modified", entro = FALSE, title = "PPtreeExt: Sub"),
             error = function(e) error_plot("PPtreeExt Sub failed")
           )
         } else {
           p3 <- tryCatch(
             ppboundMOD(data = d, test = d, meth = "MOD", entro = FALSE, entroindiv = TRUE, 
                        title = "PPtreeExt: Mul", strule = rule, tot = nrow(d)),
             error = function(e) error_plot("PPtreeExt Mul failed")
           )
         }
         
         p4 <- tryCatch({
           # Random forest needs dropping empty levels or it fails
           d_rf <- d
           d_rf$Sim <- droplevels(d_rf$Sim)
           ppbound(rule, data = d_rf, test = d_rf, meth = "RandomForest", entro = FALSE, title = "Random Forest")
         }, error = function(e) error_plot("RandomForest failed")
         )
         
         gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 4)
         
       }, error = function(e) {
         plot(0, main = paste("Error:", e$message))
       })
     })
     
   }
   
   shiny::shinyApp(ui = ui, server = server)
 }
 
 