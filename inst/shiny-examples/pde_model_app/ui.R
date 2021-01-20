ui <- fluidPage(
    
    shiny::wellPanel(
        fluidRow(column(3, tags$h5("1. Define grid"),
                        column(6, numericInput("grid_x", label = "Number of horizontal cells",value = 100, min = 20, max = 200, step = 10),
                        numericInput("grid_y", label = "Number of vertical cells", value = 150, min = 20, max = 200, step = 10)), 
                        column(6, 
                               numericInput("distance_inner_boundary", label = "Inner boundary range (m)", value = 50, min = 0, max = 500, step = 50), 
                               numericInput("distance_outer_boundary", label = "Outer boundary range (m)", value = 250, min = 200, max = 1000, step = 100)
                               )),
                 column(6, tags$h5("2. Select parameters"),
                        fluidRow(column(4, numericInput(inputId = "sigma", label = "Sigma", value = 33, min = 1, max = 100, step = 1)),
                                 column(4, numericInput(inputId = "duration", label = "Duration (days)", value = 15, min = 1, max = 90, step = 1)),
                                 column(4, numericInput(inputId = "size", label = "Release size", value = 1500, min = 500, max = 10000, step = 500))),
                        fluidRow(column(4, numericInput(inputId = "delta", label = "Dispersal", value = 150, min = 10, max = 10000, step = 10)),
                                 column(4, numericInput(inputId = "alpha", label = "Attraction", value = 0.1, min = 0.01, max = 1, step = 0.01)),
                                 column(4, numericInput(inputId = "mu", label = "Mortality", value = 0.2, min = 0.1, max = 10, step = .1))),
                        fluidRow(
                            column(4, pickerInput("release_site", label = "Release site", choices = c(1:4)))
                        )
                 ),
                 column(3, tags$h5("3. Execute"),
                        actionButton("run_sim", label = "Run simulation"))
    )
    ),
    shiny::wellPanel(
        fluidRow(radioButtons("which_plot", "Plot type:",
                     c("Grid" = "grid",
                       "Spread map" = "spread",
                       "Population decay" = "decay"))),
        fluidRow(
            column(6,plotOutput(outputId = "distPlot1", brush = brushOpts(
                id = "plot1_brush",
                resetOnNew = TRUE
            ))),
            column(6,plotOutput(outputId = "distPlot2"))
            )
    )
    
)