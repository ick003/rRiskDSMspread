# Define server logic required to draw a histogram ----
server <- function(input, output) {
    
    # Histogram of the Old Faithful Geyser Data ----
    # with requested number of bins
    # This expression that generates a histogram is wrapped in a call
    # to renderPlot to indicate that:
    #
    # 1. It is "reactive" and therefore should be automatically
    #    re-executed when inputs (input$bins) change
    # 2. Its output type is a plot
    grid2D <- 
        reactive({

        spreadR::prepare.grid( covariates=list(MRR_swarm[,c("X_m","Y_m")]),
                            N=c(input$grid_x,input$grid_y), xlim=xlimmy, ylim=ylimmy,
                            inner.xlim=inner.xlimmy, inner.ylim=inner.ylimmy)
        })
    
    simRes <- reactive(NULL)
    
    simRes <- 
        eventReactive(input$run_sim, {
            
            nSim = 1
            paramL = list()
            for(j in 1:nSim){
                tmp = list(delta=input$delta,alpha=input$alpha,nu=0,mu=input$mu,beta=1)
                paramL[[length(paramL)+1]] = tmp
                
            }
            
            prep.covars <- spreadR::prepare.covariates( grid=grid2D(),
                                               covariates=list(as.data.frame(MRR_swarm[,c("X_m","Y_m")])),
                                               sigma= input$sigma)
            
            lastDay = input$duration
            times = seq(0,lastDay,1)

            simM <- spreadR::simMozzies(params = paramL[[1]],releaseSite = grid2D_release_sites[as.numeric(input$release_site),],
                               grid=grid2D(), covariates=prep.covars, times=times, 
                               releaseSize = input$size, boundaryCond = "Neumann")
            idxKeep <- round(seq.int(1,length(times),length.out = 4))
            
            simM[idxKeep,]
        })
    
    
    output$distPlot1 <- renderPlot({
        if(input$which_plot == "grid"){
            plot( MRR_swarm[,c("X_m","Y_m")], pch=20, cex=0.5,
                  xlab='', ylab='', xlim = range(grid2D()$x.int), ylim = range(grid2D()$y.int))
            points(MRR_compounds[,c("X_m","Y_m")], pch=0, cex=1)
        abline( v=grid2D()$x.int, lwd=0.3, col='red')
        abline( h=grid2D()$y.int, lwd=0.3, col='red')
        rect(min(inner.xlimmy), min(inner.ylimmy), max(inner.xlimmy), max(inner.ylimmy), border = 'blue')
        legend("bottomright", legend=c("Swarms", "Compounds","Cell boundaries", "Inner Region"),horiz = TRUE,
               pch=c(20,0,  NA, NA), cex=0.65, lty=c(NA, NA, 1, 1), col=c("black","black","red","blue"), bg="white")
        }
        if(input$which_plot == "spread"){
            if(!is.null(simRes())){
                spreadR::mozImage( soln=simRes(), grid2D=grid2D(), pointsToAdd = MRR_compounds[,c("X_m","Y_m")], xlim=c( grid2D()$x.up, grid2D()$x.down),
                          ylim=c( grid2D()$y.up, grid2D()$y.down), mfrow=  c(2,2))
            }
        }
    })
    
    ranges2 <- reactiveValues(x = NULL, y = NULL)
    
    observe({
        brush <- input$plot1_brush
        if (!is.null(brush)) {
            ranges2$x <- c(brush$xmin, brush$xmax)
            ranges2$y <- c(brush$ymin, brush$ymax)
            
        } else {
            ranges2$x <- NULL
            ranges2$y <- NULL
        }
    })
    
    
    output$distPlot2 <- renderPlot({
        if(input$which_plot == "grid"){
            plot( MRR_swarm[,c("X_m","Y_m")], pch=20, cex=0.5,
                  xlab='', ylab='', xlim = ranges2$x, ylim = ranges2$y, 
                  sub="Define a zone on the left plot to zoom in on the right plot")
            points(MRR_compounds[,c("X_m","Y_m")], pch=0, cex=1)
            abline( v=grid2D()$x.int, lwd=0.3, col='red')
            abline( h=grid2D()$y.int, lwd=0.3, col='red')
            rect(min(inner.xlimmy), min(inner.ylimmy), max(inner.xlimmy), max(inner.ylimmy), border = 'blue')
        }
        if(input$which_plot == "spread"){
            if(!is.null(simRes())){
                if(!is.null(ranges2$x)){
                spreadR::mozImage( soln=simRes(), grid2D=grid2D(), xlim=ranges2$x, ylim=ranges2$y,pointsToAdd = MRR_compounds[,c("X_m","Y_m")], mfrow=  c(2,2), 
                                   sub="Define a zone on the left plot to zoom in on the right plot")} else{
                    spreadR::mozImage( soln=simRes(), grid2D=grid2D(), 
                                       xlim=c( grid2D()$x.up, grid2D()$x.down),
                                       ylim=c( grid2D()$y.up, grid2D()$y.down), 
                                       pointsToAdd = MRR_compounds[,c("X_m","Y_m")], 
                                       mfrow=  c(2,2), sub="Define a zone on the left plot to zoom in on the right plot")
                }
            }
        }
    })
    
}