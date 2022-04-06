
examplePlotModuleUI <- function(id) {
    ns <- shiny::NS(id)    
    options <- tagList(
        selectInput(ns("func"),"function",c("cos","sin","exp"))
    )    
    PlotModuleUI(
        NS(id,"plt"),   ## note nested NS
        options = options,
        ## width = 500,
        height = 800
    )

}

examplePlotModuleServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        
        plot_data <- shiny::reactive({
            ## prepare data
            input.func <- input$func ## we must use NS??
            x <- seq(1,10,0.1)
            y <- switch(
                input.func,
                cos = cos(x),
                sin = sin(x),
                exp = exp(x)
            )
            df <- data.frame(x=x, y=y)
            return(df)
        })
        
        require(pryr)
        ##plot.RENDER %<a-% reactive({
        plot.RENDER <- function(){    
            ## actual plotting
            df <- plot_data()
            plot(df$x, df$y)
        }

        ## nested module
        PlotModuleServer(
            "plt",
            func  = plot.RENDER,
            csvFunc = plot_data  ## downloadable data as CSV
        )        
    }
  )
}   
