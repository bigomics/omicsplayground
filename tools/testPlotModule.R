##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2021 BigOmics Analytics Sagl. All rights reserved.
##

if(interactive()) {

    library(shiny)
    ## RUN FROM root folder!
    source("shiny/global.R")  ## global variable
    source("shiny/modules/plotModules/examplePlotModule.R")  ## example
    
    load("data/example-data.pgx", verbose = 1)
    ##load("../data/GSE72056-scmelanoma.pgx", verbose = 1)
    ngs <- pgx.initialize(ngs)
    
    ui = fluidPage(
        fillRow(
            height = 600,
            examplePlotModuleUI("example"),
            dataviewTSNEPlotModuleUI("tsne1", height = c(600,800)),
            dataviewTSNEPlotModuleUI("tsne2", height = c(600,800))
        )
    )
    
    server = function(input, output, session) {
        
        ## examplePlotModuleServer("example")
        
        filterStates <- list(
            search_gene = "ETAA1",
            data_samplefilter = NULL,
            data_type = "logCPM",
            data_groupby = "condition"
        )

        dataviewTSNEPlotModuleServer(
            "tsne1", reactive(ngs), filterStates, 
        watermark = FALSE
        )
        
        filterStates2 <- list(
            search_gene = "ETAA1",
            data_samplefilter = NULL,
            data_type = "logCPM",
            data_groupby = "activated"
        )
        
        dataviewTSNEPlotModuleServer(
            "tsne2", reactive(ngs), filterStates2, 
            watermark = FALSE
        )    
    }
    
    ##shinyApp(ui, server = server, options = list(launch.browser = TRUE))
    shinyApp(ui, server = server)

}
