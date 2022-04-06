##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2021 BigOmics Analytics Sagl. All rights reserved.
##

library(shiny)
source("./global.R")  ## global variable
##source("./modules/plotModules/dataviewTSNEPlotModule.R", encoding = "UTF-8")
#source("./modules/plotModules/dataviewTSNEPlotModule2.R", encoding = "UTF-8")

load("../data/example-data.pgx",verbose=1)
ngs = pgx.initialize(ngs)

source("./modules/plotModules/PlotModule.R", encoding = "UTF-8")
source("./modules/plotModules/examplePlotModule.R", encoding = "UTF-8")

ui = fluidPage(
    fillRow(
        height = 800,
        examplePlotModuleUI("example"),
        br(), br()
    )
)

server = function(input, output, session) {
    ##examplePlotModuleServer("example", input)
    examplePlotModuleServer("example")

}

shinyApp(ui, server=server, options=list(launch.browser=TRUE))
