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
source("./modules/plotModules/dataviewTSNEPlotModule.R", encoding = "UTF-8")
source("./modules/plotModules/dataviewTSNEPlotModule2.R", encoding = "UTF-8")
source("./modules/plotModules/dataviewTSNEPlotModule3.R", encoding = "UTF-8")

ui = fluidPage(
    fillRow(
        height = 650,
        examplePlotModuleUI("example"),                
        dataviewtSNEModuleUI("tSNEPlot"),
        dataviewtSNEplotModuleUI2("tSNEPlot2"),
        dataviewtSNEplotModuleUI3("tSNEPlot3", height=c(600,800)*1)
    )
)

server = function(input, output, session) {

    colnames(ngs$Y)  ## possible group vectors

    examplePlotModuleServer("example")
    
    filterStates <- list(
        search_gene = "ETAA1",
        data_samplefilter = NULL,
        data_type = "logCPM",
        data_groupby = "condition"
    )

    dataviewtSNEModuleServer(
        "tSNEPlot", filterStates, ngs
    )
    
    dataviewtSNEplotModuleServer2(
        "tSNEPlot2", filterStates, ngs,
        label="A", imgH=600, watermark=FALSE
    )
    
    dataviewtSNEplotModuleServer3(
        "tSNEPlot3", filterStates, ngs,
        label="A", watermark=FALSE
    )
}

shinyApp(ui, server=server, options=list(launch.browser=TRUE))
