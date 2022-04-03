##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2021 BigOmics Analytics Sagl. All rights reserved.
##

source("global.R")  ## global variable
source("./modules/plotModules/dataviewTSNEPlotModule2.R", encoding = "UTF-8")

load("../data/example-data.pgx",verbose=1)
ngs = pgx.initialize(ngs)

ui = fluidPage(
    fillCol(
        fillRow(
            plotWidget("tSNEPlot"),
            br(),br(),br()
        ),
        fillRow(
            br(),br(),br(),br()
        )
    )
)

server = function(input, output, session) {

    colnames(ngs$Y)  ## possible group vectors

    filterStates <- list(
        search_gene = "ETAA1",
        data_samplefilter = NULL,
        data_type = "logCPM",
        data_groupby = "condition"
    )
    dataviewtSNEplotModule(
        "tSNEPlot", filterStates, ngs,
        label="A", imgH=400, watermark=FALSE
    )    
}

shiny::shinyApp(ui, server=server)

