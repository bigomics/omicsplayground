##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2021 BigOmics Analytics Sagl. All rights reserved.
##

source("./global.R")  ## global variable

load("../data/example-data.pgx",verbose=1)
ngs = pgx.initialize(ngs)

ui = fluidPage(
    tabView("DataView",DataViewInputs("view"),DataViewUI("view"))
)

server = function(input, output, session) {

    filterStates <- list(
        search_gene = "ETAA1",
        data_samplefilter = NULL,
        data_type = "logCPM",
        data_groupby = "condition"
    )

    inputData <- shiny::reactive(ngs)
    shiny::callModule(DataViewBoard, "view", inputData, filterStates)
}

shinyApp(ui, server=server, options=list(launch.browser=TRUE))

