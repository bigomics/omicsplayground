##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2021 BigOmics Analytics Sagl. All rights reserved.
##

library(shiny)

## RUN FROM root folder!
setwd(pkgload::pkg_path())
source("shiny/global.R")  ## global variable
source("R/00Headers.R")  ## global variable

load("data/example-data.pgx",verbose=1)
ngs = pgx.initialize(ngs)

ui = fluidPage(
    tabView("DataView", DataViewInputs("view"), DataViewUI("view"))
    ##tabView("Cluster samples",ClusteringInputs("clust"),ClusteringUI("clust"))
)

server = function(input, output, session) {
    callModule(DataViewBoard, "view", reactive(ngs))
    ##callModule( ClusteringBoard, "clust", reactive(ngs))    
}

shinyApp(ui, server=server, options=list(launch.browser=TRUE))

