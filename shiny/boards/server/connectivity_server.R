##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

ConnectivityInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

ConnectivityUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("FC correlation", shiny::uiOutput(ns("cmapCorrelation_UI"))),
            shiny::tabPanel("FC heatmap", shiny::uiOutput(ns("cmapHeatmap_UI"))),            
            shiny::tabPanel("Meta-graph", shiny::uiOutput(ns("cmapMetaAnalysis_UI"))),
            shiny::tabPanel("Experiment clustering", shiny::uiOutput(ns("cmapClustering_UI")))
        )
    )
}