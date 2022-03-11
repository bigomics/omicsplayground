##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing CorrelationBoard")

CorrelationInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

CorrelationUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    ui <- shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("Correlation",uiOutput(ns("corAnalysis_UI"))),
            shiny::tabPanel("Graph",uiOutput(ns("corGraph_UI")))            
            ## shiny::tabPanel("Functional",uiOutput(ns("corFunctional_UI"))),
            ## shiny::tabPanel("Differential",uiOutput(ns("corDiff_UI")))
        )
    )
    ui
}
