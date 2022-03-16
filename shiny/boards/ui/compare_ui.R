##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

CompareInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

CompareUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Compare", shiny::uiOutput(ns("compareScatter_UI"))),
            shiny::tabPanel("Foldchange", shiny::uiOutput(ns("FCcorrelation_UI"))),
            shiny::tabPanel("Gene Correlation", shiny::uiOutput(ns("GeneCorrelation_UI")))            
        )
    )
}