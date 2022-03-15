##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

SignatureInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

SignatureUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillRow(
        flex = c(1.5,0.05,1),
        height = 780,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Enrichment",uiOutput(ns("enplots_UI"))),
            shiny::tabPanel("Volcano plots",uiOutput(ns("volcanoPlots_UI"))),
            shiny::tabPanel("Overlap/similarity",uiOutput(ns("overlapAnalysis_UI"))),
            shiny::tabPanel("Markers",uiOutput(ns("markers_UI")))
        ),
        shiny::br(),
        shiny::tabsetPanel(
            id = ns("tabs2"),
            shiny::tabPanel("Enrichment table",uiOutput(ns("enrichmentTables_UI")))
        )
    )
}