##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing EnrichmentBoard")

EnrichmentInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

EnrichmentUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        flex = c(1.75,1),
        height = 800,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Top enriched",uiOutput(ns("topEnriched_UI"))),
            shiny::tabPanel("Plots",uiOutput(ns("subplots_UI"))),
            shiny::tabPanel("Compare",uiOutput(ns("compare_UI"))),
            shiny::tabPanel("Volcano (all)",uiOutput(ns("volcanoAll_UI"))),
            shiny::tabPanel("Volcano (methods)",uiOutput(ns("volcanoMethods_UI")))
        ),
        shiny::tabsetPanel(
            id = ns("tabs2"),
            shiny::tabPanel("Table",uiOutput(ns("tables_UI"))),
            shiny::tabPanel("Foldchange (all)",uiOutput(ns("fctable_UI"))),
            shiny::tabPanel("FDR table",uiOutput(ns("FDRtable_UI")))                       
        )
    )
}
