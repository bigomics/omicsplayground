##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing ExpressionBoard")

ExpressionInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

ExpressionUI.test <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tabsetPanel(
        shiny::tabPanel("Table",uiOutput(ns("tables_UI"))),
        shiny::tabPanel("Foldchange (all)",uiOutput(ns("fctable_UI"))),
        shiny::tabPanel("FDR table",uiOutput(ns("FDRtable_UI")))                       
    )
}

ExpressionUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        flex = c(1.5,1),
        height = 780,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Plot",uiOutput(ns("plots_UI"))),
            shiny::tabPanel("Top genes",uiOutput(ns("topgenesUI"))),
            shiny::tabPanel("Volcano (all)",uiOutput(ns("volcanoAll_UI"))),
            ## shiny::tabPanel("Volcano (all2)",uiOutput(ns("volcanoAll2_UI"))),
            shiny::tabPanel("Volcano (methods)",uiOutput(ns("volcanoMethodsUI")))
        ),
        shiny::tabsetPanel(
            id = ns("tabs2"),
            shiny::tabPanel("Table",uiOutput(ns("tables_UI"))),
            shiny::tabPanel("Foldchange (all)",uiOutput(ns("fctable_UI"))),
            shiny::tabPanel("FDR table",uiOutput(ns("FDRtable_UI")))                       
        )
    )
}
