##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

FeatureMapInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

FeatureMapUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    ui <- shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("Gene", shiny::uiOutput(ns("geneUMAP_UI"))),
            shiny::tabPanel("Geneset", shiny::uiOutput(ns("gsetUMAP_UI")))            
        )
    )
    ui
}