##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

DataViewInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

DataViewUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tabsetPanel(
        id = ns("tabs"),
        shiny::tabPanel("Plots",uiOutput(ns("plotsUI"))),
        shiny::tabPanel("QC",uiOutput(ns("countsUI"))),
        shiny::tabPanel("Counts",uiOutput(ns("genetableUI"))),
        shiny::tabPanel("Samples",uiOutput(ns("sampletableUI"))),
        shiny::tabPanel("Contrasts",uiOutput(ns("contrasttableUI"))),        
        shiny::tabPanel("Resource info",uiOutput(ns("resourceinfoUI")))
    )
}