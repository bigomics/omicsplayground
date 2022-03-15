##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

WgcnaInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

WgcnaUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    ui <- shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("WGCNA",uiOutput(ns("wgcna_UI"))),
            shiny::tabPanel("Modules",uiOutput(ns("modules_UI"))),
            shiny::tabPanel("Eigengenes",uiOutput(ns("eigen_UI"))),
            shiny::tabPanel("Intramodular",uiOutput(ns("intra_UI")))
        )
    )
    ui
}