##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

SingleCellInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

SingleCellUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        flex = c(1),
        height = 780,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("Cell type",uiOutput(ns("icp_UI"))),
            shiny::tabPanel("Mapping",uiOutput(ns("mapping_UI"))),
            shiny::tabPanel("Markers",uiOutput(ns("markersplot_UI"))),
            shiny::tabPanel("CNV",uiOutput(ns("cnaModule_UI"))),
            shiny::tabPanel("iTALK",uiOutput(ns("italk_panel_UI"))),
            shiny::tabPanel("Monocle",uiOutput(ns("monocle_panel_UI")))            
        )
    )
}
