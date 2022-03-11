##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing IntersectionBoard")

IntersectionInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

IntersectionUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Pairwise scatter",uiOutput(ns("scatterPlotMatrix_UI"))),
            shiny::tabPanel("Signature clustering",uiOutput(ns("ctClustering_UI")))
        )
    )
}