##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing BiomarkerModule")

BiomarkerInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

BiomarkerUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    ui <- shiny::fillCol(
        flex = c(1),
        height = 760,
        shiny::tabsetPanel(
            id=ns("tabs"),
            shiny::tabPanel("Importance", shiny::uiOutput(ns("pdx_biomarker_UI")))
        )
    )
    ui
}