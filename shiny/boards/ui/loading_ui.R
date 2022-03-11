##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing LoadingBoard")

LoadingInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shinyBS::tipify( shiny::actionLink(ns("module_info"), "Tutorial", icon = shiny::icon("youtube")),
               "Show more information about this module.")
    )
}

LoadingUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace

    close_session <- shiny::span()
    if(getOption("OMICS_TEST", FALSE)){
        close_session <- shiny::actionButton(ns("close"), "close")
    }

    tagList(
        close_session,
        shiny::fillCol(
            height = 750,
            shiny::tabsetPanel(
                id = ns("tabs"),
                shiny::tabPanel("Datasets",uiOutput(ns("pgxtable_UI"))),
                shiny::tabPanel("Upload data",uiOutput(ns("upload_UI")))
            )
        )
    )
}
