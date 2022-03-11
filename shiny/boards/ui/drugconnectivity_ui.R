##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing DrugConnectivityBoard")

DrugConnectivityInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

DrugConnectivityUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        flex = c(1),
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("Drug enrichment",uiOutput(ns("DSEA_enrichment_UI"))),
            shiny::tabPanel("Connectivity map (beta)",uiOutput(ns("DSEA_cmap_UI")))
            ## shiny::tabPanel("Fire plot (dev)",uiOutput(ns("fireplot_UI")))            
        )
    )
}
