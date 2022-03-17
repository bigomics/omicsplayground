##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

ConnectivityInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        tags$div(
         HTML("<b>Similar Experiments.</b> Find similar experiments by correlating their signatures.
          The main goal is to identify experiments showing similar signatures and find genes
          that are commonly up/down regulated between experiments.")
        ),
        shiny::tagList(
            shinyBS::tipify( shiny::actionLink(ns("cmap_info"), "Info", icon = shiny::icon("info-circle")),
                   "Show more information about this module"),
            shiny::hr(), shiny::br(),             
            shinyBS::tipify( shiny::selectInput(ns('cmap_contrast'),'Contrast:',
                                choices=NULL, multiple=FALSE),
                   "Select the contrast that you want to compare.",
                   placement="right", options = list(container = "body")
                   ),
            shinyBS::tipify( shiny::selectInput(ns('cmap_sigdb'),"Signature DB:", choices=NULL),
                   "Select reference signature database.",
                   placement="right", options = list(container = "body")),
            shiny::br(),
            shinyBS::tipify( shiny::actionLink(ns("cmap_options"), "Options", icon=icon("cog", lib="glyphicon")),
                   "Toggle advanced options.",
                   placement="right", options = list(container = "body")),
            shiny::br(),
            shiny::conditionalPanel(
                "input.cmap_options % 2 == 1", ns=ns,
                shinyBS::tipify(
                    shiny::checkboxInput( ns('cmap_hideclustcontrasts'),"hide cluster contrasts", TRUE),
                    "Hide cluster contrasts.",
                    placement="right", options = list(container = "body")),
                shinyBS::tipify(
                    shiny::checkboxInput(ns("cmap_abs_score"),"abs.score",TRUE),
                    "Use absolute score value",
                    placement="right", options = list(container = "body"))
            )
        )
    )
}

ConnectivityUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("FC correlation", shiny::uiOutput(ns("cmapCorrelation_UI"))),
            shiny::tabPanel("FC heatmap", shiny::uiOutput(ns("cmapHeatmap_UI"))),            
            shiny::tabPanel("Meta-graph", shiny::uiOutput(ns("cmapMetaAnalysis_UI"))),
            shiny::tabPanel("Experiment clustering", shiny::uiOutput(ns("cmapClustering_UI")))
        )
    )
}
