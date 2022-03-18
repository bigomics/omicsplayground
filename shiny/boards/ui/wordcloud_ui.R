##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

WordCloudInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

WordCloudUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        flex = c(1),
        height = 780,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("WordCloud",uiOutput(ns("wordcloud_UI")))
            ## shiny::tabPanel("Fire plot (dev)",uiOutput(ns("fireplot_UI")))            
        )
    )
}