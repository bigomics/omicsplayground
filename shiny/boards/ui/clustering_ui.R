##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

ClusteringInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

ClusteringUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillRow(
        flex = c(1.6,0.05,1),
        height = 780,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Heatmap",uiOutput(ns("hm_heatmap_UI"))),
            shiny::tabPanel("PCA/tSNE",uiOutput(ns("hm_pcaUI"))),
            shiny::tabPanel("Parallel",uiOutput(ns("hm_parcoordUI")))
        ),
        shiny::br(),
        shiny::tabsetPanel(
            id = ns("tabs2"),
            shiny::tabPanel("Annotate clusters",uiOutput(ns("hm_annotateUI"))),
            shiny::tabPanel("Phenotypes",uiOutput(ns("hm_phenoplotUI"))),
            shiny::tabPanel("Feature ranking",uiOutput(ns("hm_featurerankUI")))      
        )
    )
}
