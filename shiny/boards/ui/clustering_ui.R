##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

ClusteringInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        tags$div(
         HTML("<h3>Clustering Analysis</h3> Discover clusters of similar genes or samples using unsupervised machine learning.")
        ),
        shiny::tagList(
            shinyBS::tipify( shiny::actionLink(ns("clust_info"), "Tutorial", icon = shiny::icon("youtube")),
                   "Show more information and video tutorial about this module."),
            shiny::hr(), shiny::br(),             
            shinyBS::tipify( shiny::selectInput(ns("hm_features"),"Features:", choices=NULL, multiple=FALSE),
                   "Select a family of features.", placement="top"),
            shiny::conditionalPanel(
                "input.hm_features == '<custom>'", ns=ns,
                shinyBS::tipify( shiny::textAreaInput(ns("hm_customfeatures"), NULL, value = NULL,
                                      height = "150px", width = "100%", 
                                      rows=5, placeholder="Paste your custom gene list"),
                       "Paste a custom list of genes to be used as features.",
                       placement="bottom")
            ),
            shiny::conditionalPanel(
                "input.hm_features == '<contrast>'", ns=ns,
                tipifyR( shiny::selectInput(ns("hm_contrast"), NULL, choices=NULL),
                        "Select contrast to be used as signature.")
            ),
            shinyBS::tipify( shiny::selectInput(ns('hm_group'),'Group by:',choices=NULL),
                   "Group the samples by condition.", 
                   placement="top", options = list(container = "body")),
            shinyBS::tipify( shiny::selectInput(ns("hm_samplefilter"),"Filter samples:",
                                choices=NULL, multiple=TRUE),
                   "Filter the relevant samples for the analysis.",
                   placement="top", options = list(container = "body")),            
            shinyBS::tipify( shiny::actionLink(ns("hm_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Toggle advanced options.", placement="top"),
            shiny::br(),
            shiny::conditionalPanel(
                "input.hm_options % 2 == 1", ns=ns,
                shiny::tagList(
                           shinyBS::tipify( shiny::selectInput(ns("hm_level"),"Level:", choices=c("gene","geneset")),
                                           "Specify the level analysis: gene or geneset level.",
                                           placement="top", options = list(container = "body")),
                           shinyBS::tipify( shiny::checkboxInput(ns('hm_filterXY'),'exclude X/Y genes',FALSE),
                                           "Exclude genes on X/Y chromosomes.", 
                                           placement="top", options = list(container = "body")),
                           shinyBS::tipify( shiny::checkboxInput(ns('hm_filterMitoRibo'),
                                                                 'exclude mito/ribo genes',FALSE),
                                           "Exclude mitochondrial (MT) and ribosomal protein (RPS/RPL) genes.", 
                                           placement="top", options = list(container = "body"))
                    )
            )
        )
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
            shiny::tabPanel("Annotate clusters",
                uiOutput(ns("hm_annotateUI"))),
            shiny::tabPanel("Phenotypes",uiOutput(ns("hm_phenoplotUI"))),
            shiny::tabPanel("Feature ranking",uiOutput(ns("hm_featurerankUI")))      
        )
    )
}
