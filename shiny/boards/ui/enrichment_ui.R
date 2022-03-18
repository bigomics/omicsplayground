##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

EnrichmentInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        tags$div(
         HTML("<b>Geneset enrichment analysis.</b> Perform differential expression analysis on a geneset level,
          also called geneset enrichment analysis.")
        ),
        shiny::tagList(
            shinyBS::tipify( shiny::actionLink(ns("gs_info"), "Tutorial", icon = shiny::icon("youtube")),
                   "Show more information about this module."),
            shiny::hr(), shiny::br(),             
            shinyBS::tipify( shiny::selectInput(ns("gs_contrast"),"Contrast:", choices=NULL),
                   "Select a contrast of interest for the analysis.", placement="top"),
            shinyBS::tipify( shiny::selectInput(ns("gs_features"),"Gene set collection:", choices=NULL, multiple=FALSE),
                   "Choose a specific gene set collection for the analysis.", placement="top"),
           
            shinyBS::tipify( shiny::selectInput(ns("gs_fdr"),"FDR", c(1e-9,1e-6,1e-3,0.01,0.05,0.1,0.2,0.5,1), selected=0.2),
                           "Set the false discovery rate (FDR) threshold.", placement="top"),
                    shinyBS::tipify( shiny::selectInput(ns("gs_lfc"),"logFC threshold",
                                        choices=c(0,0.05,0.1,0.2,0.5,1,2), selected=0),
                           "Set the logarithmic fold change (logFC) threshold.",
                           placement="top"),
          
            shinyBS::tipify(shiny::actionLink(ns("gs_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Toggle advanced options.", placement="top"),
            shiny::br(),br(),
            shiny::conditionalPanel(
                "input.gs_options % 2 == 1", ns=ns, 
                shiny::tagList(
                    shinyBS::tipify(shiny::checkboxInput(ns("gs_showall"),"Show all genesets",FALSE),
                           "Enbale significant genes filtering. Display only significant genesets in the table.", 
                           placement="top", options = list(container = "body")),
                    
                    shinyBS::tipify(shiny::checkboxGroupInput(ns('gs_statmethod'),'Statistical methods:', choices=NULL),
                           "Select a method or multiple methos for the statistical test.", placement="right", options = list(container="body"))
                )
            )
        )
    )
}

EnrichmentUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        flex = c(1.75,1),
        height = 800,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Top enriched",uiOutput(ns("topEnriched_UI"))),
            shiny::tabPanel("Plots",uiOutput(ns("subplots_UI"))),
            shiny::tabPanel("Compare",uiOutput(ns("compare_UI"))),
            shiny::tabPanel("Volcano (all)",uiOutput(ns("volcanoAll_UI"))),
            shiny::tabPanel("Volcano (methods)",uiOutput(ns("volcanoMethods_UI")))
        ),
        shiny::tabsetPanel(
            id = ns("tabs2"),
            shiny::tabPanel("Table",uiOutput(ns("tables_UI"))),
            shiny::tabPanel("Foldchange (all)",uiOutput(ns("fctable_UI"))),
            shiny::tabPanel("FDR table",uiOutput(ns("FDRtable_UI")))                       
        )
    )
}
