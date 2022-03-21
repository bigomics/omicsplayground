##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


ExpressionInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        tags$div(
         HTML("<b>Differential Expression Analysis.</b> Compare expression between
        two conditions. Determine which genes are significantly downregulated or overexpressed in one of the groups.")
        ),
        shiny::tagList(
            shinyBS::tipify( shiny::actionLink(ns("gx_info"), "Tutorial", icon = shiny::icon("youtube")),
                   "Show more information about this module."),
            shiny::hr(), shiny::br(),             
            shinyBS::tipify( shiny::selectInput(ns("gx_contrast"), "Contrast:", choices=NULL),
                   "Select a contrast of interest for the analysis.", placement="top"),
            shinyBS::tipify( shiny::selectInput(ns("gx_features"),"Gene family:", choices=NULL, multiple=FALSE),
                   "Choose a specific gene family for the analysis.", placement="top"),
            shiny::fillRow( flex=c(1,1),
                    shinyBS::tipify( shiny::selectInput(ns("gx_fdr"),"FDR", choices=c(1e-9,1e-6,1e-3,0.01,0.05,0.1,0.2,0.5,1), selected=0.2),
                           "Set the false discovery rate (FDR) threshold.", placement="top"),
                    shinyBS::tipify( shiny::selectInput(ns("gx_lfc"),"logFC threshold",
                                        choices=c(0,0.1,0.2,0.5,1,2,5), selected=0),
                           "Set the logarithmic fold change (logFC) threshold.", placement="top")
                    ),
            shiny::br(),br(),br(),br(),
            shinyBS::tipify( shiny::actionLink(ns("gx_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Toggle advanced options.", placement="top"),
            shiny::br(),br(),
            shiny::conditionalPanel(
                "input.gx_options % 2 == 1", ns=ns,
                shiny::tagList(
                    shinyBS::tipify(shiny::checkboxInput(ns("gx_showall"),"show all genes", FALSE),
                           "Display all genes in the table. Disable filtering of significant genes.", 
               placement="top", options = list(container = "body")),
                    shinyBS::tipify( shiny::checkboxGroupInput(ns('gx_statmethod'),'Statistical methods:',
                                               choices=NULL, inline=TRUE),
                           "Select a method for the statistical test. To increase the statistical reliability of the Omics Playground,
                            we perform the DE analysis using commonly accepted methods in the literature, including t-test (standard, Welch),
                            limma (no trend, trend, voom), edgeR (QLF, LRT), and DESeq2 (Wald, LRT), and merge the results.",
                             placement="right", options=list(container="body"))
                )
            )
        )
    )
}

ExpressionUI.test <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tabsetPanel(
        shiny::tabPanel("Table",uiOutput(ns("tables_UI"))),
        shiny::tabPanel("Foldchange (all)",uiOutput(ns("fctable_UI"))),
        shiny::tabPanel("FDR table",uiOutput(ns("FDRtable_UI")))                       
    )
}

ExpressionUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        flex = c(1.5,1),
        height = 780,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Plot",uiOutput(ns("plots_UI"))),
            shiny::tabPanel("Top genes",uiOutput(ns("topgenesUI"))),
            shiny::tabPanel("Volcano (all)",uiOutput(ns("volcanoAll_UI"))),
            ## shiny::tabPanel("Volcano (all2)",uiOutput(ns("volcanoAll2_UI"))),
            shiny::tabPanel("Volcano (methods)",uiOutput(ns("volcanoMethodsUI")))
        ),
        shiny::tabsetPanel(
            id = ns("tabs2"),
            shiny::tabPanel("Table",uiOutput(ns("tables_UI"))),
            shiny::tabPanel("Foldchange (all)",uiOutput(ns("fctable_UI"))),
            shiny::tabPanel("FDR table",uiOutput(ns("FDRtable_UI")))                       
        )
    )
}
