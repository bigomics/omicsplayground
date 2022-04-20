##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

FeatureMapInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    bigdash::tabSettings(
        shiny::actionLink(ns("info"), "Info", icon=icon("info-circle")),
        shiny::hr(), shiny::br(),
        ## data set parameters
        shiny::selectInput(ns('sigvar'),'Show phenotype:', choices=NULL, multiple=FALSE),
        shiny::br(),
        shiny::br(),
        shiny::actionLink(ns("options"), "Options", icon=icon("cog", lib = "glyphicon")),
        shiny::br(),br(),
        shiny::conditionalPanel(
            "input.options % 2 == 1", ns=ns,
            shiny::tagList(
                tipifyR(shiny::radioButtons(ns('umap_type'),'UMAP datatype:',
                                        choices=c('logCPM','logFC'), inline=TRUE),
                        "The UMAP can be computed from the normalized log-expression (logCPM), or from the log-foldchange matrix (logFC). Clustering based on logCPM is the default, but when batch/tissue effects are present the logFC might be better."),
                tipifyR(shiny::selectInput(ns('filter_genes'),'Show genes:',
                                    choices=NULL, multiple=FALSE),
                        "Filter the genes to highlight on the map."),
                tipifyR(shiny::selectInput(ns('filter_gsets'),'Show genesets:',
                                    choices=NULL, multiple=FALSE),
                        "Filter the genesets to highlight on the map.")
            )
        )
    )
}

FeatureMapUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    ui <- shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("Gene", 
                shiny::fillCol(
                flex = c(NA,0.02,1,0.3),
                height = 1.1*800,
                tags$div(
                    HTML("<h4>Gene UMAP</h4>")
                ),
                shiny::br(),
                shiny::fillRow(
                    flex = c(1,0.03,1.2),
                    height = 0.85*800,                
                    plotWidget(ns('geneUMAP')),
                    shiny::br(),                
                    plotWidget(ns('geneSigPlots'))
                ),
                tableWidget(ns('geneTable'))
            )),
            shiny::tabPanel("Geneset",
                shiny::fillCol(
                flex = c(NA,0.02,1,0.3),
                height = 1.1*800,
                tags$div(
                    HTML("<h4>Geneset UMAP</h4>")
                ),
                shiny::br(),
                shiny::fillRow(
                    flex = c(1,0.03,1.2),
                    plotWidget(ns('gsetUMAP')),
                    shiny::br(),
                    plotWidget(ns('gsetSigPlots'))
                ),
                tableWidget(ns('gsetTable'))
            ))            
        )
    )
    ui
}