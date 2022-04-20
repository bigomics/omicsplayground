##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

CorrelationInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    bigdash::tabSettings(
        shiny::actionLink(ns("cor_info"), "Info", icon=icon("info-circle")),
        shiny::hr(), shiny::br(),             

        ## data set parameters
        withTooltip( shiny::selectInput(ns("cor_gene"),"Gene:", choices=NULL),
                "Choose a gene for the correlation analysis.", placement="top"),
        shiny::br(),
        withTooltip( shiny::selectInput(ns("cor_features"),"Filter genes:", choices=NULL, multiple=FALSE),
                        "Filter gene features.", placement="top"),
        shiny::conditionalPanel(
                    "input.cor_features == '<custom>'", ns=ns,
                    withTooltip( shiny::textAreaInput(ns("cor_customfeatures"),
                                                            NULL, value = NULL,
                                                            height = "100px", width = "100%", 
                                                            rows=5, placeholder="Paste your custom gene list"),
                                    "Paste a custom list of genes to be used as features.",
                                    placement="top")
                )
    )
}

CorrelationUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    ui <- shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("Correlation",
            shiny::fillCol(
            flex = c(NA,0.035,1),
            height = 800,
            tags$div(
                HTML("<h3>Gene Correlation Analysis</h3><b>(a)</b>
                 <b>Top-ranked correlation.</b> Top correlated features with respect to selected gene.
                <b>(b)</b> <b>Correlation table</b> of correlation and partial
                correlation with respect to selected gene. <b>(c)</b> <b>Scatter plots</b> of gene
                expression of top correlated genes.")
            ),
            shiny::br(),
            shiny::fillRow(
                flex = c(1,0.01,1.2),
                shiny::fillCol(
                    flex = c(1,1),
                    height = 720,
                    plotWidget(ns('cor_barplot')),
                    tableWidget(ns('cor_table'))
                ),	
                shiny::br(), ## spacer
                shiny::fillCol(
                    flex = c(NA,1),
                    height = 720,
                    plotWidget(ns('cor_scatter'))
                )
            )
        )
        ),
        shiny::tabPanel("Graph",
            shiny::fillCol(
            flex = c(NA,0.035,1),
            height = 800,
            tags$div(
                HTML("<h3>Gene Correlation Network</h3>Visualization of gene correlation as network or UMAP. <b>
                (a)</b> <b>Partial correlation network</b> around the selected gene. <b>(b)</b>
                 <b>Correlation UMAP</b>. Clustering of genes  colored by correlation (or covariance).")
            ),
            shiny::br(),
            shiny::fillRow(
                flex = c(1,0.05,1),
                plotWidget(ns('cor_graph')),
                shiny::br(),
                plotWidget(ns('cor_umap'))
            )
            )
        )
        )
    )
    ui
}
