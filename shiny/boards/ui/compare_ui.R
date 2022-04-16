##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

CompareInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    bigdash::tabSettings(
        withTooltip( shiny::actionLink(ns("info"), "Info", icon = shiny::icon("info-circle")),
                "Show more information about this module"),
        shiny::hr(), shiny::br(),             
        withTooltip( shiny::selectInput(ns('contrast1'),'Dataset1:',
                            choices=NULL, multiple=TRUE),
                "Select the contrast that you want to compare.",
                placement="right", options = list(container = "body")
                ),
        shiny::br(),            
        withTooltip( shiny::selectInput(ns('dataset2'),"Dataset2:", choices=NULL),
                "Select second dataset to compare.",
                placement="right", options = list(container = "body")),
        withTooltip( shiny::selectInput(ns('contrast2'),NULL, choices=NULL, multiple=TRUE),
                "Select second contrast to compare.",
                placement="right", options = list(container = "body")),
        shiny::br(),
        withTooltip( shiny::actionLink(ns("options"), "Options", icon=icon("cog", lib="glyphicon")),
                "Toggle advanced options.",
                placement="right", options = list(container = "body")),
        shiny::br(),
        shiny::conditionalPanel(
            "input.options % 2 == 1", ns=ns,
            shiny::br(),
            withTooltip( shiny::radioButtons(ns('plottype'),"Plot type:",
                                    choices=c("volcano","MA","scatter","UMAP1","UMAP2","heatmap"),
                                    selected='UMAP1', inline=TRUE),
                    "Select plot type.",
                    placement="right", options = list(container = "body")),
            shiny::br(),
            withTooltip( shiny::radioButtons(ns('hilighttype'),"Highlight genes:",
                                    choices=c("top scoring","custom"),
                                    inline=TRUE),
                    "Select highlight type.",
                    placement="right", options = list(container = "body")),
            shiny::conditionalPanel(
                "input.hilighttype == 'custom'", ns=ns,
                withTooltip( shiny::textAreaInput(ns("genelist"),NULL, value = NULL,
                                        height = "100px", width = "100%", 
                                        rows=5, placeholder="Paste your custom gene list"),
                        "Paste a custom list of genes to highlight.",
                        placement="right")
            ),
            shiny::br(),            
            withTooltip(
                shiny::radioButtons( ns('ntop'),"ntop", choices=c(10,20,40,100),
                                selected=20, inline=TRUE),
                "number of top genes to show",
                placement="right", options = list(container = "body"))
        )
    )
}

CompareUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Compare",
                tags$div(
                        HTML("<h4>Compare Expression</h4>"
                        )
                ),
                shiny::fillRow(
                    flex = c(1,0.1,1),
                    plotWidget(ns("scatter1")),
                    shiny::br(),
                    plotWidget(ns("scatter2"))
                )
            ),
            shiny::tabPanel("Foldchange",
                tags$div(
                        HTML("<h4>Compare Foldchange</h4>"
                        )
                ),
                shiny::fillRow(
                    flex = c(1,0.05,0.9),
                    plotWidget(ns("fcfcplot")),
                    shiny::br(),
                    plotWidget(ns("cumfcplot"))
                ) 
            ),
            shiny::tabPanel("Gene Correlation", 
                shiny::fillCol(
                    height = 770,
                    flex = c(NA,0.015,1),
                    tags$div(
                        HTML("<h4>Compare Correlation</h4>"
                        )
                    ),
                    shiny::br(),
                    shiny::fillRow(
                        flex = c(1.0,0.05,1),
                        shiny::fillCol(
                            flex = c(2.2,0.01,1),
                            plotWidget(ns("multibarplot")),                    
                            shiny::br(),
                            tableWidget(ns("score_table"))
                        ),
                        shiny::br(),
                        plotWidget(ns("genecorr"))
                    )
                )
            )            
        )
    )
}