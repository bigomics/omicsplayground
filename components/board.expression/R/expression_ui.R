##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


ExpressionInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    bigdash::tabSettings(
        withTooltip( shiny::actionLink(ns("gx_info"), "Tutorial", icon = shiny::icon("youtube")),
                "Show more information about this module."),
        shiny::hr(), shiny::br(),             
        withTooltip( shiny::selectInput(ns("gx_contrast"), "Contrast:", choices=NULL),
                "Select a contrast of interest for the analysis.", placement="top"),
        withTooltip( shiny::selectInput(ns("gx_features"),"Gene family:", choices=NULL, multiple=FALSE),
                "Choose a specific gene family for the analysis.", placement="top"),
        shiny::fillRow( flex=c(1,1),
                withTooltip( shiny::selectInput(ns("gx_fdr"),"FDR", choices=c(1e-9,1e-6,1e-3,0.01,0.05,0.1,0.2,0.5,1), selected=0.2),
                        "Set the false discovery rate (FDR) threshold.", placement="top"),
                withTooltip( shiny::selectInput(ns("gx_lfc"),"logFC",
                                    choices=c(0,0.1,0.2,0.5,1,2,5), selected=0),
                        "Set the logarithmic fold change (logFC) threshold.", placement="top")
                ),
        shiny::br(),br(),br(),br(),
        withTooltip( shiny::actionLink(ns("gx_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                "Toggle advanced options.", placement="top"),
        shiny::br(),br(),
        shiny::conditionalPanel(
            "input.gx_options % 2 == 1", ns=ns,
            shiny::tagList(
                withTooltip(shiny::checkboxInput(ns("gx_showall"),"show all genes", FALSE),
                        "Display all genes in the table. Disable filtering of significant genes.", 
            placement="top", options = list(container = "body")),
                withTooltip( shiny::checkboxGroupInput(ns('gx_statmethod'),'Statistical methods:',
                                            choices=NULL, inline=TRUE),
                        "Select a method for the statistical test. To increase the statistical reliability of the Omics Playground,
                        we perform the DE analysis using commonly accepted methods in the literature, including t-test (standard, Welch),
                        limma (no trend, trend, voom), edgeR (QLF, LRT), and DESeq2 (Wald, LRT), and merge the results.",
                            placement="right", options=list(container="body"))
            )
        )
    )
}

ExpressionUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        flex = c(1.5,1),
        height = 780,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Plot",
                shiny::fillCol(
                height = 390,
                flex = c(1,0.35,NA),
                    shiny::fillRow(
                        id = "plots",
                        flex=c(1,1,1,1),
                        plotWidget(ns("plots_volcano")),
                        plotWidget(ns("plots_maplot")),
                        plotWidget(ns("plots_boxplot")),
                        plotWidget(ns("plots_topfoldchange"))
                    ),
                    shiny::br(),
                    tags$div(
                        HTML("<b>Expression plots</b> associated with the selected contrast. <b>(a)</b> Volcano-plot plotting fold-change versuson
                            significance the x and y axes, respectively. <b>(b)</b> MA-plot plotting signal intensity versus fold-change on the x and y axes,
                            respectively. <b>(c)</b> Sorted barplot of the top diffentially expressed genes with largest (absolute) fold-change
                            for selected contrast. <b>(d)</b> Sorted barplot of the differential expression of the selected gene across all contrasts.")
                    )
                )
            ),
            shiny::tabPanel("Top genes",
                shiny::fillCol(
                height = 390,
                flex=c(1,NA,NA),
                plotWidget(ns("topgenes")),
                shiny::br(),
                tags$div(
                        HTML("<b>Top differentially expressed genes.</b> Expression barplots of the top most differentially
                         (both positively and negatively) expressed genes for the selected contrast.")
                )
            )),
            shiny::tabPanel("Volcano (all)",
                shiny::fillCol(
                height = 390,
                flex=c(1,NA,NA),
                plotWidget(ns("volcanoAll")),
                shiny::br(),
                tags$div(
                        HTML("<b>Volcano plot for all contrasts.</b> Simultaneous visualisation of volcano
                         plots of genes for all contrasts. Experimental contrasts with better statistical significance will
                          show volcano plots with 'higher' wings.")
                )
            )),
           
            shiny::tabPanel("Volcano (methods)",
                shiny::fillCol(
                height = 390,
                flex=c(1,NA,NA),
                plotWidget(ns("volcanoMethods")),
                shiny::br(),
                tags$div(
                        HTML("<b>Volcano plot for all statistical methods.</b> Simultaneous visualisation of volcano plots
                         of genes by multiple differential expression methods for the selected contrast.
                          Methods showing better statistical significance will show volcano plots with 'higher' wings.")
                )
            ))
        ),
        shiny::tabsetPanel(
            id = ns("tabs2"),
            shiny::tabPanel("Table",
                shiny::fillCol(
                    height = 1.15*320,
                    flex = c(NA,0.06,1),
                    tags$div(
                        HTML("<b>Differential Expression Analysis.</b> Compare expression between
                        two conditions. Determine which genes are significantly downregulated or overexpressed in one of the groups.")
                    ),
                    shiny::br(),
                    shiny::fillRow(
                        flex = c(1.6,0.07,1), 
                        tableWidget(ns("genetable")),
                        shiny::br(),
                        tableWidget(ns("gsettable"))
                    )
                )
            ),
            shiny::tabPanel("Foldchange (all)",
                shiny::fillCol(
                height = 320,
                tableWidget(ns("fctable"))
            )),
            shiny::tabPanel("FDR table",
                shiny::fillCol(
                height = 320,
                tableWidget(ns("FDRtable"))
            ))                       
        )
    )
}
