##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

EnrichmentInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    bigdash::tabSettings(
        withTooltip( shiny::actionLink(ns("gs_info"), "Tutorial", icon = shiny::icon("youtube")),
                "Show more information about this module."),
        shiny::hr(), shiny::br(),             
        withTooltip( shiny::selectInput(ns("gs_contrast"),"Contrast:", choices=NULL),
                "Select a contrast of interest for the analysis.", placement="top"),
        withTooltip( shiny::selectInput(ns("gs_features"),"Gene set collection:", choices=NULL, multiple=FALSE),
                "Choose a specific gene set collection for the analysis.", placement="top"),
        
        withTooltip( shiny::selectInput(ns("gs_fdr"),"FDR", c(1e-9,1e-6,1e-3,0.01,0.05,0.1,0.2,0.5,1), selected=0.2),
                        "Set the false discovery rate (FDR) threshold.", placement="top"),
                withTooltip( shiny::selectInput(ns("gs_lfc"),"logFC threshold",
                                    choices=c(0,0.05,0.1,0.2,0.5,1,2), selected=0),
                        "Set the logarithmic fold change (logFC) threshold.",
                        placement="top"),
        
        withTooltip(shiny::actionLink(ns("gs_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                "Toggle advanced options.", placement="top"),
        shiny::br(),br(),
        shiny::conditionalPanel(
            "input.gs_options % 2 == 1", ns=ns, 
            shiny::tagList(
                withTooltip(shiny::checkboxInput(ns("gs_showall"),"Show all genesets",FALSE),
                        "Enbale significant genes filtering. Display only significant genesets in the table.", 
                        placement="top", options = list(container = "body")),
                
                withTooltip(shiny::checkboxGroupInput(ns('gs_statmethod'),'Statistical methods:', choices=NULL),
                        "Select a method or multiple methos for the statistical test.", placement="right", options = list(container="body"))
            )
        )
    )
}

EnrichmentUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace

    tagList(
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Top enriched",
                div(
                    class = "row",
                    div(
                        class = "col-md-6",
                        plotWidget(ns("topEnriched"))
                    ),
                    div(
                        class = "col-md-6",
                        plotWidget(ns("topEnrichedFreq"))
                    )
                ),
                div(
                    HTML(
                        "<b>(a)</b> <b>Top enriched gene sets.</b> Enrichment plots of the top differentially enriched gene sets. Black vertical bars indicate the rank of genes in the gene set in the sorted list metric. The green curve corresponds to the 'running statistics' of the enrichment score.",
                        "<b>(b)</b> <strong>Gene frequency.</strong> The plot shows the number of times a gene is present in the top-N genesets sorted by frequency."
                    )
                )
            ),
            shiny::tabPanel("Plots",
                div(
                    class = "row",
                    div(
                        class = "col-md-3",
                        plotWidget(ns("subplot_volcano"))
                    ),
                    div(
                        class = "col-md-3",
                        plotWidget(ns("subplot_barplot"))
                    ),
                    div(
                        class = "col-md-3",
                        plotWidget(ns("subplot_geneplot"))
                    ),
                    div(
                        class = "col-md-3",
                        plotWidget(ns("subplot_scatter"))
                    )
                ),
                shiny::br(),
                tags$div(
                    HTML("<b>Enrichment plots</b> associated with the gene set (selected in <b>Table I</b>) and gene (selected in <b>Table II</b>).
                    <b>(a)</b> Volcano-plot showing significance versus fold-change on the y and x axes, respectively.
                    Genes in the gene set are highlighted in blue. <b>(b)</b> Barplot of the gene set enrichment in the groups.
                    <b>(c)</b> Barplot of selected gene in the groups. <b>(d)</b> Scatter plot of the enrichment versus the expression
                    of the selected geneset and gene, on the y and x axes, respectively.")
                )
            ),
            shiny::tabPanel("Compare",
                plotWidget( ns("compare")),
                tags$div(
                    HTML("<b>Enrichment across contrasts.</b> Enrichment plots for the selected gene set (in <b>Table I</b>)
                    across multiple contrasts. The figure allows to quickly compare the enrichment of a certain gene set
                    across all other comparisons.")
                )
            ),
            shiny::tabPanel("Volcano (all)",
                shiny::fillCol(
                height = 420,
                flex=c(1,NA),
                plotWidget(ns("volcanoAll")),
                tags$div(
                    HTML("<b>Volcano plots for all contrasts.</b> Simultaneous visualisation of volcano plots of gene
                    set enrichment across all contrasts. Volcano-plot are plotting enrichment score versus
                    significance on the x and y axes, respectively. Experimental contrasts showing better statistical
                    significance will show volcano plots with 'higher' wings.")
                )
            )),
            shiny::tabPanel("Volcano (methods)",
                plotWidget(ns("volcanoMethods")),
                tags$div(
                    HTML("<b>Volcano plots for all methods.</b> Simultaneous visualisation of volcano plots of gene
                    sets for different enrichment methods. Methods showing better statistical
                    significance will show volcano plots with 'higher' wings.")
                )
            )
        ),
        shiny::tabsetPanel(
            id = ns("tabs2"),
            shiny::tabPanel("Table",
                tags$div(
                    HTML("<b>Enrichment tables</b>. <b>(I)</b> Table summarizing the statistical results of the gene
                    set enrichment analysis for selected contrast. The number of stars indicate how many methods
                    identified the geneset significant. <b>(II)</b> Table showing the fold-change,
                    statistics and correlation of the genes in the selected gene set.")
                ),
                div(
                    class = "row",
                    div(
                        class = "col-md-7",
                        tableWidget(ns("gseatable"))
                    ),
                    div(
                        class = "col-md-5",
                        tableWidget(ns("genetable"))        
                    )
                )
            ),
            shiny::tabPanel("Foldchange (all)",
                plotWidget(ns("fctable"))
            ),
            shiny::tabPanel("FDR table",
                tableWidget(ns("FDRtable"))
            )                    
        )
    )
}
