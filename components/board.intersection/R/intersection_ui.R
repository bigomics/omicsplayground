##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

IntersectionInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    bigdash::tabSettings(
        withTooltip( shiny::actionLink(ns("info"), "Tutorial", icon = shiny::icon("youtube")),
                "Show more information about this module"),
        shiny::hr(), shiny::br(),             
        withTooltip( shiny::selectInput(ns('comparisons'),'Contrasts:', choices=NULL, multiple=TRUE),
                "Select the contrasts that you want to compare. If you select N=2 contrast a single scatterplot will be drawn. For N>=3 a scatterplot matrix will be drawn.",
                placement="top"),

        withTooltip( shiny::actionLink(ns("options"), "Options", icon=icon("cog", lib = "glyphicon")),
                "Toggle advanced options.", placement="top"),
        shiny::br(),br(),
        shiny::conditionalPanel(
            "input.options % 2 == 1", ns=ns,
            withTooltip( shiny::radioButtons(ns("level"),"Level:",
                                    choices=c("gene","geneset"), inline=TRUE),
                    "Select feature level: gene or geneset", placement="top"),                
            withTooltip( shiny::selectInput(ns("filter"),"Filter:", choices=NULL, multiple=FALSE),
                    "Filter features", placement="top"),
            shiny::conditionalPanel(
                "input.filter == '<custom>'", ns=ns,
                withTooltip( shiny::textAreaInput(ns("customlist"), NULL, value = NULL,
                                        rows=5, placeholder="Paste your custom gene list"),
                        "Paste a custom list of genes to highlight.", placement="bottom")
            )
        )
    )
}

IntersectionUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace

    shiny::tabsetPanel(
        id = ns("tabs1"),
        shiny::tabPanel(
            "Pairwise scatter",
            tags$div(
                HTML("<h4>Pairwise scatter & Venn diagram</h4> <b>(a)</b> <b>Pairs plot.</b> Pairwise scatterplots
                for two or more differential expression profiles for multiple selected contrasts.
                Similar profiles will show high correlation with points close to the diagonal. 
                <b>(b)</b> <b>Venn diagram</b> showing the number of overlapping genes for multiple contrasts.
                <b>(c)</b> <b>Venn table.</b> Genes in the selected overlap region.")
            ),
            div(
                class = "row",
                div(
                    class = "col-md-6",
                    plotWidget(ns("scatterPlotMatrix"))
                ),
                div(
                    class = "col-md-6",
                    plotWidget(ns("venndiagram")),
                    tableWidget(ns("venntable"))
                )
            )
        ),
        shiny::tabPanel(
            "Signature clustering",
            tags$div(
                HTML("<h4>Signature clustering</h4> <b>(a)</b> <b>Signature heatmap.</b> Similarity of the
                signatures visualized as a clustered heatmap. The top plot shows the distribution of foldchange
                values as boxplots. <b>(b)</b> <b>Contrast correlation.</b> The numeric values in the cells
                correspond to the Pearson correlation coefficient. Red corresponds to positive correlation
                    and blue to negative correlation.")
            ),
            div(
                class = "row",
                div(
                    class = "col-md-6",
                    plotWidget(ns("FoldchangeHeatmap"))
                ),
                div(
                    class = "col-md-6",
                    plotWidget(ns("ctcorrplot"))
                )
            )
        )
    )
}
