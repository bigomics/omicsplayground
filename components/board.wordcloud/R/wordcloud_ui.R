##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

WordCloudInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    bigdash::tabSettings(
        withTooltip( shiny::actionLink(ns("wc_info"), "Youtube", icon = shiny::icon("youtube") ),
                "Show more information about this module."),
        shiny::hr(), shiny::br(),             
        withTooltip( shiny::selectInput(ns("wc_contrast"),"Contrast:", choices=NULL),
                "Select the contrast corresponding to the comparison of interest.",
                placement="top"),
        withTooltip( shiny::actionLink(ns("wc_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                "Show/hide advanced options", placement="top"),
        shiny::br(),
        shiny::conditionalPanel(
            "input.wc_options % 2 == 1", ns=ns,
            shiny::tagList(
                withTooltip(shiny::checkboxInput(ns('wc_normalize'),'normalize activation matrix',TRUE),
                        "Click to 'normalize' the coloring of an activation matrices.")
            )
        )
    )
}

WordCloudUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        flex = c(1),
        height = 780,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("WordCloud",
            shiny::fillCol(
            height = 750,
            flex = c(NA,0.035,1),
            tags$div(
            HTML("<b>(a)</b> <b>Word enrichment</b>  plots for the top most significant contrasts. Black vertical bars indicate
            the position of gene sets, in the ranked enrichment scores, that contains the *keyword*.
            The green curve corresponds to 'running statistics' of the keyword enrichment score.
            <b>(b)</b> <b>Word cloud.</b> The size of the words are relative to the normalized enrichment score
            (NES) from the GSEA computation. <b>(c)</b> <b>Word t-SNE</b> of keywords extracted from the titles/descriptions
            of the genesets. <b>(d)</b> <b>Activation matrix</b> showing keyword enrichment across contrasts.
            <b>(e)</b> <b>Enrichment table</b> of keywords for selected contrast. <b>(f)</b> <b>Leading edge terms</b>
            for selected keyword.")
            ),
            shiny::br(),
            shiny::fillRow(
                height = 660,
                flex = c(3.8,1),
                shiny::fillCol(
                    flex=c(1.2,0.1,1),
                    height = 660,
                    shiny::fillRow(
                        flex = c(1.2,0.05,1,0.05,1),
                        plotWidget(ns("gseaplots")),
                        shiny::br(),
                        plotWidget(ns("wordcloud")),
                        shiny::br(),
                        plotWidget(ns("wordtsne"))
                    ),
                    shiny::br(),
                    shiny::fillRow(
                        flex=c(1,0.08,1),
                        tableWidget(ns("enrichmentTable")),
                        shiny::br(),
                        tableWidget(ns("leadingEdgeTable"))
                    )
                ),
                plotWidget(ns("actmap"))
            )
        ))
        )
    )
}
