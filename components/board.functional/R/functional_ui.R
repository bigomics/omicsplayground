##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

FunctionalInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    bigdash::tabSettings(
        withTooltip( shiny::actionLink(ns("fa_info"), "Youtube", icon = shiny::icon("youtube") ),
                "Show more information about this module."),
        shiny::hr(), shiny::br(),             
        withTooltip( shiny::selectInput(ns("fa_contrast"),"Contrast:", choices=NULL),
                "Select the contrast corresponding to the comparison of interest.",
                placement="top"),
        withTooltip( shiny::actionLink(ns("fa_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                "Show/hide advanced options", placement="top"),
        shiny::br(),
        shiny::conditionalPanel(
            "input.fa_options % 2 == 1", ns=ns,
            shiny::tagList(
                withTooltip(shiny::checkboxInput(ns('fa_filtertable'),'filter signficant (tables)',FALSE),
                        "Click to filter the significant entries in the tables.")
            )
        )
    )
}

FunctionalUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        flex = c(1),
        height = 780,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("KEGG",
                shiny::fillCol(
                height = 750,
                flex = c(NA,0.035,1),
                tags$div(
                    HTML("<b>(a)</b> <b>KEGG pathway map.</b> Genes are colored according to their upregulation (red) or downregulation (blue) in the contrast profile. <b>(b)</b> <b>Enrichment table</b> reporting enrichment score for each pathway for the selected contrast profile. <b>(c)</b> <b>Activation matrix</b> visualizing the activation levels of pathways across contrasts.")
                ),
                shiny::br(),
                shiny::fillRow(
                    flex = c(1.3,0.1,1),
                    height = 660,
                    shiny::fillCol(
                        flex = c(1.8,0.1,1),
                        height = 0.9*660,
                        plotWidget(ns("kegg_graph")),
                        shiny::br(),
                        tableWidget(ns("kegg_table"))
                    ),
                    shiny::br(),
                    plotWidget(ns("kegg_actmap"))
                )
            )),
            shiny::tabPanel("GO graph",
            shiny::fillCol(
                flex=c(NA,0.035,1),
                height = 750,
                tags$div(
                    HTML("<b>(a)</b> <b>Gene Ontology graph.</b> The graph represents the enrichment of the GO
                     terms as a tree structure. <b>(b)</b> <b>GO score table.</b> The score of a GO term is the
                      cumulative score of all higher order terms. <b>(c)</b> <b>Activation matrix</b>
                       visualizing the enrichment of GO terms across multiple contrast profiles.")
                ),
                shiny::br(),
                shiny::fillRow(
                    height = 660,
                    flex = c(1.2,1),
                    shiny::fillCol(
                        flex = c(2,1),
                        height = 0.9*660,
                        plotWidget(ns("GO_network")),
                        tableWidget(ns("GO_table"))
                    ),
                    plotWidget(ns("GO_actmap"))
                )
            ))
        )
    )
}
