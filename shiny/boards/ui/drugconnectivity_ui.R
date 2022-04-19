##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

DrugConnectivityInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    bigdash::tabSettings(
        withTooltip( shiny::actionLink(ns("dsea_info"), "Youtube", icon = shiny::icon("youtube") ),
                "Show more information about this module."),
        shiny::hr(), shiny::br(),             
        withTooltip( shiny::selectInput(ns("dsea_contrast"),"Contrast:", choices=NULL),
                "Select the contrast corresponding to the comparison of interest.",
                placement="top"),
        withTooltip( shiny::selectInput(ns('dsea_method'),"Analysis type:", choices = ""),
                "Select type of drug enrichment analysis: activity or sensitivity (if available).",
                placement="top")
    )
}

DrugConnectivityUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        flex = c(1),
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("Drug enrichment",
                shiny::fillCol(
                flex = c(NA,0.035,1),
                height = 750,            
                tags$div(
                    HTML("<b>(a)</b> <b>Drug connectivity</b> correlates your signature with known drug perturbation
                     profiles from the L1000 database. The figures show the most similar (or opposite) profiles by running
                     the GSEA algorithm on the profile correlation space. <b>(b)</b> <b>Enrichment table</b> summarizing
                     the statistical results of the drug enrichment analysis. <b>(c)</b> <b>Mechanism-of-action</b>
                     plot showing the top most frequent drug class (or target genes) having similar or opposite enrichment
                     compared to the query signature. <b>(d)</b> <b>Activation matrix</b> visualizing enrichment
                     levels of drug signatures across multiple contrast profiles.")
                ),
                shiny::br(),
                shiny::fillRow(
                    height = 660,
                    flex = c(2.6,1), 
                    shiny::fillCol(
                        flex = c(1.5,0.15,1),
                        height = 660,
                        shiny::fillRow(
                            flex=c(1.2,0.04,1),
                            plotWidget(ns("dsea_enplots")),
                            shiny::br(),
                            plotWidget(ns("dsea_moaplot"))
                        ),
                        shiny::br(),  ## vertical space
                        tableWidget(ns("dsea_table"))        
                    ),
                    plotWidget(ns("dsea_actmap"))
                )
                )
            ),
            shiny::tabPanel("Connectivity map (beta)",
                shiny::fillCol(
                flex = c(NA,0.035,1),
                height = 750,            
                tags$div(
                    HTML("<b>(a)</b> <b>Enrichment plot.</b> Enrichment of the selected drug perturbation
                     profile with your signature. <b>(b)</b> <b>Enrichment table</b> summarizing the statistical
                     results of the drug enrichment analysis. <b>(c)</b> <b>Connectivity map.</b>
                     Plot showing the top signatures as UMAP. Each point is one L1000 experiment.
                     The color corresponds to the rank correlation between the drug signatures and your selected contrast.")
                ),
                shiny::br(),
                shiny::fillRow(
                    height = 660,
                    flex = c(1,0.05,1.5),
                    shiny::fillCol(
                        flex = c(1.15,0.05,1),                    
                        plotWidget(ns("cmap_enplot")),
                        shiny::br(),
                        tableWidget(ns("cmap_table"))                    
                    ),
                    shiny::br(),
                    plotWidget(ns("dsea_cmap"))
                )
            )
        )
        )
    )
}
