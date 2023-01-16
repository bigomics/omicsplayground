##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

DrugConnectivityInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    withTooltip(
      shiny::actionLink(ns("dsea_info"), "Youtube", icon = shiny::icon("youtube")),
      "Show more information about this module."
    ),
    shiny::hr(), shiny::br(),
    withTooltip(shiny::selectInput(ns("dsea_contrast"), "Contrast:", choices = NULL),
      "Select the contrast corresponding to the comparison of interest.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("dsea_method"), "Analysis type:", choices = ""),
      "Select type of drug enrichment analysis: activity or sensitivity (if available).",
      placement = "top"
    )
  )
}

DrugConnectivityUI <- function(id) {
  ns <- shiny::NS(id)

  div(
    boardHeader(title = "Drug Connectivity", info_link = ns("dsea_info")),
    div(
      shiny::tabsetPanel(
        id = ns("tabs"),
        shiny::tabPanel(
          "Drug enrichment",
          div(class = "row",
            div(class = "col-md-10",
              div(class = "row",
                div(class = "col-md-6",
                  plotWidget(ns("dsea_enplots")),
                ),
                div(class = "col-md-6",
                  #plotWidget(ns("dsea_moaplot"))
                  drugconnectivity_plot_moa_ui(ns("dsea_moaplot"))
                )
              ),
              br(),
              tableWidget(ns("dsea_table"))
            ),
            div(class = "col-md-2",
              plotWidget(ns("dsea_actmap"))
            )
          ),
          div(
            HTML("<b>(a)</b> <b>Drug connectivity</b> correlates your signature with known drug perturbation
                  profiles from the L1000 database. The figures show the most similar (or opposite) profiles by running
                  the GSEA algorithm on the profile correlation space. <b>(b)</b> <b>Enrichment table</b> summarizing
                  the statistical results of the drug enrichment analysis. <b>(c)</b> <b>Mechanism-of-action</b>
                  plot showing the top most frequent drug class (or target genes) having similar or opposite enrichment
                  compared to the query signature. <b>(d)</b> <b>Activation matrix</b> visualizing enrichment
                  levels of drug signatures across multiple contrast profiles.")
          )
        ),
        shiny::tabPanel(
          "Connectivity map (beta)",
          shiny::fillCol(
            flex = c(NA, 0.035, 1),
            height = 750,
            shiny::fillRow(
              height = 660,
              flex = c(1, 0.05, 1.5),
              shiny::fillCol(
                flex = c(1.15, 0.05, 1),
                plotWidget(ns("cmap_enplot")),
                shiny::br(),
                tableWidget(ns("cmap_table"))
              ),
              shiny::br(),
              plotWidget(ns("dsea_cmap"))
            ),
            div(
              HTML("<b>(a)</b> <b>Enrichment plot.</b> Enrichment of the selected drug perturbation
                     profile with your signature. <b>(b)</b> <b>Enrichment table</b> summarizing the statistical
                     results of the drug enrichment analysis. <b>(c)</b> <b>Connectivity map.</b>
                     Plot showing the top signatures as UMAP. Each point is one L1000 experiment.
                     The color corresponds to the rank correlation between the drug signatures and your selected contrast.")
            )
          )
        )
      )
    )
  )
}
