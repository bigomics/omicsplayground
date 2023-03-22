##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

DrugConnectivityInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(),
    withTooltip(shiny::selectInput(ns("dsea_contrast"), "Contrast:", choices = NULL),
      "Select the contrast corresponding to the comparison of interest.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("dsea_method"), "Analysis type:", choices = ""),
      "Select type of drug enrichment analysis: activity or sensitivity (if available).",
      placement = "top"
    ),
    shiny::hr(),
    withTooltip(
      shiny::checkboxInput(
        ns("dseatable_filter"),
        "only annotated drugs",
        FALSE
      ),
      "Show only annotated drugs."
    )
  )
}

DrugConnectivityUI <- function(id) {
  ns <- shiny::NS(id)

  div(
    boardHeader(title = "Drug Connectivity", info_link = ns("dsea_info")),
    div(selector_default(ns("hide_caption"), label = "Show captions")),
    shiny::tabsetPanel(
      id = ns("tabs"),
      shiny::tabPanel(
        "Drug enrichment",
        div(
          class = "row",
          div(
            class = "col-md-9",
            div(
              class = "row",
              div(
                class = "col-md-6",
                drugconnectivity_plot_enplots_ui(ns("dsea_enplots"), label = "a")
              ),
              div(
                class = "col-md-6",
                drugconnectivity_plot_moa_ui(ns("dsea_moaplot"), label = "c")
              )
            ),
            br(),
            drugconnectivity_table_dsea_ui(
              ns("dsea_table"),
              height = c(360, TABLE_HEIGHT_MODAL),
              width = c("100%", "90%")
            )
          ),
          div(
            class = "col-md-3",
            drugconnectivity_plot_actmap_ui(ns("dsea_actmap"), label = "d")
          )
        )
      ),
      shiny::tabPanel(
        "Connectivity map (beta)",
        shiny::div(
          shiny::fillCol(
            flex = c(NA, 0.035, 1),
            height = 750,
            shiny::fillRow(
              height = 660,
              flex = c(1, 0.05, 1.5),
              shiny::fillCol(
                flex = c(1.15, 0.05, 1),
                drugconnectivity_plot_cmap_enplot_ui(ns("cmap_enplot"), label = "a"),
                shiny::br(),
                drugconnectivity_table_cmap_ui(
                  ns("cmap_table"),
                  height = c(380, TABLE_HEIGHT_MODAL),
                  width = c("100%", "90%")
                )
              ),
              shiny::br(),
              drugconnectivity_plot_cmap_dsea_ui(ns("cmap_dsea"), label = "c")
            ),
            div(
              class = "caption",
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
