##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

WordCloudInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(shiny::selectInput(ns("wc_contrast"), "Contrast:", choices = NULL),
      "Select the contrast corresponding to the comparison of interest.",
      placement = "top"
    )
  )
}

WordCloudUI <- function(id) {
  fullH <- 750
  rowH <- 660 ## row height of panel
  tabH <- 200 ## row height of panel
  tabH <- "70vh" ## row height of panel
  ns <- shiny::NS(id) ## namespace
  shiny::tabsetPanel(
    id = ns("tabs"),
    tabs <- shiny::tabPanel(
      "",
      div(
        class = "row",
        div(
          class = "col-md-4",
          wordcloud_plot_enrichment_ui(ns("gseaplots"), 0.5 * rowH)
        ),
        div(
          class = "col-md-4",
          wordcloud_plot_wordcloud_ui(ns("wordcloud"), 0.5 * rowH)
        ),
        div(
          class = "col-md-4",
          wordcloud_plot_wordtsne_ui(ns("wordtsne"), 0.5 * rowH)
        )
      ),
      div(
        class = "row",
        div(
          class = "col-md-6",
          wordcloud_table_enrichment_ui(ns("wordcloud_enrichmentTable"))
        ),
        div(
          class = "col-md-6",
          wordcloud_table_leading_edge_ui(ns("wordcloud_leadingEdgeTable"))
        )
      ),
      tags$div(
        HTML("<b>(a)</b> <b>Word enrichment</b>  plots for the top most significant contrasts. Black vertical bars indicate
                    the position of gene sets, in the ranked enrichment scores, that contains the *keyword*.
                    The green curve corresponds to 'running statistics' of the keyword enrichment score.
                    <b>(b)</b> <b>Word cloud.</b> The size of the words are relative to the normalized enrichment score
                    (NES) from the GSEA computation. <b>(c)</b> <b>Word t-SNE</b> of keywords extracted from the titles/descriptions
                    of the genesets. <b>(d)</b> <b>Enrichment table</b> of keywords for selected contrast. <b>(e)</b> <b>Leading edge terms</b>
                    for selected keyword.")
      )
    )
  )
  div(
    boardHeader(title = "Word cloud", info_link = ns("wc_info")),
    tabs
  )
}
