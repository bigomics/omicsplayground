##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
          wordcloud_table_enrichment_ui(
            ns("wordcloud_enrichmentTable"),
            height = c("35vh", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          )
        ),
        div(
          class = "col-md-6",
          wordcloud_table_leading_edge_ui(
            ns("wordcloud_leadingEdgeTable"),
            height = c("35vh", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          )
        )
      )
    )
  )
  div(
    boardHeader(title = "Word cloud", info_link = ns("wc_info")),
    tabs
  )
}
