##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

CorrelationInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),

    ## data set parameters
    withTooltip(shiny::selectInput(ns("cor_gene"), "Gene:", choices = NULL),
      "Choose a gene for the correlation analysis.",
      placement = "top"
    ),
    shiny::br(),
    withTooltip(shiny::selectInput(ns("cor_features"), "Filter genes:", choices = NULL, multiple = FALSE),
      "Filter gene features.",
      placement = "top"
    ),
    shiny::br(),
    withTooltip(
      shiny::radioButtons(ns("pcor_ntop"), "Nr of top genes to compute partial correlation.",
        c(50, 100, 250),
        selected = 100, inline = TRUE
      ),
      "Top genes",
      placement = "top"
    ),
    shiny::conditionalPanel(
      "input.cor_features == '<custom>'",
      ns = ns,
      withTooltip(
        shiny::textAreaInput(ns("cor_customfeatures"),
          NULL,
          value = NULL,
          height = "100px", width = "100%",
          rows = 5, placeholder = "Paste your custom gene list"
        ),
        "Paste a custom list of genes to be used as features.",
        placement = "top"
      )
    )
  )
}

CorrelationUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 800 ## full height of page
  rowH <- 340 ## full height of page

  tabs <- shiny::tabsetPanel(
    id = ns("tabs"),
    shiny::tabPanel(
      "Correlation",
      div(
        class = "row",
        div(
          class = "col-md-6",
          correlation_plot_table_corr_ui(ns("cor_barplot"),
            label = "a",
            height = c("30vh", "70vh"),
            width = c("auto", "100%")
          ),
        ),
        div(
          class = "col-md-6",
          correlation_plot_scattercorr_ui(ns("cor_scatter"),
            height = c(fullH - 50, 700),
            width = c("auto", "100%")
          )
        )
      )
    ),
    shiny::tabPanel(
      "Graph",
      div(
        class = "row",
        div(
          class = "col-md-6",
          correlation_plot_cor_graph_ui(ns("cor_graph"))
        ),
        div(
          class = "col-md-6",
          correlation_plot_correlation_UMAP_ui(ns("cor_umap"))
        )
      )
    )
  )
  div(
    boardHeader(title = "Correlation analysis", info_link = ns("data_info")),
    tabs
  )
}
