##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
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
            height = c(0.45 * fullH, 700),
            width = c("auto", 1200)
          ),
        ),
        div(
          class = "col-md-6",
          correlation_plot_scattercorr_ui(ns("cor_scatter"),
            height = c(fullH - 50, 760),
            width = c("auto", 900)
          )
        )
      ),
      tags$div(
        HTML("<b>(a)</b>
                <b>Top-ranked correlation.</b> Top correlated features with respect to selected gene.
                <b>(b)</b> <b>Correlation table</b> of correlation and partial
                correlation with respect to selected gene. <b>(c)</b> <b>Scatter plots</b> of gene
                expression of top correlated genes.")
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
          correlation_plot_partial_correlation_ui(ns("cor_umap"))
        )
      ),
      div(
        HTML("Visualization of gene correlation as network or UMAP. <b>
            (a)</b> <b>Partial correlation network</b> around the selected gene. <b>(b)</b>
            <b>Correlation UMAP</b>. Clustering of genes  colored by correlation (or covariance).")
      )
    )
  )
  div(
    boardHeader(title = "Correlation analysis", info_link = ns("data_info")),
    tabs
  )
}
