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
          correlation_plot_corr_ui(
            id = ns("cor_barplot"),
            title = "Top correlated genes",
            info.text = "Highest correlated genes in respect to the selected gene. The height of the bars correspond to the Pearson correlation value. The dark grey bars correspond to the 'partial correlation' which essentially corrects the correlation value for indirect effects and tries to estimate the amount of direct interaction.",
            caption = "Barplot showing the highest correlated genes with respect to the selected gene.",
            label = "",
            height = c("calc(65vh - 217px)", "70vh"),
            width = c("auto", "100%")
          ),
          correlation_table_corr_ui(
            id = ns("cor_barplot"),
            title = "Correlation table",
            info.text = "Statistical results from correlated gene pairs.",
            caption = "Correlation table of correlation and partial correlation with respect to the selected gene. ",
            label = "",
            height = c("35vh", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        ),
        div(
          class = "col-md-6",
          correlation_plot_scattercorr_ui(
            ns("cor_scatter"),
            title = "Correlation scatter plots",
            info.text = "Pairwise scatter plots for the co-expression of correlated gene pairs across the samples. The straight line correspond to the (linear) regression fit.",
            caption = "Scatter plots of gene expression of top correlated genes.",
            height = c("calc(100vh - 200px)", TABLE_HEIGHT_MODAL),
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
          correlation_plot_cor_graph_ui(
            ns("cor_graph"),
            title = "Partial correlation network",
            info.text = "Red edges correspond to negative correlation, grey edges to positive correlation. Width of the edges is proportional to the absolute partial correlation value of the gene pair.",
            caption = "Partial correlation network around the selected gene.",
            height = c("calc(100vh - 200px)", TABLE_HEIGHT_MODAL),
            width = c(700, "100%"))
        ),
        div(
          class = "col-md-6",
          correlation_plot_correlation_UMAP_ui(
            ns("cor_umap"),
            title = "Correlation UMAP",
            info.text = "Genes that are correlated are generally positioned close to each other. Red corresponds to positive correlation/covariance, blue for negative.",
            caption = "UMAP clustering of genes using covariance as distance metric and colored by correlation (or covariance). ",
            height = c("calc(100vh - 200px)", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    )
  )
  div(
    boardHeader(title = "Correlation analysis", info_link = ns("data_info")),
    tabs
  )
}
