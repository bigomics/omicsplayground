##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

CompareInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(
      shiny::selectInput(ns("contrast1"), "Dataset1:",
        choices = NULL, multiple = TRUE
      ),
      "Select the contrast that you want to compare.",
      placement = "right", options = list(container = "body")
    ),
    shiny::br(),
    withTooltip(shiny::selectInput(ns("dataset2"), "Dataset2:", choices = NULL),
      "Select second dataset to compare.",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::selectInput(ns("contrast2"), NULL, choices = NULL, multiple = TRUE),
      "Select second contrast to compare.",
      placement = "right", options = list(container = "body")
    ),
    shiny::br(),
    withTooltip(shiny::actionLink(ns("options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "right", options = list(container = "body")
    ),
    shiny::br(),
    shiny::conditionalPanel(
      "input.options % 2 == 1",
      ns = ns,
      shiny::br(),
      withTooltip(
        shiny::radioButtons(ns("plottype"), "Plot type:",
          choices = c("volcano", "MA", "scatter", "UMAP1", "UMAP2", "heatmap"),
          selected = "UMAP1", inline = TRUE
        ),
        "Select plot type.",
        placement = "right", options = list(container = "body")
      ),
      shiny::br(),
      withTooltip(
        shiny::radioButtons(ns("hilighttype"), "Highlight genes:",
          choices = c("top scoring", "custom"),
          inline = TRUE
        ),
        "Select highlight type.",
        placement = "right", options = list(container = "body")
      ),
      shiny::conditionalPanel(
        "input.hilighttype == 'custom'",
        ns = ns,
        withTooltip(
          shiny::textAreaInput(ns("genelist"), NULL,
            value = NULL,
            height = "100px", width = "100%",
            rows = 5, placeholder = "Paste your custom gene list"
          ),
          "Paste a custom list of genes to highlight.",
          placement = "right"
        )
      ),
      shiny::br(),
      withTooltip(
        shiny::radioButtons(ns("ntop"), "ntop",
          choices = c(10, 20, 40, 100),
          selected = 20, inline = TRUE
        ),
        "number of top genes to show",
        placement = "right", options = list(container = "body")
      )
    )
  )
}

CompareUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- "calc(100vh - 180px)"
  tabH  <- "70vh"
  
  tabs <- shiny::tabsetPanel(
    id = ns("tabs1"),
    shiny::tabPanel(
      "Compare expression",
      bslib::layout_column_wrap(
        width = 1/2,
        height = fullH,
        compare_plot_compare1_ui(
          id = ns("dt1"),
          width = c("auto", "100%"),
          height = c("100%", "70vh")
        ),
        compare_plot_compare2_ui(
          id = ns("dt2"),
          width = c("auto", "100%"),
          height = c("100%", "70vh")
        )
      )
    ),
    shiny::tabPanel(
      "Foldchange",
      bslib::layout_column_wrap(
        width = 1/2,
        height = fullH,
        compare_plot_fc_correlation_ui(
          id = ns("fcfcplot"),
          width = c("auto", "100%"),
          height = c("100%", "70vh")
        ),
        bslib::layout_column_wrap(
          width = 1,
          compare_plot_cum_fc1_ui(
            id = ns("cumfcplot1"),
            width = c("auto", "100%"),
            height = c("100%", "70vh"),
            label = "b"
          ),
          compare_plot_cum_fc2_ui(
            id = ns("cumfcplot2"),
            width = c("auto", "100%"),
            height = c("100%", "70vh"),
            label = "c"
          )
        )
      )
    ),
    shiny::tabPanel(
      "Gene Correlation",
      bslib::layout_column_wrap(
        width = 1/2,
        height = fullH,
        bslib::layout_column_wrap(
          width = 1,
          compare_plot_expression_ui(
            id = ns("multibarplot"),
            height = c("70%", TABLE_HEIGHT_MODAL)
          ),
          compare_table_corr_score_ui(
            id = ns("score_table"),
            height = c("30%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")            
          )
        ),
        compare_plot_gene_corr_ui(
          id = ns("genecorr"),
          height = c("100%", TABLE_HEIGHT_MODAL),          
        )
      )
    )
  )

  
  div(
    boardHeader(title = "Compare datasets", info_link = ns("info")),
    tabs
  )
}
