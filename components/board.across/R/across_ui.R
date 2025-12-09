##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

AcrossInputs <- function(id) {

  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    withTooltip(
      shiny::selectizeInput(
        ns("selected_datasets"),
        "Select datasets:",
        choices = NULL,
        multiple = TRUE,
        options = list(
          placeholder = "All datasets (default)"
        )
      ),
      "Select which datasets to query. Leave empty to query all datasets.",
      placement = "right", options = list(container = "body")
    ),
    shiny::br(),
    withTooltip(
      shiny::selectizeInput(
        ns("selected_genes"),
        "Select genes:",
        choices = NULL,
        multiple = TRUE,
        options = list(
          maxItems = 10,
          placeholder = "Type gene names..."
        )
      ),
      "Select up to 10 genes to query across selected datasets.",
      placement = "right", options = list(container = "body")
    ),
    shiny::br(),
    withTooltip(
      shiny::actionButton(
        ns("query_button"),
        label = "Query",
        class = "btn-outline-primary",
        icon = icon("search")
      ),
      "Click to query selected genes across selected datasets.",
      placement = "right"
    ),
    shiny::br(),
    shiny::br(),
    bslib::accordion(
      id = ns("across_accordion"),
      open = FALSE,
      bslib::accordion_panel(
        "Phenotype Filter",
        icon = icon("filter", lib = "glyphicon"),
        shiny::tags$small(
          shiny::tags$em("Filter samples by phenotype values shared across selected datasets.")
        ),
        shiny::br(),
        shiny::br(),
        withTooltip(
          shiny::selectizeInput(
            ns("filter_phenotype"),
            "Filter by phenotype:",
            choices = NULL,
            multiple = FALSE,
            options = list(
              placeholder = "Select phenotype..."
            )
          ),
          "Select a phenotype column to filter samples. Only phenotypes common to all selected datasets are shown.",
          placement = "right", options = list(container = "body")
        ),
        shiny::conditionalPanel(
          condition = paste0("input['", ns("filter_phenotype"), "'] != ''"),
          withTooltip(
            shiny::selectizeInput(
              ns("filter_values"),
              "Include values:",
              choices = NULL,
              multiple = TRUE,
              options = list(
                placeholder = "All values (default)"
              )
            ),
            "Select which phenotype values to include. Leave empty to include all values.",
            placement = "right", options = list(container = "body")
          )
        ),
        shiny::hr(),
        withTooltip(
          shiny::selectizeInput(
            ns("color_by_phenotype"),
            "Color plots by:",
            choices = c("dataset" = "dataset"),
            selected = "dataset",
            multiple = FALSE
          ),
          "Select a phenotype to color plots by. Default is coloring by dataset.",
          placement = "right", options = list(container = "body")
        )
      ),
      bslib::accordion_panel(
        "Options",
        icon = icon("cog", lib = "glyphicon"),
        withTooltip(
          shiny::radioButtons(
            ns("plot_scale"),
            "Scale:",
            choices = c("linear", "log2"),
            selected = "linear",
            inline = TRUE
          ),
          "Select scale for expression values.",
          placement = "right", options = list(container = "body")
        )
      )
    ),
    shiny::br(),
    shiny::hr(),
    shiny::tags$small(
      shiny::tags$strong("Database info:"),
      shiny::br(),
      shiny::textOutput(ns("db_info"), inline = TRUE)
    )
  )
}

AcrossUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- "calc(100vh - 181px)"
  tabH <- "70vh"

  tabs <- shiny::tabsetPanel(
    id = ns("tabs1"),
    shiny::tabPanel(
      "Expression by Sample",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        bs_alert("View gene expression values across all samples from multiple datasets stored in the TileDB database."),
        bslib::layout_columns(
          col_widths = 12,
          height = fullH,
          across_plot_barplot_ui(
            id = ns("barplot"),
            title = "Expression by Sample",
            info.text = "Bar plot showing expression values for selected genes across all samples. Samples are colored by their source dataset.",
            info.methods = "Expression values are retrieved from a TileDB database containing counts from multiple PGX files. Values are displayed as raw counts or log2-transformed based on the scale setting.",
            info.extra_link = NULL,
            width = c("auto", "100%"),
            height = c("100%", "75vh")
          )
        )
      )
    ),
    shiny::tabPanel(
      "Expression by Dataset",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        bs_alert("Compare gene expression distributions across different datasets using boxplots."),
        bslib::layout_columns(
          col_widths = 12,
          height = fullH,
          across_plot_boxplot_ui(
            id = ns("boxplot"),
            title = "Distribution by Dataset",
            info.text = "Boxplot showing the distribution of expression values for selected genes across different datasets.",
            info.methods = "Expression values are grouped by source dataset and displayed as boxplots. Jittered points can be overlaid to show individual sample values.",
            info.extra_link = NULL,
            width = c("auto", "100%"),
            height = c("100%", "75vh")
          )
        )
      )
    ),
    shiny::tabPanel(
      "Data Table",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        bs_alert("View and download the raw expression data for selected genes."),
        bslib::layout_columns(
          col_widths = 12,
          height = fullH,
        across_table_data_ui(
          id = ns("datatable"),
          title = "Expression Data",
            height = c("100%", "75vh"),
            width = c("auto", "100%")
          )
        )
      )
    )
  )

  div(
    boardHeader(title = "Across Datasets", info_link = ns("info")),
    tabs
  )
}

