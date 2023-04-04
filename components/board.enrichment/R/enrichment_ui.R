##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

EnrichmentInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(shiny::selectInput(ns("gs_contrast"), "Contrast:", choices = NULL),
      "Select a contrast of interest for the analysis.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("gs_features"), "Gene set collection:", choices = NULL, multiple = FALSE),
      "Choose a specific gene set collection for the analysis.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("gs_fdr"), "FDR", c(1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1), selected = 0.2),
      "Set the false discovery rate (FDR) threshold.",
      placement = "top"
    ),
    withTooltip(
      shiny::selectInput(ns("gs_lfc"), "logFC threshold",
        choices = c(0, 0.05, 0.1, 0.2, 0.5, 1, 2), selected = 0
      ),
      "Set the logarithmic fold change (logFC) threshold.",
      placement = "top"
    ),
    withTooltip(shiny::actionLink(ns("gs_options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(), br(),
    shiny::conditionalPanel(
      "input.gs_options % 2 == 1",
      ns = ns,
      shiny::tagList(
        withTooltip(shiny::checkboxInput(ns("gs_showall"), "Show all genesets", FALSE),
          "Enbale significant genes filtering. Display only significant genesets in the table.",
          placement = "top", options = list(container = "body")
        ),
        withTooltip(shiny::checkboxGroupInput(ns("gs_statmethod"), "Statistical methods:", choices = NULL),
          "Select a method or multiple methos for the statistical test.",
          placement = "right", options = list(container = "body")
        ),
        withTooltip(shiny::checkboxInput(ns("gs_top10"), "top 10 gene sets", FALSE),
          "Display only top 10 differentially enirhced gene sets (positively and negatively) in the <b>enrihcment analysis</b> table.",
          placement = "top", options = list(container = "body")
        ),
      )
    )
  )
}

EnrichmentUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 800
  rowH <- 420 ## row height of panels
  imgH <- 340 ## height of images
  tabV <- "70vh" ## height of tables
  tabH <- 340 ## row height of panels
  tabH <- "80vh" ## height of tables

  tabs <- tagList(
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "Top enriched",
        div(
          class = "row",
          div(
            class = "col-md-6",
            enrichment_plot_top_enrich_gsets_ui(
              ns("topEnriched"),
              height = c(imgH, TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          ),
          div(
            class = "col-md-6",
            enrichment_plot_freq_top_gsets_ui(
              ns("topEnrichedFreq"),
              height = c(imgH, TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),
      shiny::tabPanel(
        "Plots",
        div(
          class = "row",
          div(
            class = "col-md-3",
            enrichment_plot_volcano_ui(
              ns("subplot_volcano"),
              height = c(imgH, 700),
              width = c("auto", "100%")
            )
          ),
          div(
            class = "col-md-3",
            enrichment_plot_barplot_ui(
              ns("subplot_barplot"),
              height = c(imgH, 700),
              width = c("auto", 900)
            )
          ),
          div(
            class = "col-md-3",
            enrichment_plot_geneplot_ui(
              ns("subplot_geneplot"),
              height = c(imgH, 700),
              width = c("auto", 900)
            )
          ),
          div(
            class = "col-md-3",
            enrichment_plot_scatter_ui(
              ns("subplot_scatter"),
              height = c(imgH, 700),
              width = c("auto", 900)
            )
          )
        )
      ),
      shiny::tabPanel(
        "Compare",
        enrichment_plot_compare_ui(
          ns("compare"),
          height = c(imgH, 450),
          width = c("auto", 1500)
        )
      ),
      shiny::tabPanel(
        "Volcano (all)",
        shiny::fillCol(
          height = 420,
          flex = c(1, NA),
          enrichment_plot_volcanoall_ui(
            ns("volcanoAll"),
            height = c(imgH, 450),
            width = c("auto", "100%")
          )
        )
      ),
      shiny::tabPanel(
        "Volcano (methods)",
        enrichment_plot_volcanomethods_ui(
          ns("volcanoMethods"),
          height = c(imgH, 450),
          width = c("auto", "100%")
        )
      )
    ),
    shiny::tabsetPanel(
      id = ns("tabs2"),
      shiny::tabPanel(
        "Table",
        div(
          class = "row",
          div(
            class = "col-md-7",
            enrichment_table_enrichment_analysis_ui(
              ns("gseatable"),
              width = c("100%", "100%"),
              height = c(285, TABLE_HEIGHT_MODAL)
            )
          ),
          div(
            class = "col-md-5",
            enrichment_table_genes_in_geneset_ui(
              ns("genetable"),
              height = c(285, TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            )
          )
        )
      ),
      shiny::tabPanel(
        "Foldchange (all)",
        enrichment_table_gset_enrich_all_contrasts_ui(
          ns("fctable"),
          height = c(295, TABLE_HEIGHT_MODAL),
          width = c("100%", "100%")
        )
      ),
      shiny::tabPanel(
        "FDR table",
        enrichment_table_n_sig_gsets_ui(
          ns("FDRtable"),
          height = c(295, TABLE_HEIGHT_MODAL),
          width = c("100%", "100%")
        )
      )
    )
  )
  div(
    boardHeader(title = "Geneset enrichment", info_link = ns("gs_info")),
    tabs
  )
}
