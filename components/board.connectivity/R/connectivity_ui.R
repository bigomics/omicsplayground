##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

ConnectivityInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    withTooltip(
      shiny::actionLink(ns("cmap_info"), "Info", icon = shiny::icon("info-circle")),
      "Show more information about this module"
    ),
    shiny::hr(), shiny::br(),
    withTooltip(
      shiny::selectInput(ns("cmap_contrast"), "Contrast:",
        choices = NULL, multiple = FALSE
      ),
      "Select the contrast that you want to compare.",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::selectInput(ns("cmap_sigdb"), "Signature DB:", choices = NULL),
      "Select reference signature database.",
      placement = "right", options = list(container = "body")
    ),
    shiny::br(),
    withTooltip(shiny::actionLink(ns("cmap_options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "right", options = list(container = "body")
    ),
    shiny::br(),
    shiny::conditionalPanel(
      "input.cmap_options % 2 == 1",
      ns = ns,
      withTooltip(
        shiny::checkboxInput(ns("cmap_hideclustcontrasts"), "hide cluster contrasts", TRUE),
        "Hide cluster contrasts.",
        placement = "right", options = list(container = "body")
      ),
      withTooltip(
        shiny::checkboxInput(ns("cmap_abs_score"), "abs.score", TRUE),
        "Use absolute score value",
        placement = "right", options = list(container = "body")
      )
    )
  )
}

ConnectivityUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  tagList(
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "FC correlation",
        tags$div(
          HTML(paste(
            "<b>(a)</b> <b>FC scatter plots.</b> Scatter plots of gene expression foldchange values between two contrasts. Foldchanges that are similar show high correlation, i.e. are close to the diagonal.",
            "<b>(b)</b> <b>Similarity scores.</b> Normalized enrichment scores (NES) and Pearson correlation (rho) of reference profiles with respect to the currently selected contrast. The top 100 up/down genes are considered for the calculation of rho or NES. The score is calculated as rho^2*NES.",
            "<b>(c)</b> <b>Pairs plot.</b> Pairwise scatterplots of two differential expression profiles for selected contrasts. Similar profiles will show high correlation with points close to the diagonal."
          ))
        ),
        div(
          class = "row",
          div(
            class = "col-md-6",
            connectivity_plot_cmap_FCFCplots_ui(ns("cmap_FCFCplots"))
          ),
          div(
            class = "col-md-6",
            tableWidget(ns("connectivityScoreTable"))
          )
        ),
        plotWidget(ns("cmapPairsPlot"))
      ),
      shiny::tabPanel(
        "FC heatmap",
        tags$div(
          HTML(paste(
            "<b>(a)</b> <b>Meta-foldchange.</b> The barplot visualizes the cumulative foldchange between the top-10 most similar profiles.",
            "<b>(b)</b> <b>Meta-enrichment.</b> The barplot visualizes the cumulative enrichment of the top-10 most similar profiles.",
            "<b>(c)</b> <b>Connectivity Heatmap.</b> Similarity of the contrasts profiles as a heatmap. Contrasts that are similar will be clustered close together."
          ))
        ),
        div(
          class = "row",
          div(
            class = "col-md-6",
            plotWidget(ns("cumFCplot"))
          ),
          div(
            class = "col-md-6",
            plotWidget(ns("cumEnrichmentPlot"))
          )
        ),
        shiny::br(),
        plotWidget(ns("connectivityHeatmap"))
      ),
      shiny::tabPanel(
        "Meta-graph",
        tags$div(
          HTML(paste(
            "<b>(a)</b> <b>Leading-edge graph.</b> Network of shared leading-edge genes between top-N most similar signatures. The edge width corresponds to the number of signatures that share that pair of genes in their top differentially expressed genes.",
            "<b>(b)</b> <b>Enrichment graph.</b> Network of shared enrichmed genesets between top-N most similar signatures. The edge width corresponds to the number of signatures that share that pair of genesets in their top enriched genesets."
          ))
        ),
        div(
          class = "row",
          div(
            class = "col-md-6",
            plotWidget(ns("leadingEdgeGraph"))
          ),
          div(
            class = "col-md-6",
            plotWidget(ns("enrichmentGraph"))
          )
        )
      ),
      shiny::tabPanel(
        "Experiment clustering",
        tags$div(
          HTML(paste(
            "<b>(a)</b> <b>Connectivity Map.</b> The CMap shows the similarity of the contrasts as a t-SNE plot. Contrasts that are similar will be clustered close together, contrasts that are different are placed farther away.",
            "<b>(b)</b> <b>Similarity scores.</b> Normalized enrichment scores (NES) and Pearson correlation (rho) of reference profiles with respect to the currently selected contrast. The top 100 up/down genes are considered for the calculation of rho or NES. The score is calculated as rho^2*NES. "
          ))
        ),
        div(
          class = "row",
          div(
            class = "col-md-6",
            plotWidget(ns("connectivityMap"))
          ),
          div(
            class = "col-md-6",
            tableWidget(ns("connectivityScoreTable2"))
          )
        )
      )
    )
  )
}
