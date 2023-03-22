##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

ConnectivityInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(
      shiny::selectInput(ns("contrast"), "Contrast:",
        choices = NULL, multiple = FALSE
      ),
      "Select the contrast that you want to compare.",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::selectInput(ns("sigdb"), "Signature DB:", choices = NULL),
      "Select reference signature database.",
      placement = "right", options = list(container = "body")
    ),
    shiny::selectInput(ns("connectivityScoreTable_qsig"), "threshold (padj)",
      c(0.01, 0.05, 0.2, 1),
      selected = 1
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
      withTooltip(
        shiny::checkboxInput(ns("hideclustcontrasts"), "hide cluster contrasts", TRUE),
        "Hide cluster contrasts.",
        placement = "right", options = list(container = "body")
      ),
      withTooltip(
        shiny::checkboxInput(ns("abs_score"), "abs.score", TRUE),
        "Use absolute score value",
        placement = "right", options = list(container = "body")
      )
    )
  )
}

ConnectivityUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

    tabs <- shiny::tabsetPanel(
      id = ns("tabs1"), 

      ## ---------------------------- panels1 ------------------------------------------
      shiny::tabPanel(
        "FC correlation",
        div(
          class = "row",
          div(
            class = "col-md-6",
            connectivity_plot_FCFCplots_ui(
              ns("FCFCplots"), label = "a",
              height = c(350, 600),
              width = c("auto", 1280)
            ),
            connectivity_table_similarity_scores_ui(
              ns("connectivityScoreTable"),
              height = c("10vh", TABLE_HEIGHT_MODAL),
              width = c("auto", "90%"),
              label = "b"
            )
          ),
          div(
            class = "col-md-6",
            connectivity_plot_scatterPlot_ui(
              ns("scatterPlot"), label = "c")
            )
        ),
        tags$div(
          class = "caption",

          HTML(paste(
            "<b>(a)</b> <b>FC-FC scatter plots.</b> Scatter plots of fold-changes values between two contrasts.
            Fold-changes that are similar show high correlation, i.e. are close to the diagonal.",
            "<b>(b)</b> <b>Similarity scores.</b> Normalized enrichment scores
            (NES) and Pearson correlation (rho) of reference profiles with respect
            to the currently selected contrast. The top 100 up/down genes are
            considered for the calculation of rho or NES. The score is calculated
            as rho^2*NES.",
            "<b>(c)</b> <b>Interactive FC-FC scatter plot.</b> Scatterplots of two
            logFC differential expression profiles for selected contrasts. Similar
            profiles will show high correlation with points close to the diagonal."
          ))

        )
      ),
      shiny::tabPanel(
        "FC heatmap",
        div(
          class = "row",
          div(
            class = "col-md-7",
            connectivity_plot_cumFCplot_ui(ns("cumFCplot"),
              label = "a",
              height = c(300, 600),
              width = c("auto", 1300)
            )
          ),
          div(
            class = "col-md-5",
            connectivity_plot_cumEnrichmentPlot_ui(ns("cumEnrichmentPlot"),
              label = "b",
              height = c(300, 600),
              width = c("auto", 1000)
            )
          )
        ),
        shiny::br(),
        connectivity_plot_connectivityHeatmap_ui(ns("connectivityHeatmap"), label = "c"),
        tags$div(
          class = "caption",
          HTML(paste(
            "<b>(a)</b> <b>Meta-foldchange.</b> The barplot visualizes the
            cumulative foldchange between the top-10 most similar profiles.",
            "<b>(b)</b> <b>Meta-enrichment.</b> The barplot visualizes the
            cumulative enrichment of the top-10 most similar profiles.",
            "<b>(c)</b> <b>Connectivity Heatmap.</b> Similarity of the contrasts
            profiles as a heatmap. Contrasts that are similar will be clustered
            close together."
          ))
        )
      ),
      shiny::tabPanel(
        "Meta-graph",
        div(
          class = "row",
          div(
            class = "col-md-6",
            connectivity_plot_leadingEdgeGraph_ui(ns("leadingEdgeGraph"), label = "a")
          ),
          div(
            class = "col-md-6",
            connectivity_plot_enrichmentGraph_ui(ns("enrichmentGraph"), label = "b")
          )
        ),
        tags$div(
          class = "caption",
          HTML(paste(
            "<b>(a)</b> <b>Leading-edge graph.</b> Network of shared leading-edge genes between top-N most similar signatures. The edge width corresponds to the number of signatures that share that pair of genes in their top differentially expressed genes.",
            "<b>(b)</b> <b>Enrichment graph.</b> Network of shared enrichmed genesets between top-N most similar signatures. The edge width corresponds to the number of signatures that share that pair of genesets in their top enriched genesets."
          ))
        )
      ),
      shiny::tabPanel(
        "Experiment clustering",
        div(
          class = "row",
          div(
            class = "col-md-6",
            connectivity_plot_connectivityMap_ui(ns("connectivityMap"), label = "a")
          ),
          div(
            class = "col-md-6",
            connectivity_table_similarity_scores_ui(
              ns("connectivityScoreTable2"),
              height = c(660, TABLE_HEIGHT_MODAL),
              width = c("auto", "90%")
            )
          )
        ),
        tags$div(
          class = "caption",
          HTML(paste(
            "<b>(a)</b> <b>Connectivity Map.</b> The CMap shows the similarity of the contrasts as a t-SNE plot. Contrasts that are similar will be clustered close together, contrasts that are different are placed farther away.",
            "<b>(b)</b> <b>Similarity scores.</b> Normalized enrichment scores (NES) and Pearson correlation (rho) of reference profiles with respect to the currently selected contrast. The top 100 up/down genes are considered for the calculation of rho or NES. The score is calculated as rho^2*NES. "
          ))
        )
      )
    )
  div(
    boardHeader(title = "Similar experiments", info_link = ns("info")),
    div(selector_default(ns("hide_caption"), label = "Show captions")),
    tabs
  )
}
