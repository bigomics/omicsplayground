##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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

      ## ---------------------------- panel1 ------------------------------------------
      shiny::tabPanel(
        "FC correlation",
        bslib::layout_column_wrap(
          width = 1/2,
          height = "calc(100vh - 180px)",
          bslib::layout_column_wrap(
            width = 1,
            connectivity_plot_FCFCplots_ui(
              ns("FCFCplots"),
              label = "a",
              title = "FC scatter plots",
              info.text = "The FC-FC scatter plot provides a pairwise scatterplot of logFC fold-change profiles for the selected contrasts. The main purpose of this panel is to identify similarity or dissimilarity between selected contrasts. The scatter plot is interactive and shows information of each gene by hovering over it with the mouse.",
              caption = "Foldchange scatterplot of the selected contrast against a selected pairwise comparison from the public database.",
              height = c("50%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            connectivity_table_similarity_scores_ui(
              ns("connectivityScoreTable"),
              title = "Similarity scores",
              info.text = "Normalized enrichment scores (NES) and Pearson correlation (rho) of reference profiles with respect to the currently selected contrast. The top 100 up/down genes are considered for the calculation of rho or NES. The score is calculated as rho^2*NES. Highlighting a specific dataset will change the FC-FC scatterplot accordingly.",
              caption = "Table displaying the most correlated gene expression profiles from public datasets to the selected contrast. Select a profile to visualise the corresponding scatter plot.",
              height = c("50%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%"),
              label = "b"
            )
          ),
          connectivity_plot_scatterPlot_ui(
            ns("scatterPlot"),
            title = "FC-FC scatterplot",
            info.text = "Scatter plots of gene expression foldchange values between two contrasts. Foldchanges that are similar show high correlation, i.e. are close to the diagonal. You can switch to enrichment type plots in the plot settings.",
            caption = "Scatter plots displaying the public profiles most correlated to the selected contrast by fold-change.",
            label = "c",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      ),

      ## ---------------------------- panel2 ------------------------------------------
      shiny::tabPanel(
        "FC heatmap",
        bslib::layout_column_wrap(
          width = 1,
          height = "calc(100vh - 180px)",
          bslib::layout_column_wrap(
            width = 1,
            height = "35%",
            style = htmltools::css(grid_template_columns = "7fr 5fr"),
            connectivity_plot_cumFCplot_ui(
              id = ns("cumFCplot"),
              title = "Cumulative foldchange",
              info.text = "The barplot visualizes the cumulative foldchange between the top-10 most similar profiles. Genes that are frequently shared with high foldchange will show a higher cumulative score. You can choose between signed or absolute foldchange in the options.",
              caption = "Barplot visualising the most dysregulated genes across the 10 most similar gene expression profiles to the queried contrast.",
              label = "a",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            connectivity_plot_cumEnrichmentPlot_ui(
              id = ns("cumEnrichmentPlot"),
              title = "Cumulative enrichment",
              info.text = "Gene sets that are frequently shared with high enrichment will show a higher cumulative scores. You can choose between signed or absolute enrichment in the options.",
              caption = "The barplot visualizes the cumulative enrichment of the top-10 most similar profiles.",
              label = "b",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          ),
          bslib::layout_column_wrap(
            width = 1,
            height = "65%",
            connectivity_plot_connectivityHeatmap_ui(
              id = ns("connectivityHeatmap"),
              title = "Connectivity Heatmap",
              info.text = "Contrasts that are similar will be clustered close together.",
              caption = "Heatmap displaying the similarity of the selected contrast with gene expression profiles from public datasets",
              label = "c",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),
      ## ----------------------------- panel 3 -------------------------
      shiny::tabPanel(
        "Meta-graph",
        bslib::layout_column_wrap(
          width = 1/2,
          height = "calc(100vh - 180px)",
          connectivity_plot_leadingEdgeGraph_ui(
            id = ns("leadingEdgeGraph"),
            title = "Leading-edge graph",
            info.text = "The edge width corresponds to the number of signatures that share that pair of genes in their top differentially expressed genes. ",
            caption = "Network of shared leading-edge genes between top-N most similar signatures.",
            label = "a",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          ),
          connectivity_plot_enrichmentGraph_ui(
            id = ns("enrichmentGraph"),
            title = "Enrichment graph",
            info.text = "The edge width corresponds to the number of signatures that share that pair of genesets in their top enriched genesets. In the plot options you can set the threshold the edges.",
            caption = "Network of shared enriched genesets between top-N most similar signatures.",
            label = "b",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      ),
      shiny::tabPanel(
        "Experiment clustering",
        bslib::layout_column_wrap(
          width = 1/2,
          height = "calc(100vh - 180px)",
          connectivity_plot_connectivityMap_ui(
            ns("connectivityMap"),
            title = "Connectivity map",
            info.text = "The Connectivity Map shows the similarity of the contrasts profiles as a t-SNE plot. Contrasts that are similar will be clustered close together, contrasts that are different are placed farther away.",
            caption = "Connectivity map placing a selected pairwise comparison expression profile in the context of  collected public datasets.",
            label = "a",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          ),
          connectivity_table_similarity_scores_ui(
            ns("connectivityScoreTable2"),
            title = "Similarity scores",
            info.text = "Normalized enrichment scores (NES) and Pearson correlation (rho) of reference profiles with respect to the currently selected contrast. The top 100 up/down genes are considered for the calculation of rho or NES. The score is calculated as rho^2*NES.",
            caption = "Table showing the similarity scores of the gene expression profile of the selected contrast with those obtained from public reference datasets.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    )
  div(
    boardHeader(title = "Similar experiments", info_link = ns("info")),
    tabs
  )
}
