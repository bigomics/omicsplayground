##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ConnectivityInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
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
    shiny::br(),
    bslib::accordion(
      id = ns("compare_accordion"),
      open = FALSE,
      bslib::accordion_panel(
        "Options",
        icon = icon("cog", lib = "glyphicon"),
        withTooltip(
          shiny::checkboxInput(ns("hideclustcontrasts"), "hide cluster contrasts", TRUE),
          "Hide cluster contrasts.",
          placement = "right", options = list(container = "body")
        ),
        withTooltip(
          shiny::checkboxInput(ns("abs_score"), "abs.score", TRUE),
          "Use absolute score value",
          placement = "right", options = list(container = "body")
        ),
        shiny::hr(),
        shiny::tags$head(shiny::tags$style("#cmap-genelist.form-control {font-size:11px !important;padding:3px;height:200px;}")),
        withTooltip(
          shiny::textAreaInput(ns("genelist"), tspan("Select genes:"),
            rows = 15
          ),
          "Paste a list of genes that defines your signature. By default, the top50 most (absolute) differentially expressed genes (by logFC) are chosen for the selected comparison."
        ),
        shiny::radioButtons(
          ns("genelist_ntop"),
          tspan("ngenes:"),
          choices = c(10, 50, 100),
          sel = 50, inline = TRUE
        )
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
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 181px)",
        bs_alert("Compare different experiments by correlating their fold-change signatures. Highly correlated logFC signatures suggest similar experiments."),
        bslib::layout_columns(
          col_widths = c(6, 6),
          height = "calc(100vh - 181px)",
          bslib::layout_columns(
            col_widths = 12,
            connectivity_plot_FCFCplots_ui(
              ns("FCFCplots"),
              label = "a",
              title = "FC scatter plots",
              info.text = "Scatter plots of gene expression foldchange values between two contrasts. Foldchanges that are similar show high correlation, i.e. are close to the diagonal. You can switch to enrichment type plots in the plot settings.",
              caption = "Scatter plots displaying the public profiles most correlated to the selected contrast by fold-change.",
              height = c("50%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            connectivity_table_similarity_scores_ui(
              ns("connectivityScoreTable"),
              title = "Similarity scores",
              info.text = "Connectivity score (score), Pearson correlation (rho), normalized enrichment scores (NES), odd ratio (odd.ratio), and cosine distance (tau) of reference profiles with respect to the currently selected contrast. The top 100 up/down genes are considered for the calculations. The score is calculated as abs(rho)*NES*odd.ratio*abs(tau). Highlighting a specific dataset will change the FC-FC scatterplot accordingly.",
              caption = "Table displaying the most correlated gene expression profiles from public datasets to the selected contrast. Select a profile to visualise the corresponding scatter plot.",
              height = c("50%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%"),
              label = "b"
            )
          ),
          connectivity_plot_scatterPlot_ui(
            ns("scatterPlot"),
            title = "FC-FC scatterplot",
            info.text = "The FC-FC scatter plot provides a pairwise scatterplot of logFC fold-change profiles for the selected contrasts. The main purpose of this panel is to identify similarity or dissimilarity between selected contrasts. The scatter plot is interactive and shows information of each gene by hovering over it with the mouse.",
            caption = "Foldchange scatterplot of the selected contrast against a selected pairwise comparison from the public database.",
            label = "c",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    ),

    ## ---------------------------- panel2 ------------------------------------------

    shiny::tabPanel(
      "FC Heatmap",
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 181px)",
        bs_alert("Compare the fold-change of similar signatures across different experiments."),
        bslib::layout_columns(
          col_widths = c(3, 9),
          height = "calc(100vh - 180px)",
          connectivity_plot_connectivityMap_ui(
            ns("connectivityMap"),
            title = "Connectivity map",
            info.text = "The Connectivity Map shows the similarity of logFC signatures as a t-SNE plot. Signatures that are similar will be clustered close together, signatures that are different are placed farther away.",
            caption = "Connectivity map placing a selected pairwise comparison expression profile in the context of  collected public datasets.",
            label = "a",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          ),
          connectivity_plot_connectivityHeatmap_ui(
            id = ns("connectivityHeatmap"),
            title = "Connectivity Heatmap",
            info.text = "Contrasts that are similar will be clustered close together.",
            caption = "Heatmap displaying the logFC of the selected contrast with most similar gene expression profiles from public datasets",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        ),
        bslib::layout_columns(
          col_widths = 12,
          connectivity_table_foldchange_ui(
            ns("connectivityFoldchangeTable"),
            title = "Fold-change table",
            info.text = "Gene expression fold-changes (log2FC) of similar signatures across different experiments.",
            caption = "Gene expression fold-changes (log2FC) of similar signatures across different experiments.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    ),

    ## ---------------------------- panel3 ------------------------------------------
    shiny::tabPanel(
      "Meta-network",
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 181px)",
        bslib::layout_columns(
          col_widths = c(6, 6),
          height = "35%",
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
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        ),
        bslib::layout_columns(
          col_widths = c(6, 6),
          height = "35%",
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
        )
      )
    ),

    ## ---------------------------- panel4 ------------------------------------------
    shiny::tabPanel(
      "AI Summary",
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 181px)",
        AISummaryCardUI(
          ns("connectivityAISummary"),
          title = "AI Connectivity Summary",
          info.text = "AI-generated summary of the connectivity analysis results for the selected contrast.",
          caption = "AI-generated connectivity summary.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    )
  )

  ## returned UI object
  div(
    boardHeader(title = "Similar experiments", info_link = ns("info")),
    tabs
  )
}
