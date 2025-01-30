##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

CorrelationInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),

    ## data set parameters
    withTooltip(shiny::selectInput(ns("gene"), tspan("Gene:"), choices = NULL),
      "Choose a gene for the correlation analysis.",
      placement = "top"
    ),
    shiny::br(),
    withTooltip(shiny::selectInput(ns("cor_filter"), tspan("Filter genes:"), choices = NULL, multiple = FALSE),
      "Filter gene features.",
      placement = "top"
    ),
    # shiny::conditionalPanel(
    #   "input.cor_filter == '<custom>'",
    #   ns = ns,
      withTooltip(
        shiny::textAreaInput(ns("cor_customfeatures"),
          NULL,
          value = NULL,
          height = "100px", width = "100%",
          rows = 5
        ),
        "Paste a custom list of genes to be used as features.",
        placement = "top"
      ),
    #),
    shiny::br(), shiny::br(), shiny::br(), shiny::br(),
    withTooltip(shiny::actionLink(ns("adv_options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(), shiny::br(),
    # shiny::conditionalPanel(
    #   "input.adv_options % 2 == 1",
    #   ns = ns,
      withTooltip(
        shiny::radioButtons(ns("pcor_ntop"), tspan("Nr. of genes to compute partial correlation."),
          c(50, 100, 250),
          selected = 100, inline = TRUE
      ),
      "Number of top genes to compute partial correlation",
      placement = "top"
    )
    # ),
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
      bslib::layout_columns(
        col_widths = c(6, 6),
        height = "calc(100vh - 180px)",
        bslib::layout_columns(
          col_widths = 12,
          correlation_plot_barplot_ui(
            id = ns("cor_barplot"),
            title = "Top correlated features",
            info.text = "Barplot of the highest correlated features in respect to the selected {Gene}. The dark grey bars correspond to the 'partial correlation', which relation is computed using the amount of genes specified under Settings > Options {Nr. of genes to compute partial correlation}.",
            info.methods = "The barplot displays the Pearson correlation value (computed using the core R stats package). The 'partial correlation' corrects the correlation value for indirect effects and tries to estimate the amount of direct interaction, it does so using the glasso algorithm (from the glasso R package [1]) to estimate a sparse inverse covariance matrix and derive partial correlations.",
            info.references = list(
              list(
                "Friedman J (2019) glasso: Graphical Lasso: Estimation of Gaussian Graphical Models",
                "https://doi.org/10.32614/CRAN.package.glasso"
              )
            ),
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#correlation-analyses",
            caption = "Barplot showing the highest correlated feature with respect to the selected feature.",
            label = "",
            height = c("50%", "70vh"),
            width = c("auto", "100%")
          ),
          correlation_table_corr_ui(
            id = ns("cor_table"),
            title = "Correlation table",
            info.text = "Statistical results from correlated gene pairs.",
            caption = "Correlation table of correlation and partial correlation with respect to the selected gene. ",
            label = "",
            height = c("50%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        ),
        correlation_plot_scattercorr_ui(
          ns("cor_scatter"),
          title = "Correlation scatter plots",
          info.text = "Scatter plots of the co-expression of correlated gene pairs between the selected {Gene} and the top genes correlated to it (as seen on Top correlated genes plot) across the samples. The straight line correspond to the (linear) regression fit. The samples can be colored using the {Color by} plot setting and the layout of the scatter plots can be configured by using the {Layout} and {Swap XY-axes} plot settings.",
          info.methods = "See Top correlated features",
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#correlation-analyses",
          caption = "Scatter plots of gene expression of top correlated genes.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ),
    shiny::tabPanel(
      "Graph",
      bslib::layout_columns(
        col_widths = c(6, 6),
        height = "calc(100vh - 180px)",
        correlation_plot_cor_graph_ui(
          ns("cor_graph"),
          title = "Partial correlation network",
          info.text = "Partial correlation network graph of the selected {Gene}. The plot has different settings to control the {radius} (each of each node in the network), {pcor treshold} (minimum value of the partial correlation coefficient that must be met for an edge to be displayed between two nodes) and {layout} (algorithm to position the nodes). The 'partial correlation' is computed using the amount of genes specified under Settings > Options {Nr. of genes to compute partial correlation}.",
          info.methods = "Network graph is computed from the partial correlation, which corrects the correlation value (computed using the core R stats package using Pearson correlation) for indirect effects and tries to estimate the amount of direct interaction, it does so using the glasso algorithm (from the glasso R package [1]) to estimate a sparse inverse covariance matrix and derive partial correlations. Grey edges correspond to positive correlation, red edges correspond to negative correlation. The width of the edge is proportional to the absolute partial correlation value of the gene pair.",
          info.references = list(
            list(
              "Friedman J (2019) glasso: Graphical Lasso: Estimation of Gaussian Graphical Models",
              "https://doi.org/10.32614/CRAN.package.glasso"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#correlation-analyses",
          caption = "Partial correlation network around the selected gene.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),
        correlation_plot_correlation_UMAP_ui(
          ns("cor_umap"),
          title = "Correlation UMAP",
          info.text = "UMAP clustering of genes using covariance as distance metric and colored by correlation or covariance, depending on the plot setting {color by}.",
          info.methods = "Uniform Manifold Approximation and Projection (UMAP) is a non-linear dimensionality reduction method that enables visualization of high-dimensional data in a low-dimensional space (it has been computed using the umap R package [1]). Genes that are correlated are generally positioned close to each other. Red corresponds to positive correlation/covariance, blue for negative.",
          info.references = list(
            list(
              "Konopka T (2023) umap: Uniform Manifold Approximation and Projection",
              "https://doi.org/10.32614/CRAN.package.umap"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
          caption = "UMAP clustering of genes using covariance as distance metric and colored by correlation (or covariance). ",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    )
  )

  ## full page
  div(
    boardHeader(title = "Correlation analysis", info_link = ns("data_info")),
    tabs
  )
}
