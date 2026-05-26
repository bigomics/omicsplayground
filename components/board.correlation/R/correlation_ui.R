##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

CorrelationInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    withTooltip(shiny::selectInput(ns("gene"), tspan("Gene:"), choices = NULL),
      "Choose a gene for the correlation analysis.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("cor_filter"), tspan("Filter genes:"), choices = NULL, multiple = FALSE),
      "Filter gene features.",
      placement = "top"
    ),
    shiny::conditionalPanel(
      "input.cor_filter == '<custom>'",
      ns = ns,
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
    ),
    shiny::br(),
    bslib::accordion(
      id = ns("cor_accordion"),
      open = FALSE,
      bslib::accordion_panel(
        "Options",
        icon = icon("cog", lib = "glyphicon"),
        withTooltip(
          shiny::radioButtons(ns("pcor_ntop"), tspan("Nr. of genes to compute partial correlation."),
            c(50, 100, 250),
            selected = 100, inline = TRUE
          ),
          "Number of top genes to compute partial correlation",
          placement = "top"
        )
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
      bslib::layout_columns(
        col_widths = c(6, 6),
        height = "calc(100vh - 181px)",
        bslib::layout_columns(
          col_widths = 12,
          correlation_plot_barplot_ui(
            id = ns("cor_barplot"),
            title = "Top correlated features",
            info.text = "Barplot of features exhibiting highest correlation with the selected {Gene}. The dark blue bars show the partial correlation coefficient. The green bars show the canonical Pearson's correlation coefficient. The number of features (N) used to compute partial correlation can be set in Settings > Options {Nr. of genes to compute partial correlation}.",
            info.methods = "Partial correlation coefficient between any two features across samples is computed using the Glasso algorithm from the glasso R package. Specifically, canonical Pearson's correlation between the selected {Gene} and all other available features is first computed on the log2-transformed (and normalized, by default) data using the core R stats package. The top N features are then identified as those exhibiting the highest average squared correlation value. A covariance matrix is calculated on the log2-transformed and normalized transposed data matrix (samples x top N features). Glasso is then applied on the covariance matrix to estimate a sparse inverse covariance matrix. A lasso (L1) regularization of 0.01 is applied. Using the function cov2cor from the Matrix R package, the partial correlation coefficient matrix is derived from the glasso-estimated inverse covariance matrix. Compared to canonical Pearson's correlation, partial correlation coefficients measure effect size and direction of the linear relationship between two variables while adjusting for potential confounding effects (indirect effects) of other variables, therefore providing an estimate of potential direct interactions.",
            info.references = list(
              list(
                "Friedman J (2019) glasso: Graphical Lasso: Estimation of Gaussian Graphical Models",
                "https://doi.org/10.32614/CRAN.package.glasso"
              )
            ),
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#correlation-analyses",
            caption = "Barplot showing the Pearson's and partial correlation coefficients for the features exhibiting highest correlation with selected feature.",
            label = "",
            height = c("50%", "70vh"),
            width = c("auto", "100%")
          ),
          correlation_table_corr_ui(
            id = ns("cor_table"),
            title = "Correlation table",
            info.text = "Statistical results from canonical Pearson's correlation and partial correlation analysis between features. The column 'cor' reports the canonical Pearson's correlation. The column 'pcor' reports the partial correlation coefficients derived from the glasso-estimated inverse covariance matrix. For details on how top correlated features are inferred, please refer to the information provided for the 'Top correlated features' barplot.",
            caption = "",
            label = "",
            height = c("50%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        ),
        correlation_plot_scattercorr_ui(
          ns("cor_scatter"),
          title = "Correlation scatter plots",
          info.text = "Scatter plots of log2-transformed (and normalized, by default) expression values of the selected {Gene} and the top correlated features across samples. Top correlated features are also displayed in the 'Top correlated features' barplot. In each scatter plot, each dot corresponds to a sample. The straight line corresponds to the linear regression fit. The samples can be colored using the {Color by} plot setting. The layout of the scatter plots can be configured by using the {Layout} and {Swap XY-axes} plot settings. Optionally, a correlation line (r=1) can also be displayed in the plot.",
          info.methods = "For details on how top correlated features are inferred, please refer to the information provided for the 'Top correlated features' barplot.",
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#correlation-analyses",
          caption = "Scatter plots of log2-transformed (and normalized, by default) expression values of the selected {Gene} and the top correlated features across samples",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ),
    shiny::tabPanel(
      "Graph",
      bslib::layout_columns(
        col_widths = c(6, 6),
        height = "calc(100vh - 181px)",
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
