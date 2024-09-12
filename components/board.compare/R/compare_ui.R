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
    withTooltip(
      shiny::actionButton(
        ns("compare_button"),
        label = "Update",
        class = "btn-outline-primary",
        icon = icon("refresh")
      ),
      "Click to update the comparison plots.",
      placement = "right"
    ),
    shiny::br(),
    shiny::br(),
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
          choices = c("volcano", "MA", "scatter", "UMAP1", "UMAP2"),
          selected = "UMAP1", inline = TRUE
        ),
        "Select plot type.",
        placement = "right", options = list(container = "body")
      ),
      shiny::br(),
      withTooltip(
        shiny::radioButtons(ns("hilighttype"), tspan("Highlight genes:"),
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
            rows = 5
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
  tabH <- "70vh"

  tabs <- shiny::tabsetPanel(
    id = ns("tabs1"),
    shiny::tabPanel(
      "Compare expression",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        bs_alert("Compare different experiments by correlating their fold-change signatures. Highly correlated logFC signatures suggest similar experiments."),
        bslib::layout_columns(
          col_widths = c(6, 6),
          height = fullH,
          compare_plot_compare1_ui(
            id = ns("dataset1"),
            title = "Dataset 1",
            info.text = "Signature plot for the loaded dataset selected comparison (under {Dataset1}). This plot displays a UMAP clustering by default, however under option settings {Plot type} different methods can be used. Additionally a custom set of genes can be highlighted using {Highlight genes} or the {ntop} scoring. Multiple contrasts can be selected.",
            info.methods = "Uniform Manifold Approximation and Projection (UMAP) is a non-linear dimensionality reduction method that enables visualization of high-dimensional data in a low-dimensional space (it has been computed using the umap R package [1]). Genes that are correlated are generally positioned close to each other. Red corresponds to positive correlation/covariance, blue for negative. Bland-Altman (MA) plot displays mean intensity versus fold-change. Volcano plot displays fold-change versus significance. Scatter plot displays the expression of the genes for the comparison. The heatmap is generated using the ComplexHeatmap R/Bioconductor package [2] on scaled log-expression values (z-score) using euclidean distance and Ward linkage using the fastcluster R package [3]. ",
            info.references = list(
              list(
                "Konopka T (2023) umap: Uniform Manifold Approximation and Projection",
                "https://doi.org/10.32614/CRAN.package.umap"
              ),
              list(
                "Gu Z (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics.",
                "https://doi.org/10.1093/bioinformatics/btw313"
              ),
              list(
                "Müllner D (2013). “fastcluster: Fast Hierarchical, Agglomerative Clustering Routines for R and Python.” Journal of Statistical Software, 53(9), 1–18.",
                "https://doi.org/10.18637/jss.v053.i09"
              )
            ),
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
            width = c("auto", "100%"),
            height = c("100%", "70vh")
          ),
          compare_plot_compare2_ui(
            id = ns("dataset2"),
            title = "Dataset 2",
            info.text = "Signature plot for the selected dataset and comparison (under {Dataset2}). This plot displays a UMAP clustering by default, however under option settings {Plot type} different methods can be used. Additionally a custom set of genes can be highlighted using {Highlight genes} or the {ntop} scoring. Multiple contrasts can be selected.",
            info.methods = "See Dataset 1.",
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
            width = c("auto", "100%"),
            height = c("100%", "70vh")
          )
        )
      )
    ),
    shiny::tabPanel(
      "Foldchange",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        bs_alert("Compare signatures by plotting their fold-changes as pairwise scatterplots. Highly correlated logFC signatures suggest similar experiments. Data comparison between datasets of the same species is done via gene symbol and datasets from diff species is done using the human ortholog"),
        bslib::layout_columns(
          col_widths = c(6, 6),
          height = fullH,
          compare_plot_fcfc_ui(
            id = ns("fcfcplot"),
            title = "FC Correlation",
            info.text = "Scatter plot of the fold-change correlation for the loaded dataset contrast (selected under {Dataset1}) and the second selected dataset and contrast (under {Dataset2}). A custom set of genes can be highlighted using {Highlight genes} or the {ntop} scoring option settings. Multiple contrasts can be selected.",
            info.methods = "Pearson Pairwise correlation is used to compute the correlation between the datasets (using the base stats R package). Highly correlated logFC signatures suggest similar experiments. Data comparison between datasets of the same species is done via gene symbol and datasets from different species is done using the human ortholog. Scatters that are similar show high correlation, i.e. are close to the diagonal.",
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#correlation-analyses",
            width = c("auto", "100%"),
            height = c("100%", "70vh")
          ),
          bslib::layout_columns(
            col_widths = 12,
            compare_plot_cum_fc1_ui(
              id = ns("cumfcplot1"),
              title = "Foldchange (Dataset 1)",
              info.text = "Barplot showing the cumulative fold changes for the loaded dataset contrast (selected under {Dataset1}). Multiple contrasts can be selected.",
              width = c("auto", "100%"),
              height = c("100%", "70vh"),
              label = "b"
            ),
            compare_plot_cum_fc2_ui(
              id = ns("cumfcplot2"),
              title = "Foldchange (Dataset 2)",
              info.text = "Barplot showing the cumulative fold changes for the second selected dataset and contrast (under {Dataset2}). Multiple contrasts can be selected.",
              width = c("auto", "100%"),
              height = c("100%", "70vh"),
              label = "c"
            )
          )
        )
      )
    ),
    shiny::tabPanel(
      "Gene Correlation",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        bs_alert("If the samples are exactly the same in your two datasets, you can plot their gene expression and find highly correlated features, e.g. for genes and proteins."),
        bslib::layout_columns(
          col_widths = c(6, 6),
          height = fullH,
          bslib::layout_columns(
            col_widths = 12,
            compare_plot_expression_ui(
              id = ns("multibarplot"),
              title = "Expression",
              info.text = "Barplots of expression values for the loaded dataset contrast (selected under {Dataset1}) and the second selected dataset and contrast (under {Dataset2}). Multiple contrasts can be selected.",
              height = c("70%", TABLE_HEIGHT_MODAL)
            ),
            compare_table_corr_score_ui(
              id = ns("score_table"),
              height = c("30%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          ),
          compare_plot_genecorr_ui(
            id = ns("genecorr"),
            title = "Gene correlation",
            info.text = "Scatter plots of gene expression correlation between the loaded dataset contrast (selected under {Dataset1}) and the second selected dataset and contrast (under {Dataset2}). Scatters that are similar show high correlation, i.e. are close to the diagonal.",
            height = c("100%", TABLE_HEIGHT_MODAL)
          )
        )
      )
    )
  )


  div(
    boardHeader(title = "Compare datasets", info_link = ns("info")),
    tabs
  )
}
