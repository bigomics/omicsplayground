##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

LasagnaInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::selectInput(ns("contrast"), "Select comparison", choices = NULL),
    shiny::selectInput(ns("layers"), "Show layers:", choices = NULL, multiple = TRUE),
    shiny::conditionalPanel(
      "input.layers && input.layers.indexOf('gset') > -1",
      ns = ns,
      shiny::selectInput(ns("gsfilter"), "Geneset filter:", choices = NULL)
    ),
    shiny::br(),
    shiny::actionButton(ns("updateplots"), "Update plots", size = "xs", icon = icon("refresh")),
    shiny::br(),
    shiny::br(),
    bslib::accordion(
      id = ns("clust_options"),
      open = FALSE,
      bslib::accordion_panel(
        "Cluster options",
        icon = icon("cog", lib = "glyphicon"),
        shiny::tagList(
          shiny::radioButtons(ns("clustmethod"), "Cluster method:",
            choices = c("pca", "tsne", "umap"), selected = "pca", inline = TRUE
          )
        )
      )
    ),
    bslib::accordion(
      id = ns("mpartite_options"),
      open = FALSE,
      bslib::accordion_panel(
        "Network options",
        icon = icon("cog", lib = "glyphicon"),
        shiny::tagList(
          ##shiny::checkboxInput(ns("top50"),"top 50",TRUE),
          shiny::checkboxInput(ns("consensus"),"consensus",FALSE),          
          shiny::checkboxInput(ns("sp_weight"),"SP weighting",FALSE),
          shiny::sliderInput(ns("minrho"),"Edge threshold:",0,0.95,0.5,0.05),
          shiny::hr(),
          shiny::radioButtons(ns("node_value"), "Node value:",
            choices = c("logFC", "rho"), selected = "logFC", inline = TRUE
          ),
          shiny::hr(),
          shiny::radioButtons(ns("ntop"), "Number of nodes:", c(50, 1000), inline = TRUE)
        )
      )
    )
  )
}

my_navset_card_tab <- function(...) {
  htmltools::tagAppendAttributes(
    bslib::navset_card_tab(
      ...,
      tags$style(HTML("@media (min-width: 1200px) {.root_navset { height: calc(100vh - 36px); }}"))
    ),
    class = "root_navset border-0"
  )
}

MPARTITE_INFO <- "The <b>Multi-partite graph</b> shows the correlation structure between multiple sets of features. The color of the edges correspond to positive (purple) and negative (yellow) correlation. Thicker edges mean higher correlation. The sizes of the circles represent the page-rank centrality of the feature. The log2FC is indicated for the chosen comparison. The node color corresponds to up (red) and down (blue) regulation."

NETWORK_INFO <- "Multi-type network"

LasagnaUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    boardHeader(title = "LASAGNA", info_link = ns("info")),
    shiny::tabsetPanel(
      # my_navset_card_tab(
      id = ns("tabs"),
      # title = "LASAGNA",

      ## ----------------------------------------------------------------
      shiny::tabPanel(
        ## bslib::nav_panel(
        "Multi-layer model",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto", 1),
          bs_alert(HTML("<b>LASAGNA</b> is a stacked layer model for multi-omics where each layer corresponds to a data type. The acronym stands for a <u>L</u>ayered <u>A</u>pproach to <u>S</u>imultaneous <u>A</u>nalysis of <u>G</u>enomic and <u>N</u>etwork <u>A</u>ssociations'.")),
          bslib::layout_columns(
            col_widths = c(6, 6),
            height = "calc(100vh - 180px)",
            mofa_plot_lasagna3D_ui(
              ns("lasagna"),
              title = "Multi-layer model",
              info.text = "Layered Approach to Simultaneous Analysis of Genomic and Network Associations ('LASAGNA'). The LASAGNA plot is a stacked layer plot to visualize multi-omics data. Specifically, each layer shows a data type-specific UMAP. LASAGNA just shows the datatype-specific UMAPs overlayed.",
              info.references = list(list("Melville J (2024). “uwot: The Uniform Manifold Approximation and Projection (UMAP) Method for Dimensionality Reduction.”, CRAN.", "https://doi.org/10.32614/CRAN.package.uwot")),
              caption = "Layered Approach to Simultaneous Analysis of Genomic and Network Associations ('LASAGNA'). The LASAGNA plot is a stacked layer plot to visualize data type-specific UMAPs overlayed.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_plot_lasagna_clustering_ui(
              ns("clusters"),
              title = "Feature UMAP per datatype",
              info.text = "Feature-level UMAP clustering per data type. Visually explore signatures of distinct datatypes across the same set of samples. Feature-level clustering enables assessment of how the distinct data types/modalities define distinct (functional) groups. This analysis may reveal that distinct data types capture different heterogeneities in the data, potentially associated with unique biological functions. On the other end, similar clustering patterns between distinct data types may indicate shared regulation.",
              info.methods = "Feature-level clustering is determined by Uniford Manifold Approximation and Projection (UMAP) applied to each data type separately. UMAP is computed using the uwot R package. The colors in the UMAP reflect the low-to-high correlation with the selected comparison to explore the impact of different conditions.",
              info.references = list(list("Melville J (2024). “uwot: The Uniform Manifold Approximation and Projection (UMAP) Method for Dimensionality Reduction.”, CRAN.", "https://doi.org/10.32614/CRAN.package.uwot")),
              caption = "Feature-level UMAP clustering per data type. Visually explore signatures of distinct datatypes across the same set of samples. Feature-level clustering is determined by UMAP applied to each data type separately. The colors in the UMAP reflect the low-to-high correlation with the selected comparison ",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),

      ## ##----------------------------------------------------------------
      shiny::tabPanel(
        "Multi-partite graph",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          row_heights = c("auto", 1),
          bs_alert(HTML(MPARTITE_INFO)),
          bslib::layout_columns(
            # col_widths = c(8,4),
            col_widths = c(12),
            height = "calc(100vh - 180px)",
            mofa_plot_lasagna_partite_ui(
              ns("lasagnaPartite"),
              title = "Multi-partite graph",
              ## caption = NULL,
              info.text = MPARTITE_INFO,
              info.references = NULL,
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),
      ## end tabPanel
      ## ##----------------------------------------------------------------
      shiny::tabPanel(
        "Multi-type network",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          row_heights = c("auto", 1),
          bs_alert(HTML(NETWORK_INFO)),
          bslib::layout_columns(
            col_widths = c(12),
            height = "calc(100vh - 180px)",
            mofa_plot_lasagna_network_ui(
              ns("lasagnaNetwork"),
              title = "Multi-type network",
              info.text = NETWORK_INFO,
              info.references = NULL,
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ) ## end tabPanel
    ) ## end tabsetPanel
  ) ## end div
}
