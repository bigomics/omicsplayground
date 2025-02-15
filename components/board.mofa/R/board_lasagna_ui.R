##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

LasagnaInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    ## data set parameters
    shiny::selectInput(ns("contrast"), "Select comparison", choices = NULL),
    shiny::br(),
    bslib::accordion(
      id = ns("data_type_accordion"),
      open = FALSE,
      bslib::accordion_panel(
        "Options",
        icon = icon("cog", lib = "glyphicon"),
        shiny::tagList(
          shiny::selectInput(ns("ngenes"), tspan("Number genes:"),
            choices = c(500, 1000, 2000, 4000, 8000),
            selected = 1000
          ),
          shiny::selectInput(ns("cutheight"), "Merge cut height",
            choices = c(0.05, 0.10, 0.25, 0.5, 0.9, 0.999),
            selected = 0.25
          )
        )
      )
    )
  )
}

my_navset_card_tab <- function(...) {
  htmltools::tagAppendAttributes(
    bslib::navset_card_tab(
      ... ,
      tags$style(HTML("@media (min-width: 1200px) {.root_navset { height: calc(100vh - 36px); }}"))
    ),
    class = "root_navset border-0"
  )
}


LasagnaUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

##  shiny::div(
  ##boardHeader(title = "LASAGNA", info_link = ns("info")),
  ##  shiny::tabsetPanel(
  #bslib::navset_card_tab(     
  my_navset_card_tab(
    id = ns("tabs"),
    #title = boardHeader(title = "LASAGNA", info_link = ns("info")),
    title = "LASAGNA",
    height = "calc(100vh - 50px)",
    ##----------------------------------------------------------------
    ##shiny::tabPanel(
    bslib::nav_panel(      
      "Multi-layer model",
      bslib::layout_columns(
        col_widths = 12,
        row_heights = c("auto",1),
        bs_alert(HTML("<b>LASAGNA</b> is a stacked layer model for multi-omics where each layer corresponds to a data type. The acronym stands for a <u>L</u>ayered <u>A</u>pproach to <u>S</u>imultaneous <u>A</u>nalysis of <u>G</u>enomic and <u>N</u>etwork <u>A</u>ssociations'.")),
        bslib::layout_columns(
          col_widths = c(6,6),
          mofa_plot_lasagna_ui(
            ns("lasagna"),
            title = "Multi-layer model",
            info.text = "Layered Approach to Simultaneous Analysis of Genomic and Network Associations ('LASAGNA'). The LASAGNA plot is a stacked layer plot to visualize multi-omics data. Specifically, each layer shows a data type-specific UMAP. LASAGNA just shows the datatype-specific UMAPs overlayed.",
            caption = "Layered Approach to Simultaneous Analysis of Genomic and Network Associations ('LASAGNA'). The LASAGNA plot is a stacked layer plot to visualize data type-specific UMAPs overlayed.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          ),
          mofa_plot_clustering_ui(
            ns("clusters"),
            title = "Feature UMAP per datatype",
            info.text = "Feature-level UMAP clustering per data type. Visually explore signatures of distinct datatypes across the same set of samples. Feature-level clustering is determined by Uniford Manifold Approximation and Projection (UMAP) applied to each data type separately. UMAP is computed using the uwot R package. Feature-level clustering enables assessment of how the distinct data types/modalities define distinct (functional) groups. This analysis may reveal that distinct data types capture different heterogeneities in the data, potentially associated with unique biological functions. On the other end, similar clustering patterns between distinct data types may indicate shared regulation. The colors in the UMAP reflect the low-to-high correlation with the selected comparison to explore the impact of different conditions.",
            info.references = list(list("Melville J (2024). “uwot: The Uniform Manifold Approximation and Projection (UMAP) Method for Dimensionality Reduction.”, CRAN.", "https://doi.org/10.32614/CRAN.package.uwot")),
            caption = "Feature-level UMAP clustering per data type. Visually explore signatures of distinct datatypes across the same set of samples. Feature-level clustering is determined by UMAP applied to each data type separately. The colors in the UMAP reflect the low-to-high correlation with the selected comparison ",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    ),

      
    ##----------------------------------------------------------------
    ## shiny::tabPanel(
    bslib::nav_panel(      
      "Path scoring",
      bslib::layout_columns(
        col_widths = 12,
        #height = "calc(100vh - 181px)",
        ## bs_alert(HTML("<b>LASAGNA</b> is a stacked layer model for multi-omics where each layer corresponds to a data type. The acronym stands for 'a Layered Approach to Simultaneous Analysis of Genomic and Network Associations'.")),
        bslib::layout_columns(
          col_widths = bslib::breakpoints(
            xxxl = c(6, 6),
            xl = c(12, 12),              
            sm = c(12, 12)
          ),
          bslib::layout_columns(
            col_widths = bslib::breakpoints(
              xxxl = c(6, 6, 12),
              xl = c(5, 3, 4),              
              sm = c(12, 12, 12)
            ),
            mofa_plot_lasagnaSP_ui(
              ns("lasagnaSP"),
                title = "Parallel coordinates",
              info.text = "LASAGNA is a acronym of 'Layered Approach to Simultaneous Analysis of Genomic and Network Associations.'",
              caption = "Scale independence and mean connectivity plots to determine the soft threshold.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_plot_lasagnaSP_scores_ui(
              ns("lasagnaSP"),
              title = "Path scores",
              info.text = "Path scores",
              caption = ""
            ),
            mofa_plot_lasagnaSP_frequency_ui(
              ns("lasagnaSP"),
              title = "Term frequency",
              info.text = "Frequency of features",
              caption = "Term frequency in top ranked shortest path solutions."
            )
          ),
          mofa_table_lasagnaSP_ui(
            ns("lasagnaSP"),
            title = "Pathscore table",
            info.text = "Clustering of features",
            caption = "Each datatype affinity matrix captures the pairwise similarities between samples, highlighting high similarities among samples within the same datatype.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    )
    
  )
  ##)
}
