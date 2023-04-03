##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ClusteringInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    withTooltip(shiny::selectInput(ns("hm_features"), "Features:", choices = NULL, multiple = FALSE),
      "Select a family of features.",
      placement = "top"
    ),
    withTooltip(
      shiny::radioButtons(ns("hm_clustmethod"), "Layout:",
        c("default", "tsne", "pca", "umap"),
        inline = TRUE
      ),
      "Choose the layout method for clustering to visualise.",
    ),
    shiny::conditionalPanel(
      "input.hm_features == '<custom>'",
      ns = ns,
      withTooltip(
        shiny::textAreaInput(ns("hm_customfeatures"), NULL,
          value = NULL,
          height = "150px", width = "100%",
          rows = 5, placeholder = "Paste your custom gene list"
        ),
        "Paste a custom list of genes to be used as features.",
        placement = "bottom"
      )
    ),
    shiny::conditionalPanel(
      "input.hm_features == '<contrast>'",
      ns = ns,
      tipifyR(
        shiny::selectInput(ns("hm_contrast"), NULL, choices = NULL),
        "Select contrast to be used as signature."
      )
    ),
    withTooltip(shiny::selectInput(ns("hm_group"), "Group by:", choices = NULL),
      "Group the samples by condition.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(
      shiny::selectInput(ns("hm_samplefilter"), "Filter samples:",
        choices = NULL, multiple = TRUE
      ),
      "Filter the relevant samples for the analysis.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(shiny::actionLink(ns("hm_options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(),
    shiny::conditionalPanel(
      "input.hm_options % 2 == 1",
      ns = ns,
      shiny::tagList(
        withTooltip(shiny::selectInput(ns("hm_level"), "Level:", choices = c("gene", "geneset")),
          "Specify the level analysis: gene or geneset level.",
          placement = "top", options = list(container = "body")
        ),
        withTooltip(shiny::checkboxInput(ns("hm_filterXY"), "exclude X/Y genes", FALSE),
          "Exclude genes on X/Y chromosomes.",
          placement = "top", options = list(container = "body")
        ),
        withTooltip(
          shiny::checkboxInput(
            ns("hm_filterMitoRibo"),
            "exclude mito/ribo genes", FALSE
          ),
          "Exclude mitochondrial (MT) and ribosomal protein (RPS/RPL) genes.",
          placement = "top", options = list(container = "body")
        )
      )
    )
  )
}

ClusteringUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 800 ## full height of page
  rowH  <- 350

  fullH <- "80vh" ## full height of full page
  rowH  <- "39.2vh"
  
  div(
    class = "row",
    ## h4("Cluster Samples"),
    boardHeader(title = "Cluster Samples", info_link = ns("board_info")),
    div(
      class = "col-md-7",
      shiny::tabsetPanel(
        id = ns("tabs1"),
        shiny::tabPanel(
          "Heatmap",
          clustering_plot_splitmap_ui(
            id = ns("splitmap"),
            label = "a",
            height = c(fullH,"80vh"),
            width = "100%"
          )
        ),
        shiny::tabPanel(
          "PCA/tSNE",
          clustering_plot_clustpca_ui(
            ns("PCAplot"),
            label = "",
            height = c(fullH, "70vh"),
            parent = ns
          ),
          tags$div(
            class = "caption",
            HTML("<b>PCA/tSNE plot.</b> The plot visualizes the similarity in expression of
                          samples as a scatterplot in reduced dimension (2D or 3D).
                          Samples that are similar are clustered near to each other, while samples with different
                          expression are positioned farther away. Groups of samples with similar profiles
                          will appear as <i>clusters</i> in the plot.")
          )
        ),
        shiny::tabPanel(
          "Parallel",
          shinyjqui::jqui_sortable(
              bslib::layout_column_wrap(
                 width = 1,                 
                 clustering_plot_parcoord_ui(
                     id = ns("parcoord"),
                     label = "a",
                     width = c("100%", "100%"),
                     height = c(rowH, 600)
                 ),
                 clustering_table_parcoord_ui(
                     id = ns("parcoord"),
                     label = "a",
                     width = c("100%", "100%"),
                     height = c(rowH, 600)
                 )
              ) ## layout   
          ) ## sortable
        )
    )),
    div(
      class = "col-md-5",
      shiny::tabsetPanel(
        id = ns("tabs2"),
        shiny::tabPanel(
          "Annotate clusters",
          clustering_plot_clusterannot_ui(
            id = ns("plots_clustannot"),
            label = "a",
            height = c(380, 650),
            width = c("100%", "100%")
          ),
          clustering_table_clustannot_ui(
            ns("tables_clustannot"),
            height = c(280, TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        ),
        shiny::tabPanel(
          "Phenotypes",
          clustering_plot_phenoplot_ui(
            id = ns("clust_phenoplot"),
            label = "",
            height = c(fullH, 650)
          )
        ),
        shiny::tabPanel(
          "Feature ranking",
          clustering_plot_featurerank_ui(
            id = ns("clust_featureRank"),
            label = "",
            height = c(fullH, 650),
            width = c("auto", "100%")
          )
        )
      )
    )
  )
}
