##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

SNFInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    ## data set parameters
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
          )
        )
      )
    )
  )
}

SNFUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    boardHeader(title = "Similarity Network Fusion", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ##----------------------------------------------------------------
      shiny::tabPanel(
        "SNF Clustering",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          bs_alert(HTML("<b>Similarity Network Fusion</b> (SNF) is a network-based method for multi-omics integration by taking multiple views of a network and fusing them together to construct an overall status matrix (Wang et al., 2014). By integrating multiple datatypes, SNF improves the clustering of the samples.")),
          bslib::layout_columns(
            col_widths = bslib::breakpoints(
              xxxl = c(7, 5),
              xl = c(7, 5),              
              sm = c(12, 12)
            ),
            mofa_plot_snf_ui(
              ns("snf_affinity"),
              title = "SNF affinity matrices / t-SNE",
              info.text = "Similarity Network Fusion fuses multiple views (data types)  together to construct an overall integrated matrix. The heatmaps display sample correlation of pairwise Euclidean distances and the final integrated affinity matrix. The learned matrix can then be used for multiple analyses, including clustering, and classification.",
              info.methods = "Prior to SNF, missing values (if any) are imputed using SVD2. For each data type, pairwise Pearson correlation distances are computed between all pairs of data points. Affinity matrices are then calculated from these distance matrices, using number of neighbors K=10-30 and hyperparameter alpha=0.5. Similarity Network Fusion then fuses the matrices together to construct an overall integrated matrix.",
              info.references = list(list("Wang B, Mezlini A, Demir F, Fiume M, Zu T, Brudno M, Haibe-Kains B, Goldenberg A (2014). “Similarity Network Fusion: a fast and effective method to aggregate multiple data types on a genome wide scale.” Nature Methods.", "https://www.nature.com/articles/nmeth.2810")),
              caption = "Each datatype affinity matrix captures the pairwise similarities between samples, highlighting high similarities among samples within the same datatype.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            bslib::layout_columns(
              col_widths = 12,
              mofa_plot_snf_heatmap_ui(
                ns("snf_heatmap"),
                title = "SNF heatmap",
                info.text = "Clustering of SNF-integrated multi-omics data. Heatmap of normalized multi-omics data. The SNF clusters capture multi-omic features exhibiting similar behavior. Therefore, the heatmap is well versed to enable assessment of samples' clustering driven by multiple data types/modalities.",
                info.methods = "Heatmap of normalized multi-omics data. Pearson correlation is used as similarity measure. The clustering and number of clusters was determined using the Louvain method.",
                info.references = list(list("Wang B, Mezlini A, Demir F, Fiume M, Zu T, Brudno M, Haibe-Kains B, Goldenberg A (2014). “Similarity Network Fusion: a fast and effective method to aggregate multiple data types on a genome wide scale.” Nature Methods.", "https://www.nature.com/articles/nmeth.2810")),
                ##caption = "Clustering of SNF-integrated multi-omics data. Heatmap of normalized multi-omics data. Samples' clustering is driven by Similarity Network Fusion (SNF) method. Euclidean distance is used as similarity measure.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              mofa_plot_snfgraph_ui(
                ns("snf_cluster"),
                title = "Graph clustering of samples",
                info.text = "Multi-omics clustering of samples using SNF. The color of the edges correspond to the datatype specific correlation strength. The threshold an label type can be chosen in the plot options.",
                info.methods = "Edge colors correspond to datatypes.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              )
            )
          )
        )
      )

    )
  )
}
