##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

MGseaInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    ## data set parameters
    shiny::selectInput(ns("contrast"), "Select contrast", choices = NULL)
  )
}

MGseaUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    boardHeader(title = "Multi-Omics GSEA",
                info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ##----------------------------------------------------------------
      shiny::tabPanel(
        "multiGSEA",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          bs_alert(HTML("<b>MultiGSEA</b> combines pathway enrichment on multiple omics layers to create a robust composite multi-omics pathway enrichment measure.")),
          bslib::layout_columns(
            height = "calc(100vh - 181px)",
            col_widths = bslib::breakpoints(
              sm = c(12, 12, 12, 12),
              xl = c(5, 3, 4, 7, 5),              
              xxxl = c(4, 3, 5, 4)
            ),
            mgsea_plot_enrichment_ui(
              ns("menrichment"),
              title = "multiGSEA enrichment",
              info.text = "Upon selection of a contrast from the 'Contrast' drop-down menu on the right, an ordered bar plot of cumulative absolute enrichment score (x-axis) of each gene set is shown. Please select a gene set from the 'multiGSEA scores' table to view a Volcano plot of Effect-size (log2FC) vs Significance (-log10q) for each feature mapped in the selected gene set.",
              caption = "Multi-omics gene set enrichment analysis.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            ),            
            mofa_plot_mgsea_ui(
              ns("mgsea_plot"),
              title = "MultiGSEA plot",
              info.text = "MultiGSEA plot. Scatter plot of enrichment scores of two omics data types. Pathway/genesets that are enriched in multiple datatypes are expected to exhibit higher multi.score values. To aid improved visualization a small random noise is added as standard normal values to the enrichment scores. Pathways enriched in both datatypesa are colored in red.",
              caption = "Scatter plot of enrichment scores of two omics data types. Pathway/genesets that are enriched in multiple datatypes are expected to exhibit higher multi.score values. To aid improved visualization a small random noise is added as standard normal values to the enrichment scores.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_plot_pathwayheatmap_ui(
              ns("pathwayheatmap"),
              title = "Pathway heatmap",
              info.text = "Heatmap of hierarchical clustering of z-score of features mapped in the selected gene set from the 'multiGSEA scores' table. Metadata variables present in the uploaded samples.csv are shown in the annotation bars at the top of the heatmap. A 'prefix' annotates each feature's datatype.",
              #info.methods = ".....",
              #info.references = list(list("X et al")),
              caption = "Integrated Multi-omics heatmap.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_table_mgsea_ui(
              ns("mgsea_table"),
              title = "multiGSEA scores",
              info.text = "Table of multi-omics genesets and associated statistics. Geneset enrichment is performed using rank-correlation. A raw multi-omics score is first computed as the absolute geometric mean across datatypes' genesets. The 'sign' of the gset rank correlation-based enrichment is also determined (-1 for negative, +1 for positive correlations) and row-wise standard deviation computed. A 'multi-omics p-value' is also computed using Stouffer method for p-value integration. Multi-omics p-value are then corrected for multiple testing errors to derive 'multi-omics q-values'. The listed 'multi.score' column is computed as an integrated score using the following formula: 'raw multi-omics score' * -log(multi-omics p-value) * 'multi-omics sign'. Pathways/genesets that are enriched in multiple datatypes are expected to have higher multi.score values.",
              caption = "Table of multi-omics genesets and associated statistics.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            ),
            mofa_plot_pathbank_ui(
              ns("pathbank_pathway"),
              title = "Multi-omics pathway",
              info.text = "Multi-omics analyses ultimately aim to overlay distinct datatypes into integrated cellular pathways where transripts, proteins and metabolites may exert biological effects. Please select a gene set from the 'multiGSEA scores' table to view the associated pathway if available.",
              caption = "Pathways that integrate proteomics and metabolomics data types in a single pathway diagram.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            )                        
          )
        )
      )

      
    )
  )
}
