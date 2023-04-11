##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

EnrichmentInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(shiny::selectInput(ns("gs_contrast"), "Contrast:", choices = NULL),
      "Select a contrast of interest for the analysis.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("gs_features"), "Gene set collection:", choices = NULL, multiple = FALSE),
      "Choose a specific gene set collection for the analysis.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("gs_fdr"), "FDR", c(1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1), selected = 0.2),
      "Set the false discovery rate (FDR) threshold.",
      placement = "top"
    ),
    withTooltip(
      shiny::selectInput(ns("gs_lfc"), "logFC threshold",
        choices = c(0, 0.05, 0.1, 0.2, 0.5, 1, 2), selected = 0
      ),
      "Set the logarithmic fold change (logFC) threshold.",
      placement = "top"
    ),
    withTooltip(shiny::actionLink(ns("gs_options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(), br(),
    shiny::conditionalPanel(
      "input.gs_options % 2 == 1",
      ns = ns,
      shiny::tagList(
        withTooltip(shiny::checkboxInput(ns("gs_showall"), "Show all genesets", FALSE),
          "Enbale significant genes filtering. Display only significant genesets in the table.",
          placement = "top", options = list(container = "body")
        ),
        withTooltip(shiny::checkboxGroupInput(ns("gs_statmethod"), "Statistical methods:", choices = NULL),
          "Select a method or multiple methos for the statistical test.",
          placement = "right", options = list(container = "body")
        ),
        withTooltip(shiny::checkboxInput(ns("gs_top10"), "top 10 gene sets", FALSE),
          "Display only top 10 differentially enirhced gene sets (positively and negatively) in the <b>enrihcment analysis</b> table.",
          placement = "top", options = list(container = "body")
        ),
      )
    )
  )
}

EnrichmentUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 800
  rowH <- 420 ## row height of panels
  imgH <- 340 ## height of images
  tabV <- "70vh" ## height of tables
  tabH <- 340 ## row height of panels
  tabH <- "80vh" ## height of tables

  tabs <- tagList(
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "Top enriched",
        div(
          class = "row",
          div(
            class = "col-md-6",
            enrichment_plot_top_enrich_gsets_ui(
              ns("topEnriched"),
              title = "Top enriched gene sets",
              info.text = "This plot shows the top enriched gene sets for the selected comparison in the Contrast settings. Black vertical bars indicate the rank of genes in the gene set in the sorted list metric. The green curve corresponds to the 'running statistics' of the enrichment score (ES). The more the green ES curve is shifted to the upper left of the graph, the more the gene set is enriched in the first group. Conversely, a shift of the ES curve to the lower right, corresponds to more enrichment in the second group.",
              caption = "Gene set enrichment plots of the top differentially enriched gene sets. ",
              height = c("calc(45vh - 120px)", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          ),
          div(
            class = "col-md-6",
            enrichment_plot_freq_top_gsets_ui(
              ns("topEnrichedFreq"),
              title = "Frequency in top gene sets",
              info.text = "The plot shows the number of times a gene is present in the top-N genesets sorted by frequency. Genes that are frequently shared among the top enriched gene sets may suggest driver genes.",
              caption = "Gene frequency plot indicating the most recurring genes across the most correlated gene sets.",
              height = c("calc(45vh - 120px)", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),
      shiny::tabPanel(
        "Plots",
        div(
          class = "row",
          div(
            class = "col-md-3",
            enrichment_plot_volcano_ui(
              ns("subplot_volcano"),
              title = "Volcano plot",
              info.text = "Volcano-plot showing significance versus fold-change on the y and x axes, respectively. Genes in the gene set that is selected from the enrichment analysis Table I are highlighted in blue.",
              caption = "Volcano-plot showing significance versus fold-change with genes from the selected gene set highlighted.",
              height = c("calc(45vh - 120px)", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          ),
          div(
            class = "col-md-3",
            enrichment_plot_barplot_ui(
              ns("subplot_barplot"),
              title = "Enrichment barplot",
              info.text = "An enrichment barplot per sample group for the gene set that is selected from the enrichment analysis Table I. Samples can be ungrouped in the barplot by selecting ungroup samples from the plot Settings.",
              caption = "Barplot of the selected gene set in the phenotypic groups. ",
              height = c("calc(45vh - 120px)", TABLE_HEIGHT_MODAL),
              width = c("auto", 900)
            )
          ),
          div(
            class = "col-md-3",
            enrichment_plot_geneplot_ui(
              ns("subplot_geneplot"),
              title = "Expression geneplot",
              info.text = "An expression barplot per sample group for the gene that is selected from the genes Table II. Samples can be ungrouped in the barplot by selecting ungroup samples from the plot Settings.",
              caption = "Barplot of the selected gene in the phenotypic groups. ",
              height = c("calc(45vh - 120px)", TABLE_HEIGHT_MODAL),
              width = c("auto", 900)
            )
          ),
          div(
            class = "col-md-3",
            enrichment_plot_scatter_ui(
              ns("subplot_scatter"),
              title = "Enrichment vs. expression",
              info.text = "A scatter plot of enrichment scores versus expression values across the samples for the gene set selected from the enrichment analysis Table I and the gene selected from the genes Table II.",
              caption = "Scatter plot of the selected gene set enrichment scores versus the selected gene expression values by sample.",
              height = c("calc(45vh - 120px)", TABLE_HEIGHT_MODAL),
              width = c("auto", 900)
            )
          )
        )
      ),
      shiny::tabPanel(
        "Compare",
        bslib::layout_column_wrap(
           width = 1,
           enrichment_plot_compare_ui(
              ns("compare"),
              title = "Enrichment of geneset across multiple contrasts",
              info.text = "Under the Compare tab, enrichment profiles of the selected geneset in enrichment Table I can be visualised against all available contrasts.",
              caption = "Enrichment plots for the selected gene set (in Table I) across multiple contrasts.",
              height = c("calc(45vh - 120px)", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
           )
        )
      ),
      shiny::tabPanel(
        "Volcano (all)",
        bslib::layout_column_wrap(
           width = 1,
           enrichment_plot_volcanoall_ui(
             id = ns("volcanoAll"),
             title = "Volcano plots for all contrasts",
             info.text = "Under the Volcano (all) tab, the platform simultaneously displays multiple volcano plots for gene sets across all contrasts. This provides users an overview of the statistics across all comparisons. By comparing multiple volcano plots, the user can immediately see which comparison is statistically weak or strong based on the 'height' of the 'wings'.",
             caption = "Simultaneous visualisation of volcano plots of gene set enrichment across all contrasts.",
             height = c("calc(45vh - 120px)", TABLE_HEIGHT_MODAL),
             width = c("auto", "100%")
           )
         )
      ),
      shiny::tabPanel(
        "Volcano (methods)",
        bslib::layout_column_wrap(
           width = 1,
           enrichment_plot_volcanomethods_ui(
             ns("volcanoMethods"),
             title = "Volcano plots for all methods",
             info.text = "The Volcano (methods) panel displays the volcano plots provided by different enrichment calculation methods. This provides users an quick overview of the sensitivity of the statistical methods at once. Methods showing better statistical significance will show volcano plots with 'higher' wings.",
             caption = "Simultaneous visualisation of volcano plots of gene sets for different enrichment methods.",
             height = c("calc(45vh - 120px)", TABLE_HEIGHT_MODAL),
             width = c("auto", "100%")
           )
         )
      )
    ),
    shiny::tabsetPanel(
      id = ns("tabs2"),
      shiny::tabPanel(
        "Table",
        div(
          class = "row",
          div(
            class = "col-md-7",
            enrichment_table_enrichment_analysis_ui(
              ns("gseatable"),
              title = "Enrichment analysis",
              info.text = "Similar to the differential gene expression analysis, users can perform differential expression analysis on a geneset level that is referred as gene set enrichment analysis. To ensure statistical reliability, the platform performs the gene set enrichment analysis using multiple methods, including Spearman rank correlation , GSVA , ssGSEA , Fisher exact test , GSEA , camera and fry . The combined result from the methods is displayed in this table, where for each geneset the meta.q corresponds to the highest q val. The number of stars indicates how many methods detected a significant correlation.",
              caption = "Table summarizing the statistical results of the gene set enrichment analysis for selected contrast. ",
              height = c("calc(45vh - 120px)", TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            )
          ),
          div(
            class = "col-md-5",
            enrichment_table_genes_in_geneset_ui(
              ns("genetable"),
              title = "Genes in gene set",
              info.text = "By clicking on a gene set in the table I, it is possible to see the gene list of that gene set in this table. By clicking on a gene in this table, users can check the expression status of the gene for the selected contrast in the Expression barplot and its correlation to the gene set in the Gene to gene set correlation scatter plot under the Plots section.",
              caption = "Table showing the fold-change, statistics and correlation of the genes in the selected gene set.",
              height = c("calc(45vh - 120px)", TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            )
          )
        )
      ),
      shiny::tabPanel(
        "Enrichment (all)",
        bslib::layout_column_wrap(
          width = 1,
          enrichment_table_gset_enrich_all_contrasts_ui(
            ns("fctable"),
            title = "Gene set enrichment for all contrasts",
            info.text = "The column `fc.var` corresponds to the variance of the gene set across all contrasts.",
            caption = "The Enrichment (all) panel reports the gene set enrichment for all contrasts in the selected dataset.",
            height = c("calc(45vh - 120px)", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          )
        )
      ),
      shiny::tabPanel(
        "FDR table",
        bslib::layout_column_wrap(
          width = 1,
          enrichment_table_n_sig_gsets_ui(
            ns("FDRtable"),
            title = "Number of significant gene sets",
            info.text = "Using the table the user can determine which statistical methods perform better for a particular contrast.",
            caption = "The FDR table panel reports the number of significant gene sets at different FDR thresholds, for all contrasts and all methods. ",
            height = c("calc(45vh - 120px)", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          )
        )
      )
    )
  )
  div(
    boardHeader(title = "Geneset enrichment", info_link = ns("gs_info")),
    tabs
  )
}
