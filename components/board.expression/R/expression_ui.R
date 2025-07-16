##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ExpressionInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    withTooltip(shiny::selectInput(ns("gx_contrast"), "Contrast:", choices = NULL),
      "Select a contrast of interest for the analysis.",
      placement = "top"
    ),
    withTooltip(
      shiny::selectInput(ns("gx_features"), tspan("Gene family:"),
        choices = NULL, multiple = FALSE
      ),
      "Choose a specific gene family for the analysis.",
      placement = "top"
    ),
    bslib::layout_column_wrap(
      width = 1/2,
      withTooltip(
        selectInput(ns("gx_fdr"),
          "FDR",
          choices = c(1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1),
          selected = 0.2
        ),
        "Set the false discovery rate (FDR) or P-value threshold.",
        placement = "top"
      ),
      withTooltip(
        shiny::selectInput(ns("gx_lfc"), "logFC",
          choices = c(0, 0.1, 0.2, 0.5, 1, 2, 5), selected = 0
        ),
        "Set the logarithmic fold change (logFC) threshold.",
        placement = "top"
      )
    ),
    shiny::br(),
    bslib::accordion(
      id = ns("gx_accordion"),
      open = FALSE,
      bslib::accordion_panel(
        "Options",
        icon = icon("cog", lib = "glyphicon"),
        withTooltip(shiny::checkboxInput(ns("gx_showall"), tspan("Show all genes"), FALSE),
          "Display all genes in the table. Disable filtering of significant genes.",
          placement = "top", options = list(container = "body")
        ),
        withTooltip(shiny::checkboxInput(ns("show_pv"), "Show p-values", FALSE),
          "Show p-values in the table. WARNING: Nominal p-values are NOT corrected for multiple testing errors. We do not advice their use.",
          placement = "top", options = list(container = "body")
        ),
        br(),
        withTooltip(
          shiny::checkboxGroupInput(ns("gx_statmethod"), "Statistical methods:",
            choices = NULL, inline = TRUE
          ),
          "Select a method for the statistical test. To increase the statistical reliability of the Omics Playground,
            we perform the DE analysis using commonly accepted methods in the literature, including t-test (standard,
            Welch), limma (no trend, trend, voom), edgeR (QLF, LRT), and DESeq2 (Wald, LRT), and merge the results.",
          placement = "right", options = list(container = "body")
        ),
        withTooltip(
          shiny::selectInput(
            inputId = ns("pval_cap"),
            label = "Significance cap",
            choices = c("1e-12", "1e-20", "Uncapped")
          ),
          "Significance cap",
          placement = "right", options = list(container = "body")
        )
      )
    )
  )
}

ExpressionUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- "calc(100vh - 125px)" ## full height of page (minus header)
  halfH <- "calc(50vh - 125px)" ## half height of page

  tabs1 <- shiny::tabsetPanel(
    id = ns("tabs1"),
    shiny::tabPanel(
      "Overview",
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        height = halfH,
        expression_plot_volcano_ui(ns("plots_volcano"),
          label = "a",
          title = "Volcano plot",
          info.text = "Volcano plot of genes for the selected {Contrast} displaying fold-change versus significance. By selecting a specific gene under the Differential expression analysis table it will be highlighted. Imputed genes are displayed with a cross. All other genes are displayed with a filled circle. Similarly, if a geneset is selected under the Gene sets with gene table it will be highlighted.",
          info.methods = "Statistical significance assessed using three independent statistical methods: DESeq2 (Wald test) [1], edgeR (QLF test) [2] and limma-trend [3]. The maximum q-value of the three methods is taken as aggregate q-value, which corresponds to taking the intersection of significant genes from all three tests.",
          info.references = list(
            list(
              "Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550.",
              "https://doi.org/10.1186/s13059-014-0550-8"
            ),
            list(
              "Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), 139-140.",
              "https://doi.org/10.1093/bioinformatics/btp616"
            ),
            list(
              "Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7)",
              "https://doi.org/10.1093/nar/gkv007"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#statistical-testing",
          caption = "Volcano-plot displaying fold-change versus significance.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),
        expression_plot_maplot_ui(
          id = ns("plots_maplot"),
          title = "MA plot",
          info.text = "MA plot of genes for the selected {Contrast} displaying fold-change (M-values) versus the mean intensity (A-values). Imputed genes are displayed with a cross. All other genes are displayed with a filled circle. By selecting a specific gene under the Differential expression analysis table it will be highlighted. Similarly, if a geneset is selected under the Gene sets with gene table it will be highlighted.",
          info.methods = "See Volcano plot",
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#statistical-testing",
          caption = "MA-plot displaying signal intensity versus fold-change.",
          label = "b",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),
        expression_plot_barplot_ui(
          id = ns("plots_barplot"),
          title = "Differential expression",
          info.text = "Barplot displaying the expression for the gene selected under the Differential expression analysis table for the conditions of the selected {Contrast}. Under the plot settings it is possible to control the scale with {log scale}, show the expression of other groups {show others} and group/ungroup the conditions {grouped}.",
          caption = "Barplot displaying the expression for the gene selected under the Differential expression analysis table for the conditions of the selected {Contrast}.",
          label = "c",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),
        expression_plot_topfoldchange_ui(
          id = ns("plots_topfoldchange"),
          title = "Fold change by comparison",
          info.text = "Barplot displaying the fold change summary across all contrasts for the gene selected in the Differential expression analysis table.",
          caption = "Sorted barplot of the differential expression of the selected gene across all contrasts.",
          label = "d",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ),
    shiny::tabPanel(
      "Top features",
      bslib::layout_columns(
        col_widths = 12,
        height = halfH,
        expression_plot_topgenes_ui(
          id = ns("topgenes"),
          title = "Expression of top differentially expressed features",
          info.text = "Barplot displaying the expression for the top differentially expressed features for the conditions of the selected {Contrast}. Under the plot settings it is possible to control the scale with {log scale}, show the expression of other groups {show others} and group/ungroup the conditions {grouped}.",
          caption = "Expression barplots of the top most differentially expressed genes for the selected comparison.",
          label = "a",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ),
    shiny::tabPanel(
      "Volcano by comparison",
      bslib::layout_columns(
        col_widths = 12,
        height = halfH,
        expression_plot_volcanoAll_ui(
          id = ns("volcanoAll"),
          title = "Volcano plots for all contrasts",
          info.text = "Volcano plot of genes for all contrasts displaying fold-change versus significance. The plots can be colored by using the {Color up/down regulated} plot setting; also the plots can be scaled using the {scale per plot} plot setting.",
          info.methods = "Statistical significance assessed using three independent statistical methods: DESeq2 (Wald test) [1], edgeR (QLF test) [2] and limma-trend [3]. The maximum q-value of the three methods is taken as aggregate q-value, which corresponds to taking the intersection of significant genes from all three tests. By comparing multiple volcano plots, it can immediately be seen which comparison is statistically weak or strong.",
          info.references = list(
            list(
              "Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550.",
              "https://doi.org/10.1186/s13059-014-0550-8"
            ),
            list(
              "Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), 139-140.",
              "https://doi.org/10.1093/bioinformatics/btp616"
            ),
            list(
              "Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7)",
              "https://doi.org/10.1093/nar/gkv007"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#statistical-testing",
          caption = "Simultaneous visualisation of volcano plots of genes for all comparisons.",
          label = "a",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ),
    shiny::tabPanel(
      "Volcano by method",
      bslib::layout_columns(
        col_widths = 12,
        height = halfH,
        expression_plot_volcanoMethods_ui(
          id = ns("volcanoMethods"),
          title = "Volcano plots for all methods",
          info.text = "Volcano plots of genes for selected {Contrast} of all the statistical methods displaying fold-change versus significance. The plots can be colored by using the {Color up/down regulated} plot setting; also the plots can be scaled using the {scale per plot} plot setting.",
          info.methods = "Statistical significance assessed using three independent statistical methods: DESeq2 (Wald test) [1], edgeR (QLF test) [2] and limma-trend [3]. The maximum q-value of the three methods is taken as aggregate q-value, which corresponds to taking the intersection of significant genes from all three tests. By comparing multiple volcano plots, it can immediately be seen which method captures better signal. Methods showing better statistical significance will show volcano plots with 'higher wings'.",
          info.references = list(
            list(
              "Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550.",
              "https://doi.org/10.1186/s13059-014-0550-8"
            ),
            list(
              "Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), 139-140.",
              "https://doi.org/10.1093/bioinformatics/btp616"
            ),
            list(
              "Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7)",
              "https://doi.org/10.1093/nar/gkv007"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#statistical-testing",
          caption = "Simultaneous visualisation of volcano plots of genes by multiple differential expression methods for the selected contrast. ",
          label = "a",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ), ## end upper tabPanel
    shiny::tabPanel(
      "FC-FC comparison",
      bslib::layout_columns(
        col_widths = 12,
        height = halfH,
        expression_plot_fc_fc_ui(
          id = ns("fc_fc"),
          title = "FC-FC comparison",
          info.text = "Compare custom FC with baseline FC",
          caption = "FC-FC comparison: Compare custom FC with baseline FC",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ) ## end upper tabPanel
  ) ## end of tabs1

  tabs2 <- shiny::tabsetPanel(
    id = ns("tabs2"),
    shiny::tabPanel(
      "Table",
      bslib::layout_columns(
        col_widths = c(7, 5),
        height = halfH,
        expression_table_genetable_ui(
          ns("genetable"),
          title = "Differential expression analysis",
          info.text = "The table shows the results of the statistical tests. Omics Playground performs DE analysis using four commonly accepted methods, namely, T-test (standard, Welch), limma (no trend, trend, voom), edgeR (QLF, LRT), and DESeq2 (Wald, LRT), and combines the statistical results using a meta.q value that represents the highest q value among the methods. The number of stars indicate how many methods identified significant. The table is interactive (scrollable, clickable); users can sort by logFC, meta.q, or average expression in either conditions. The column pct.missingness reports the percentage of samples (within the selected contrast) in which a feature is missing in the original counts/abundance data matrix. A logFC value of 999 indicates NA (i.e., FC could not be computed due to feature missingness).",
          caption = "Table showing the significant results of the differential expression analysis on the selected contrast.",
          width = c("100%", "100%"),
          height = c("100%", TABLE_HEIGHT_MODAL)
        ),
        expression_table_gsettable_ui(
          ns("gsettable"),
          title = "Gene sets with gene",
          info.text = "By clicking on a gene in the Table I, it is possible to see which genesets contain that gene in this table, and check the differential expression status in other comparisons from the Gene in contrasts plot under the Plots tab.",
          caption = "Table indicating all the gene sets that contain the gene highlighted in the differential expression table.",
          width = c("100%", "100%"),
          height = c("100%", TABLE_HEIGHT_MODAL)
        )
      ) ## end of layout
    ),
    shiny::tabPanel(
      "Foldchange (all)",
      bslib::layout_columns(
        col_widths = 12,
        height = halfH,
        expression_table_fctable_ui(
          ns("fctable"),
          title = "Gene fold changes for all contrasts",
          info.text = "Differential expression (fold-change) across all contrasts. The column `rms.FC` corresponds to the root-mean-square fold-change across all contrasts.",
          caption = "The Foldchange (all) table reports the gene fold changes for all contrasts in the selected dataset.",
          width = c("100%", "100%"),
          height = c("100%", TABLE_HEIGHT_MODAL)
        )
      )
    ),
    shiny::tabPanel(
      "FDR table",
      bslib::layout_columns(
        col_widths = 12,
        height = halfH,
        expression_table_FDRtable_ui(
          ns("FDRtable"),
          title = "Number of significant genes",
          info.text = "Use this table to quickly see which methods are more sensitive. The left part of the table (in blue) correspond to the number of significant down-regulated genes, the right part (in red) correspond to the number of significant overexpressed genes.",
          caption = "This table reports the number of significant genes at different FDR thresholds for all contrasts and methods.",
          width = c("100%", "100%"),
          height = c("100%", TABLE_HEIGHT_MODAL)
        )
      )
    )
  ) ## end tabs2 = bottom tabsetPanel

  div(
    boardHeader(title = "Differential expression", info_link = ns("gx_info")),
    bslib::layout_columns(
      col_widths = 12,
      height = fullH,
      gap = "0px",
      tabs1,
      tabs2
    )
  )
}
