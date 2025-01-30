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
    withTooltip(shiny::selectInput(ns("gs_features"), tspan("Gene set collection:"), choices = NULL, multiple = FALSE),
      "Choose a specific gene set collection for the analysis.",
      placement = "top"
    ),
    shiny::fillRow(
      flex = c(1, 1),
      withTooltip(
        shiny::selectInput(ns("gs_fdr"), "FDR", c(1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1), selected = 0.2),
        "Set the false discovery rate (FDR) threshold.",
        placement = "top"
      ),
      withTooltip(
        shiny::selectInput(ns("gs_lfc"), "logFC",
          choices = c(0, 0.05, 0.1, 0.2, 0.5, 1, 2), selected = 0
        ),
        "Set the logarithmic fold change (logFC) threshold.",
        placement = "top"
      )
    ),
    shiny::br(), shiny::br(),
    shiny::br(), shiny::br(),
    withTooltip(shiny::actionLink(ns("gs_options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(), br(),
    ## shiny::conditionalPanel(
    ##   "input.gs_options % 2 == 1",
    ##   ns = ns,
    shiny::tagList(
      withTooltip(shiny::checkboxInput(ns("gs_showall"), tspan("Show all genesets"), FALSE),
        "Enbale significant genes filtering. Display only significant genesets in the table.",
        placement = "top", options = list(container = "body")
      ),
      withTooltip(shiny::checkboxGroupInput(ns("gs_statmethod"), "Statistical methods:", choices = NULL),
        "Select a method or multiple methos for the statistical test.",
        placement = "right", options = list(container = "body")
      ),
      withTooltip(shiny::checkboxInput(ns("gs_top10"), tspan("top 10 gene sets"), FALSE),
        "Display only top 10 differentially enirhced gene sets (positively and negatively) in the enrihcment analysis table.",
        placement = "top", options = list(container = "body")
      )
    )
    ## ),
  )
}

EnrichmentUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- "calc(100vh - 120px)" ## full height of page (minus header)
  halfH <- "calc(50vh - 120px)" ## half height of page

  tabs1 <- shiny::tabsetPanel(
    id = ns("tabs1"),
    shiny::tabPanel(
      "Enrichment",
      bslib::layout_columns(
        col_widths = c(5, 4, 3),
        height = halfH,
        enrichment_plot_top_enrich_gsets_ui(
          ns("topEnriched"),
          title = "Top enriched gene sets",
          info.text = "Plot showing the top enriched genesets for the selected {Contrast}. If nothing is selected on the Enrichment analysis table, the top genesets plots will be shown, while if a geneset is selected on the table, it will be displayed alone.",
          info.methods = "Geneset enrichment is performed using CAMERA [1], GSEA [2], ssGSEA [3], fGSEA [4], GSVA [5] and fry [6]. The q-values yielded by the different methods are then combined into a meta-q value, where the meta value corresponds to the maximum. Black vertical bars indicate the rank of genes in the gene set in the sorted list metric. The green curve corresponds to the 'running statistics' of the enrichment score (ES). The more the green ES curve is shifted to the upper left of the graph, the more the gene set is enriched in the first group. Conversely, a shift of the ES curve to the lower right, corresponds to more enrichment in the second group.",
          info.references = list(
            list(
              "Wu, D., & Smyth, G. K. (2012). Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic acids research, 40(17), e133-e133.",
              "https://doi.org/10.1093/nar/gks461"
            ),
            list(
              "Mootha, V. K., Lindgren, C. M., Eriksson, K. F., Subramanian, A., Sihag, S., Lehar, J., ... & Groop, L. C. (2003). PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. Nature genetics, 34(3), 267-273.",
              "https://doi.org/10.1038/ng1180"
            ),
            list(
              "Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Moody, S. E., Dunn, I. F., ... & Hahn, W. C. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature, 462(7269), 108-112.",
              "https://doi.org/10.1038/nature08460"
            ),
            list(
              "Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N., & Sergushichev, A. (2016). Fast gene set enrichment analysis. biorxiv, 060012.",
              "https://doi.org/10.1101/060012"
            ),
            list(
              "Hänzelmann, S., Castelo, R., & Guinney, J. (2013). GSVA: gene set variation analysis for microarray and RNA-seq data. BMC bioinformatics, 14, 1-15.",
              "https://doi.org/10.1186/1471-2105-14-7"
            ),
            list(
              "Ritchie, M. E., Phipson, B., Wu, D. I., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic acids research, 43(7), e47-e47.",
              "https://doi.org/10.1093/nar/gkv007"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#functional-analyses",
          caption = "Gene set enrichment plots of the top differentially enriched gene sets. ",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),
        enrichment_plot_volcano_ui(
          ns("subplot_volcano"),
          title = "Volcano plot",
          info.text = "Volcano-plot showing significance versus fold-change for the selected {Contrast}. By selecting a geneset on the Enrichment analysis table, the genes from it will be highlighted.",
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
              "Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7).",
              "https://doi.org/10.1093/nar/gkv007"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#statistical-testing",
          caption = "Volcano-plot showing significance versus fold-change with genes from the selected gene set highlighted.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),
        enrichment_plot_barplot_ui(
          ns("subplot_barplot"),
          title = "Enrichment barplot",
          info.text = "Enrichment barplot of selected {Contrast} for the gene set that is selected on the Enrichment analysis table. Samples can be ungrouped using the {ungroup samples} option in plot settings.",
          info.methods = "Geneset enrichment is performed using CAMERA [1], GSEA [2], ssGSEA [3], fGSEA [4], GSVA [5] and fry [6]. The q-values yielded by the different methods are then combined into a meta-q value, where the meta value corresponds to the maximum. Black vertical bars indicate the rank of genes in the gene set in the sorted list metric. The green curve corresponds to the 'running statistics' of the enrichment score (ES). The more the green ES curve is shifted to the upper left of the graph, the more the gene set is enriched in the first group. Conversely, a shift of the ES curve to the lower right, corresponds to more enrichment in the second group.",
          info.references = list(
            list(
              "Wu, D., & Smyth, G. K. (2012). Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic acids research, 40(17), e133-e133.",
              "https://doi.org/10.1093/nar/gks461"
            ),
            list(
              "Mootha, V. K., Lindgren, C. M., Eriksson, K. F., Subramanian, A., Sihag, S., Lehar, J., ... & Groop, L. C. (2003). PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. Nature genetics, 34(3), 267-273.",
              "https://doi.org/10.1038/ng1180"
            ),
            list(
              "Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Moody, S. E., Dunn, I. F., ... & Hahn, W. C. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature, 462(7269), 108-112.",
              "https://doi.org/10.1038/nature08460"
            ),
            list(
              "Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N., & Sergushichev, A. (2016). Fast gene set enrichment analysis. biorxiv, 060012.",
              "https://doi.org/10.1101/060012"
            ),
            list(
              "Hänzelmann, S., Castelo, R., & Guinney, J. (2013). GSVA: gene set variation analysis for microarray and RNA-seq data. BMC bioinformatics, 14, 1-15.",
              "https://doi.org/10.1186/1471-2105-14-7"
            ),
            list(
              "Ritchie, M. E., Phipson, B., Wu, D. I., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic acids research, 43(7), e47-e47.",
              "https://doi.org/10.1093/nar/gkv007"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#functional-analyses",
          caption = "Barplot of the selected gene set in the phenotypic groups. ",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", 900)
        )
      )
    ),
    shiny::tabPanel(
      tspan("Gene expression"),
      bslib::layout_columns(
        col_widths = c(3, 4, 5),
        height = halfH,
        enrichment_plot_geneplot_ui(
          ns("subplot_geneplot"),
          title = "Expression plot",
          info.text = "An expression barplot per sample group for the gene that is selected from the genes Table II. Samples can be ungrouped in the barplot by selecting ungroup samples from the plot Settings.",
          caption = "Barplot of the selected gene in the phenotypic groups. ",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", 900)
        ),
        enrichment_plot_scatter_ui(
          ns("subplot_scatter"),
          title = "Enrichment vs. expression",
          info.text = "Scatter plot of enrichment scores versus expression values for selected {Comparison} for the gene set selected from the Enrichment analysis table and the gene selected from the genes table.",
          info.methods = "See Enrichment barplot",
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#functional-analyses",
          caption = "Scatter plot of the selected gene set enrichment scores versus the selected gene expression values by sample.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", 900)
        ),
        enrichment_plot_freq_top_gsets_ui(
          ns("topEnrichedFreq"),
          title = "Most frequent genes",
          info.text = "Barchart showing the number of times a gene is present in the top-N genesets for the selected {Contrast} sorted by frequency. Genes that are frequently shared among the top enriched gene sets may suggest driver genes.",
          info.methods = "See Top enriched gene sets",
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#functional-analyses",
          caption = "Gene frequency plot indicating the most recurring genes across the most correlated gene sets.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ),
    shiny::tabPanel(
      "Enrichment by comparison",
      bslib::layout_columns(
        col_widths = 12,
        height = halfH,
        enrichment_plot_compare_ui(
          ns("compare"),
          title = "Enrichment of geneset across multiple comparisons",
          info.text = "Plot showing the the selected geneset (on Enrichment analysis table) for all available contrasts.",
          info.methods = "Geneset enrichment is performed using CAMERA [1], GSEA [2], ssGSEA [3], fGSEA [4], GSVA [5] and fry [6]. The q-values yielded by the different methods are then combined into a meta-q value, where the meta value corresponds to the maximum. Black vertical bars indicate the rank of genes in the gene set in the sorted list metric. The green curve corresponds to the 'running statistics' of the enrichment score (ES). The more the green ES curve is shifted to the upper left of the graph, the more the gene set is enriched in the first group. Conversely, a shift of the ES curve to the lower right, corresponds to more enrichment in the second group.",
          info.references = list(
            list(
              "Wu, D., & Smyth, G. K. (2012). Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic acids research, 40(17), e133-e133.",
              "https://doi.org/10.1093/nar/gks461"
            ),
            list(
              "Mootha, V. K., Lindgren, C. M., Eriksson, K. F., Subramanian, A., Sihag, S., Lehar, J., ... & Groop, L. C. (2003). PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. Nature genetics, 34(3), 267-273.",
              "https://doi.org/10.1038/ng1180"
            ),
            list(
              "Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Moody, S. E., Dunn, I. F., ... & Hahn, W. C. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature, 462(7269), 108-112.",
              "https://doi.org/10.1038/nature08460"
            ),
            list(
              "Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N., & Sergushichev, A. (2016). Fast gene set enrichment analysis. biorxiv, 060012.",
              "https://doi.org/10.1101/060012"
            ),
            list(
              "Hänzelmann, S., Castelo, R., & Guinney, J. (2013). GSVA: gene set variation analysis for microarray and RNA-seq data. BMC bioinformatics, 14, 1-15.",
              "https://doi.org/10.1186/1471-2105-14-7"
            ),
            list(
              "Ritchie, M. E., Phipson, B., Wu, D. I., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic acids research, 43(7), e47-e47.",
              "https://doi.org/10.1093/nar/gkv007"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#functional-analyses",
          caption = "Enrichment plots for the selected gene set (in Table I) across multiple comparisons.",
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
        enrichment_plot_volcanoall_ui(
          id = ns("volcanoAll"),
          title = "Volcano plots for all contrasts",
          info.text = "Volcano plot of genesets for all contrasts displaying fold-change versus significance. The plots can be scaled using the {scale per method} plot setting.",
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
              "Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7).",
              "https://doi.org/10.1093/nar/gkv007"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#statistical-testing",
          caption = "Simultaneous visualisation of volcano plots of gene set enrichment across all contrasts.",
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
        enrichment_plot_volcanomethods_ui(
          ns("volcanoMethods"),
          title = "Volcano plots for all methods",
          info.text = "Volcano plot of genesets for the selected {Contrast} displaying fold-change versus significance. The plots can be scaled using the {scale per method} plot setting.",
          info.methods = "Statistical significance assessed using multiple statistical methods: DESeq2 (Wald test) [1], edgeR (QLF test) [2] and limma-trend [3]. By comparing multiple volcano plots, it can immediately be seen which method is statistically weak or strong.",
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
              "Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7).",
              "https://doi.org/10.1093/nar/gkv007"
            )
          ),
          info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#statistical-testing",
          caption = "Simultaneous visualisation of volcano plots of gene sets for different enrichment methods.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    )
  )

  tabs2 <- shiny::tabsetPanel(
    id = ns("tabs2"),
    shiny::tabPanel(
      "Table",
      bslib::layout_columns(
        col_widths = c(8, 4),
        height = halfH,
        enrichment_table_enrichment_analysis_ui(
          ns("gseatable"),
          title = "Enrichment analysis",
          info.text = "Similar to the differential gene expression analysis, users can perform differential expression analysis on a geneset level that is referred as gene set enrichment analysis. To ensure statistical reliability, the platform performs the gene set enrichment analysis using multiple methods, including Spearman rank correlation , GSVA , ssGSEA , Fisher exact test , GSEA , camera and fry . The combined result from the methods is displayed in this table, where for each geneset the meta.q corresponds to the highest q val. The number of stars indicates how many methods detected a significant correlation.",
          caption = "Table summarizing the statistical results of the gene set enrichment analysis for selected contrast. ",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("100%", "100%")
        ),
        enrichment_table_genes_in_geneset_ui(
          ns("genetable"),
          title = "Genes in gene set",
          info.text = "By clicking on a gene set in the table I, it is possible to see the gene list of that gene set in this table. By clicking on a gene in this table, users can check the expression status of the gene for the selected contrast in the Expression barplot and its correlation to the gene set in the Gene to gene set correlation scatter plot under the Plots section.",
          caption = "Table showing the fold-change, statistics and correlation of the genes in the selected gene set.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("100%", "100%")
        )
      )
    ),
    shiny::tabPanel(
      "Enrichment (all)",
      bslib::layout_columns(
        col_widths = 12,
        height = halfH,
        enrichment_table_gset_enrich_all_contrasts_ui(
          ns("fctable"),
          title = "Gene set enrichment for all contrasts",
          info.text = "The column `fc.var` corresponds to the variance of the gene set across all contrasts.",
          caption = "The Enrichment (all) panel reports the gene set enrichment for all contrasts in the selected dataset.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("100%", "100%")
        )
      )
    ),
    shiny::tabPanel(
      "FDR table",
      bslib::layout_columns(
        col_widths = 12,
        height = halfH,
        enrichment_table_n_sig_gsets_ui(
          ns("FDRtable"),
          title = "Number of significant gene sets",
          info.text = "Using the table the user can determine which statistical methods perform better for a particular contrast.",
          caption = "The FDR table panel reports the number of significant gene sets at different FDR thresholds, for all contrasts and all methods. ",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("100%", "100%")
        )
      )
    )
  )

  div(
    boardHeader(title = "Geneset enrichment", info_link = ns("gs_info")),
    bslib::layout_columns(
      col_widths = 12,
      height = fullH,
      gap = "0px",
      tabs1,
      tabs2
    )
  )
}
