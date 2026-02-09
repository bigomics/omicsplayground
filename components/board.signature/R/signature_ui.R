##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

style0 <- "font-size: 0.9em; color: #24A; background-color: #dde6f0; border-style: none; padding:0; margin-top: -15px;"

SignatureInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    withTooltip(
      shiny::selectInput(ns("type"),
        label = "Signature type:",
        ##  choices = c("contrast", "<custom>", "hallmark", "KEGG")
        choices = c("<custom>", "contrast")
      ),
      "Specify the type of signature of an interest. Users can choose between custom signature, a contrast profile.",
      placement = "top", options = list(container = "body")
    ),
    shiny::br(),
    shiny::conditionalPanel(
      "input.type == '<custom>'",
      ns = ns,
      shiny::tagList(
        withTooltip(
          shiny::textAreaInput(ns("genelist"), tspan("Genes:"),
            value = NULL,
            rows = 12, placeholder = "Paste your gene list"
          ),
          "Paste a list of genes that defines your signature.",
          placement = "top",
          options = list(container = "body")
        ),
        shiny::br(),
        withTooltip(
          shiny::actionButton(ns("example2"), NULL, style = style0),
          "Use the list of genes involved in apoptosis as a signature."
        ),
        withTooltip(
          shiny::actionButton(ns("example3"), NULL, style = style0),
          "Use the list of genes involved in cell cycle as a signature."
        ),
        withTooltip(
          shiny::actionButton(ns("example1"), NULL, style = style0),
          "Use the list of genes involved in immune checkpoint as a signature."
        ),
        withTooltip(
          shiny::actionButton(
            ns("compute_button"),
            label = "Compute",
            class = "btn-outline-primary",
            icon = icon("refresh")
          ),
          "Click to compute.",
          placement = "right"
        )
      )
    ),
    shiny::conditionalPanel(
      "input.type != '<custom>'",
      ns = ns,
      withTooltip(
        shiny::selectInput(ns("feature"), "Signature:",
          choices = NULL, selected = NULL
        ),
        "Select a specific signature.",
        placement = "top",
        options = list(container = "body")
      )
    )
  )
}

SignatureUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  left.panel <- shiny::tabsetPanel(
    id = ns("tabs1"),
    ## ----------------------------- volcano panel  ------------------
    shiny::tabPanel(
      "Volcano plots",
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 181px)",
        bs_alert("Overlay your custom list of genes on top of the volcano plots for each comparison. You can enter your list of genes on the right."),
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          signature_plot_volcano_ui(
            ns("volcanoPlots"),
            title = "Volcano plots",
            info.text = "Volcano plot of the geneset enrichment. If no contrast is selected on the Enrichment by contrasts table, a volcano for each contrast is shown, if a contrast is selected on the table, then that one contrast alone is displayed. By default, the genes selected on {Genes} settings are highlighted, if a certain gene is to be highlighted, select it from the Genes in signature table. The plot can be colored by using the {Color up/down regulated} plot setting.",
            info.methods = "Statistical testing of differential enrichment of genesets is performed using an aggregation of multiple statistical methods: Fisher’s exact test, fGSEA [1], Camera [2] and GSVA/limma [3] [4]. The maximum q-value of the selected methods is taken as aggregate meta.q value, which corresponds to taking the intersection of significant genes from all tests. As each method uses different estimation parameters (NES for GSEA, odd-ratio for fisher, etc.) for the effect size, for consistency, the average log fold-change of the genes in the geneset as sentinel value is taken. For positive enrichment, genes of the query signature would fall on the upper right of the volcano plot, for negative enrichment, on the upper left.",
            info.references = list(
              list(
                "Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N., & Sergushichev, A. (2016). Fast gene set enrichment analysis. biorxiv, 060012.",
                "https://doi.org/10.1101/060012"
              ),
              list(
                "Wu, D., & Smyth, G. K. (2012). Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic acids research, 40(17), e133-e133.",
                "https://doi.org/10.1093/nar/gks461"
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
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#statistical-testing",
            caption = "Volcano plots visualising the test signature in all available contrasts.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          )
        )
      )
    ),
    ## ----------------------------- enrichment panel  ------------------
    shiny::tabPanel(
      "Enrichment",
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 181px)",
        bs_alert("This panel shows your custom list of genes on top of the GSEA enrichment plots for each comparison. Enter your list of genes in the right box."),
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          signature_plot_enplots_ui(
            ns("enplots"),
            title = "Enrichment plots",
            info.text = "Plot showing the top enriched genes. If nothing is selected on the Enrichment by contrasts table, all contrast plots will be shown, while if a contrast is selected on the table, it will be displayed alone. The genes selected on {Genes} settings are highlighted.",
            info.methods = "Statistical testing of differential enrichment of genesets is performed using an aggregation of multiple statistical methods: Fisher’s exact test, fGSEA [1], Camera [2] and GSVA/limma [3] [4]. The maximum q-value of the selected methods is taken as aggregate meta.q value, which corresponds to taking the intersection of significant genes from all tests. As each method uses different estimation parameters (NES for GSEA, odd-ratio for fisher, etc.) for the effect size, for consistency, the average log fold-change of the genes in the geneset as sentinel value is taken. Positive enrichment means that this particular contrast shows similar expression changes as the query signature.",
            info.references = list(
              list(
                "Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N., & Sergushichev, A. (2016). Fast gene set enrichment analysis. biorxiv, 060012.",
                "https://doi.org/10.1101/060012"
              ),
              list(
                "Wu, D., & Smyth, G. K. (2012). Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic acids research, 40(17), e133-e133.",
                "https://doi.org/10.1093/nar/gks461"
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
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#statistical-testing",
            caption = "Gene set enrichment plots indicating the type of correlation of the test signature with the available contrast profiles.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          )
        )
      )
    ),
    ## ----------------------------- overlap panel ------------------
    shiny::tabPanel(
      "Overlap/similarity",
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 180px)",
        bs_alert("This panel compares other gene sets with your custom list of genes to find similar genesets. Similarity is measured using Fisher's test."),
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c(1, 1),
          signature_plot_overlap_ui(
            ns("overlapScorePlot"),
            title = "Signature overlap scores",
            info.text = "Barplot of the overlap scores for the loaded data against multiple public genesets databases. The plot displays by default the top 60 overlap scores, it can be changed using the {Number of features} plot setting, also the feature names can be toggled with the {Show feature names} plto setting.",
            info.methods = "The overlap score of the gene set combines the odds ratio and significance (q-value) of the Fisher's test.",
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#statistical-testing",
            caption = "The plot shows the gene sets most correlated with the test signature.",
            width = c("auto", "100%"),
            height = c("40%", TABLE_HEIGHT_MODAL)
          ),
          signature_table_overlap_ui(
            ns("overlapTable"),
            title = "Overlap with other signatures",
            info.text = "Under the Overlap/similarity tab, users can find the similarity of their gene list with all the gene sets and pathways in the platform, including statistics such as the total number of genes in the gene set (K), the number of intersecting genes between the list and the gene set (k), the overlapping ratio of k/K, logarithm of the odds ratio (log.OR), as well as the p and q values by the Fisher’s test for the overlap test.",
            caption = "The table indicates the gene sets available in the platform that are most correlated with the tested signature.",
            height = c("60%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    ),
    ## ----------------------------- panel markers ------------------
    shiny::tabPanel(
      "Markers",
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 181px)",
        bs_alert("The markers plot shows the expression levels of the tested genes in the dataset samples as a colored t-SNE plot in red (highly expressed) and light grey (low expressed). The first figure shows the single-sample enrichment of your signature list in red (upregulation) and blue (downregulation)."),
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          signature_plot_markers_ui(
            ns("markers"),
            title = "Markers plot",
            info.text = "Scatter plot displaying the expression levels of the tested genes in the dataset samples as a colored t-SNE plot. The plots can be sorted by correlation, probability or name using the {Sort by} plot setting, the layout can also be modified using the {Layout} plot setting.",
            info.methods = "Using the gene list specified under the settings {Genes}, a t-SNE plot (using the Rtsne R package [1]) is produced of samples for each gene. The samples are colored with respect to the high expression (in red) or low expression (in grey) of that particular gene. The first figure shows the single-sample enrichment of your signature list in red (upregulation) and blue (downregulation).",
            info.references = list(
              list(
                "Krijthe JH (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using Barnes-Hut Implementation. R package version 0.17.",
                "https://doi.org/10.32614/CRAN.package.Rtsne"
              )
            ),
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
            caption = "t-SNE plot showing the expression levels of the tested genes in each of the dataset samples.",
            height = c("100%", TABLE_HEIGHT_MODAL)
          )
        )
      )
    )
  )

  right.panel <- shiny::tabsetPanel(
    id = ns("tabs2"),
    shiny::tabPanel(
      "Enrichment table",
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 181px)",
        signature_table_enrich_by_contrasts_ui(
          ns("enrichmentContrastTable"),
          title = "Enrichment by contrasts",
          info.text = "Enrichment scores of query signature across all contrasts. The table summarizes the enrichment statistics of the gene list in all contrasts using the GSEA algorithm. The NES corresponds to the normalized enrichment score of the GSEA analysis.",
          caption = "Table showing the overall enrichment scores of the tested signature in the available contrasts.",
          height = c("40%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),
        signature_table_genes_in_signature_ui(
          ns("enrichmentGeneTable"),
          title = "Genes in signature",
          info.text = "Genes of the current signature corresponding to the selected contrast. Genes are sorted by decreasing (absolute) fold-change.",
          caption = "Table indicating the expression levels of the genes of the tested signature in the available contrasts.",
          height = c("60%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ),
    shiny::tabPanel(
      "AI Summary",
      bslib::layout_columns(
        col_widths = 12,
        height = "calc(100vh - 181px)",
        AISummaryCardUI(
          ns("aiSummary"),
          title = "AI Summary",
          info.text = "AI-generated summary of the gene signature overlap and enrichment analysis. The summary interprets the top overlapping gene sets and key genes to identify the biological function of the signature.",
          caption = "AI-generated interpretation of the signature analysis results.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    )
  )

  div(
    boardHeader(title = "Test signatures", info_link = ns("info")),
    bslib::layout_columns(
      col_widths = c(8, 4),
      left.panel,
      right.panel
    )
  )
}
