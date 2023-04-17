##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

style0 <- "font-size: 0.9em; color: #24A; background-color: #dde6f0; border-style: none; padding:0; margin-top: -15px;"

SignatureInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::tags$head(shiny::tags$style("#sig-genelistUP.form-control {font-size:11px !important;padding:3px;height:200px;}")),
    shiny::hr(), shiny::br(),
    withTooltip(
      shiny::textAreaInput(ns("genelistUP"), "Genes:",
        value = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA",
        rows = 15, placeholder = "Paste your gene list"
      ),
      "Paste a list of genes that defines your signature.",
      placement = "top",
      options = list(container = "body")
    ),
    shiny::br(),
    withTooltip(
      shiny::actionButton(ns("example2"), "[apoptosis] ", style = style0),
      "Use the list of genes involved in apoptosis as a signature."
    ),
    withTooltip(
      shiny::actionButton(ns("example3"), "[cell_cycle] ", style = style0),
      "Use the list of genes involved in cell cycle as a signature."
    ),
    withTooltip(
      shiny::actionButton(ns("example1"), "[immune_chkpt] ", style = style0),
      "Use the list of genes involved in immune checkpoint as a signature."
    ),
    shiny::br(), br(),
    withTooltip(shiny::actionLink(ns("options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(),
    shiny::conditionalPanel(
      "input.options % 2 == 1",
      ns = ns,
      shiny::tagList(
        withTooltip(
          shiny::selectInput(ns("type"),
            label = "Signature type:",
            choices = c("<custom>", "contrast", "hallmark", "KEGG")
          ),
          "Specify the type of signature of an interest. Users can choose between custom signature, a contrast profile, or some predefined gene sets including Hallmark and KEGG pathways.",
          placement = "top", options = list(container = "body")
        ),
        shiny::conditionalPanel(
          "input.type != '<custom>'",
          ns = ns,
          withTooltip(
            shiny::selectInput(ns("feature"), "Signature:",
              choices = "<custom>", selected = "<custom>"
            ),
            "Select a specific signature group.",
            placement = "top",
            options = list(container = "body")
          )
        )
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
      bslib::layout_column_wrap(
        width = 1,
        height = "calc(100vh - 190px)",
        ##height = "100%",        
        signature_plot_volcano_ui(
          ns("volcanoPlots"),
          title = "Volcano plots",
          info.text = "For positive enrichment, genes of the query signature would fall on the upper right of the volcano plot, for negative enrichment, on the upper left.",
          caption = "Volcano plots visualising the test signature in all available contrasts.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("100%", "100%")
        )
      )
    ),
    ## ----------------------------- enrichment panel  ------------------        
    shiny::tabPanel(
      "Enrichment",
      bslib::layout_column_wrap(
        width = 1,
        height = "calc(100vh - 190px)",
        signature_plot_enplots_ui(
          ns("enplots"),
          title = "Enrichment plots",
          info.text = "Enrichment of the query signature in all constrasts. Positive enrichment means that this particular contrast shows similar expression changes as the query signature.",
          caption = "Gene set enrichment plots indicating the type of correlation of the test signature with the available contrast profiles.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("100%", "100%")
        )
      )
    ),
    ## ----------------------------- overlap panel ------------------        
    shiny::tabPanel(
      "Overlap/similarity",
      bslib::layout_column_wrap(
        width = 1,
        height = "calc(100vh - 190px)",
        signature_plot_overlap_ui(
          ns("overlapScorePlot"),
          title = "Signature overlap scores",
          info.text = "The vertical axis shows the overlap score of the gene set which combines the odds ratio and significance (q-value) of the Fisher's test.",
          caption = "The plot shows the gene sets most correlated with the test signature.",
          width = c("auto", "100%"),
          height = c("50%", TABLE_HEIGHT_MODAL)
        ),
        signature_table_overlap_ui(
          ns("overlapTable"),
          title = "Overlap with other signatures",
          info.text = "Under the Overlap/similarity tab, users can find the similarity of their gene list with all the gene sets and pathways in the platform, including statistics such as the total number of genes in the gene set (K), the number of intersecting genes between the list and the gene set (k), the overlapping ratio of k/K, logarithm of the odds ratio (log.OR), as well as the p and q values by the Fisher’s test for the overlap test.",
          caption = "The table indicates the gene sets available in the platform that are most correlated with the tested signature.",
          height = c("50%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ),
    ## ----------------------------- panel markers ------------------
    shiny::tabPanel(
      "Markers",
      bslib::layout_column_wrap(
        width = 1,
        height = "calc(100vh - 190px)",
        signature_plot_markers_ui(
          ns("markers"),
          title = "Markers plot",
          info.text = "After uploading a gene list, the Markers section produces a t-SNE plot of samples for each gene, where the samples are colored with respect to the upregulation (in red) or downregulation (in blue) of that particular gene.",
          caption = "t-SNE plot showing the expression levels of the tested genes in each of the dataset samples.",
          height = c("100%", TABLE_HEIGHT_MODAL)
        )
      )
    )
  )

  right.panel <- shiny::tabsetPanel(
    id = ns("tabs2"),
    shiny::tabPanel(
      "Enrichment table",
      bslib::layout_column_wrap(
        width = 1,
        height = "calc(100vh - 190px)",
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
    )
  )
  
  div(
    boardHeader(title = "Test signatures", info_link = ns("info")),
    bslib::layout_column_wrap(
      width = 1,
      ## height = "calc(100vh - 190px)",
      style = htmltools::css(grid_template_columns = "8fr 4fr"),
      left.panel,
      right.panel
    )
  )
}
