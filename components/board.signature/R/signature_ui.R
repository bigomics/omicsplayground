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

  fullH <- 700 ## full height of page
  fullH <- "70vh" ## full height of page
  halfH <- "35vh" ## full height of page

  tabs <- div(
    class = "row",
    div(
      class = "col-md-8",
      shiny::tabsetPanel(
        id = ns("tabs1"),
        ## ----------------------------- volcano panel  ------------------                
        shiny::tabPanel(
          "Volcano plots",
          signature_plot_volcano_ui(
            ns("volcanoPlots"),
            height = c(fullH, TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          )
        ),
        ## ----------------------------- enrichment panel  ------------------        
        shiny::tabPanel(
          "Enrichment",
          signature_plot_enplots_ui(
            ns("enplots"),
            height = c(fullH, TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          )
        ),
        ## ----------------------------- overlap panel ------------------        
        shiny::tabPanel(
          "Overlap/similarity",
          bslib::layout_column_wrap(
            width = 1,
            signature_plot_overlap_ui(
              ns("overlapScorePlot"),
              width = c("auto", "100%"),
              height = c(halfH, TABLE_HEIGHT_MODAL)
            ),
            signature_table_overlap_ui(
              ns("overlapTable"),
              height = c(halfH, TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        ),
        ## ----------------------------- panel markers ------------------
        shiny::tabPanel(
          "Markers",
          signature_plot_markers_ui(
            ns("markers"),
            height = c(fullH, TABLE_HEIGHT_MODAL)
          )
        )
      )
    ),
    div(
      class = "col-md-4",
      shiny::tabsetPanel(
        id = ns("tabs2"),
        shiny::tabPanel(
          "Enrichment table",
          signature_table_enrich_by_contrasts_ui(
            ns("enrichmentContrastTable"),
            height = c(230, TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          ),
          signature_table_genes_in_signature_ui(
            ns("enrichmentGeneTable"),
            height = c(360, TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    )  ## col
  ) ## tabs as row
  
  div(
    boardHeader(title = "Test signatures", info_link = ns("info")),
    tabs
  )
}
