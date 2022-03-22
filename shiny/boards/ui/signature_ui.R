##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

style0 = "font-size: 0.9em; color: #24A; background-color: #dde6f0; border-style: none; padding:0; margin-top: -15px;"

SignatureInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        tags$div(
         HTML("<b>Signature Analysis.</b> Users can test their gene signature by
                calculating an enrichment score. Upload your own gene list, or select
                a contrast which then takes the top differentially expressed genes as
                signature."
            )
        ),
        shiny::tagList(
            shiny::tags$head(shiny::tags$style("#sig-genelistUP.form-control {font-size:11px !important;padding:3px;height:200px;}")),
            shinyBS::tipify( shiny::actionLink(ns("info"), "Tutorial", icon = shiny::icon("youtube")),
                   "Show more information about this module"),
            shiny::hr(), shiny::br(),
            shinyBS::tipify(shiny::textAreaInput(ns("genelistUP"), "Genes:", value = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA",
                                 rows=15, placeholder="Paste your gene list"),
                   "Paste a list of signature genes.", placement="top",
                   options = list(container = "body")),
            shiny::br(),
            shinyBS::tipify(shiny::actionButton(ns("example2"),"[apoptosis] ", style=style0),
                   "Use the list of genes involved in apoptosis as a signature."),
            shinyBS::tipify(shiny::actionButton(ns("example3"),"[cell_cycle] ", style=style0),
                   "Use the list of genes involved in cell cycle as a signature."),
            shinyBS::tipify(shiny::actionButton(ns("example1"),"[immune_chkpt] ", style=style0),
                   "Use the list of genes involved in immune checkpoint as a signature."),
            shiny::br(),br(),
            shinyBS::tipify( shiny::actionLink(ns("options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Toggle advanced options.", placement="top"),
            shiny::br(),
            shiny::conditionalPanel(
                "input.options % 2 == 1", ns=ns,
                shiny::tagList(
                    shinyBS::tipify(shiny::selectInput(ns("type"), label="Signature type:",
                                       choices=c("<custom>","contrast","hallmark","KEGG")),
                           "Specify the type of signature of an interest. Users can choose between custom signature, a contrast profile, or some predefined gene sets including Hallmark and KEGG pathways.",
                           placement="top", options = list(container = "body")),
                    shiny::conditionalPanel(
                        "input.type != '<custom>'", ns=ns,
                        shinyBS::tipify(shiny::selectInput(ns("feature"),"Signature:",
                                           choices="<custom>", selected="<custom>"),
                               "Select a specific signature group.", placement="top",
                               options = list(container = "body"))
                    )
                )
            )
        )
    )
}

SignatureUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillRow(
        flex = c(1.5,0.05,1),
        height = 780,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Enrichment",uiOutput(ns("enplots_UI"))),
            shiny::tabPanel("Volcano plots",uiOutput(ns("volcanoPlots_UI"))),
            shiny::tabPanel("Overlap/similarity",uiOutput(ns("overlapAnalysis_UI"))),
            shiny::tabPanel("Markers",uiOutput(ns("markers_UI")))
        ),
        shiny::br(),
        shiny::tabsetPanel(
            id = ns("tabs2"),
            shiny::tabPanel("Enrichment table",uiOutput(ns("enrichmentTables_UI")))
        )
    )
}