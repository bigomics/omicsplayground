##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

style0 = "font-size: 0.9em; color: #24A; background-color: #dde6f0; border-style: none; padding:0; margin-top: -15px;"

SignatureInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    bigdash::tabSettings(
        shiny::tags$head(shiny::tags$style("#sig-genelistUP.form-control {font-size:11px !important;padding:3px;height:200px;}")),
        withTooltip( shiny::actionLink(ns("info"), "Tutorial", icon = shiny::icon("youtube")),
                "Show more information about this module"),
        shiny::hr(), shiny::br(),
        withTooltip(shiny::textAreaInput(ns("genelistUP"), "Genes:", value = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA",
                                rows=15, placeholder="Paste your gene list"),
                "Paste a list of signature genes.", placement="top",
                options = list(container = "body")),
        shiny::br(),
        withTooltip(shiny::actionButton(ns("example2"),"[apoptosis] ", style=style0),
                "Use the list of genes involved in apoptosis as a signature."),
        withTooltip(shiny::actionButton(ns("example3"),"[cell_cycle] ", style=style0),
                "Use the list of genes involved in cell cycle as a signature."),
        withTooltip(shiny::actionButton(ns("example1"),"[immune_chkpt] ", style=style0),
                "Use the list of genes involved in immune checkpoint as a signature."),
        shiny::br(),br(),
        withTooltip( shiny::actionLink(ns("options"), "Options", icon=icon("cog", lib = "glyphicon")),
                "Toggle advanced options.", placement="top"),
        shiny::br(),
        shiny::conditionalPanel(
            "input.options % 2 == 1", ns=ns,
            shiny::tagList(
                withTooltip(shiny::selectInput(ns("type"), label="Signature type:",
                                    choices=c("<custom>","contrast","hallmark","KEGG")),
                        "Specify the type of signature of an interest. Users can choose between custom signature, a contrast profile, or some predefined gene sets including Hallmark and KEGG pathways.",
                        placement="top", options = list(container = "body")),
                shiny::conditionalPanel(
                    "input.type != '<custom>'", ns=ns,
                    withTooltip(shiny::selectInput(ns("feature"),"Signature:",
                                        choices="<custom>", selected="<custom>"),
                            "Select a specific signature group.", placement="top",
                            options = list(container = "body"))
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
            shiny::tabPanel("Enrichment",
                shiny::fillCol(
                height = 800,
                flex = c(NA,0.03,1),
                tags$div(
                    HTML("<b>Enrichment plots.</b> Enrichment of the query signature in all constrasts.
                     Positive enrichment means that this particular contrast
                    shows similar expression changes as the query signature."
                    )
                ),
                shiny::br(),
                plotWidget(ns("enplots"))
            )),
            shiny::tabPanel("Volcano plots",
                shiny::fillCol(
                height = 800,
                flex = c(NA,0.03,1),
                tags$div(
                    HTML("
                        <b>Volcano plots.</b> Visualization of the query signature on the volcano plots of all constrasts.
                        For positive enrichment, genes of the query signature would fall on the upper right of the volcano plot,
                        for negative enrichment, on the upper left."
                    )
                ),
                shiny::br(),
                plotWidget(ns("volcanoPlots"))
            )
            ),
            shiny::tabPanel("Overlap/similarity",
                shiny::fillCol(
                flex = c(NA,0.035,1,0.04,1),
                height = 800,
                tags$div(
                    HTML("
                        <b>Overlap/Similarity table.</b><b>(a)</b> Top overlapping gene sets with selected signature.
                        The vertical axis shows the overlap score of the gene set which combines the odds ratio and
                        significance (q-value) of the Fisher's test. <b>(b)</b> Table summarizing the results
                        of the Fishers's test for overlap. The column \'common genes\' reports the
                        shared gene in order of largest fold-change."
                    )
                ),
                shiny::br(),
                plotWidget(ns("overlapScorePlot")),
                shiny::br(),
                tableWidget(ns("overlapTable"))
            )),
            shiny::tabPanel("Markers",
                shiny::fillCol(
                flex = c(NA,0.025,1),
                height = 800,
                tags$div(
                    HTML("
                        <b>Markers t-SNE plot</b>. T-SNE plot for each gene, where the dot (corresponding to samples)
                        are colored depending on the upregulation (in red) or downregulation (in blue) of that particular gene."
                    )
                ),
                shiny::br(),            
                plotWidget(ns("markers"))
            ))
        ),
        shiny::br(),
        shiny::tabsetPanel(
            id = ns("tabs2"),
            shiny::tabPanel("Enrichment table",
                shiny::fillCol(
                flex = c(NA,0.04,1.0,0.04,1.5), ## width = 600,
                height = 800,
                tags$div(
                    HTML("
                        <b>Enrichment of query signature across all contrasts.</b> <b>(a)</b> Enrichment scores across all
                        contrasts for the selected query signature . The NES corresponds to the normalized enrichment
                        score of the GSEA analysis. <b>(b)</b> Genes in the query signature sorted by decreasing (absolute)
                        fold-change corresponding to the selected contrast."
                    )
                ),
                shiny::br(),
                plotWidget(ns("enrichmentContrastTable")),
                shiny::br(),
                plotWidget(ns("enrichmentGeneTable"))
            ))
        )
    )
}