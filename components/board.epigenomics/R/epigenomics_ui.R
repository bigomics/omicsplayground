## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.

EpigenomicsInputs <- function(id) {

  ns <- shiny::NS(id)

  bigdash::tabSettings(

    withTooltip(
      shiny::selectizeInput(
        inputId = ns("select_chromosome"),
        label = "Chromosome:",
        choices = NULL,
        multiple = TRUE
      ),
      "Select a chromosome of interest.",
      placement = "top"
    ),

    withTooltip(
      shiny::selectizeInput(
        inputId = ns("search_gene"),
        label = tspan("Gene:"),
        choices = NULL,
        options = list(maxOptions = 1001)
      ),
      "Type a gene of interest.",
      placement = "top"
    ),

    withTooltip(
      shiny::selectInput(
        inputId = ns("data_samplefilter"),
        label = "Filter out samples:",
        choices = NULL,
        multiple = TRUE
      ),
      "Remove specified samples",
      placement = "top"
    ),

    withTooltip(
      shiny::selectInput(
        inputId = ns("select_pheno"),
        label = "Select phenotype:",
        choices = NULL
      ),
      "Select phenotype to visualize.",
      placement = "top"
    )
  )
}

EpigenomicsUI <- function(id) {
  ns <- shiny::NS(id)

  div(
    boardHeader(title = "Epigenomics", info_link = ns("board_info")),
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "Beta Ideograms",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          epigenomics_plot_methylIdeogram_ui(
            id = ns("methylIdeogram"),
            title = "Methylation Ideogram",
            label = "a",
            caption = "Chromosome ideogram showing average beta values per genomic region.",
            info.text = "Displays methylation beta values across chromosomes' ideogram.",
            info.methods = "Beta values are used. If needed, beta values are computed from M values. CpG probes are mapped to chromosomal coordinates using 'hg19' reference genome. The gray line represents the average beta value per region, per chromosome. Specifically, the per-CpG mean is first calculated across all available samples. Note that samples can optionally be filtered out using the drop-down menu. All CpGs mapped within each 1Mb window are identified and the arithmetic mean is calculated. Therefore, effectively, the average of the per-sample-averaged beta values is calculated within each 1Mb chromosomal bin. To aid visualization of dense 450K and EPIC arrays, each 1Mb bin mean is then smoothed using a probe-count weighted LOESS fit with a fixed local neighborhood of 5 bins (~5 Mb), so that denser regions anchor the curve more strongly than sparse ones. The span is set dynamically (loess_bins / total_bins_on_chromosome), with default loess_bins = 5). Thus, local neighbourhood always covers ~5 adjacent bins (~5 Mb) regardless of chromosome length. Bins with more probes carry more weight in the local fit, downweighting sparse  regions. The LOESS-predicted values at each bin midpoint are clamped to [0,1] and are connected as a continuous line. Points in blue are CpG probes with beta values <= 0.2 (hypomethylated). Points in yellow are CpG probes with beta values within the 0.2 - 0.8 range. Points in red are CpG probes with beta values >= 0.8 (hypermethylated). The panel placed below each chromosome ideogram shows a barplot of the number of CpG probes in the 1Mb region.",
            info.references = NULL,
            info.extra_link = NULL,
            height = c("100%", "600px"),
            width = c("100%", "100%")
          )
        )
      )
    )
  )
}
