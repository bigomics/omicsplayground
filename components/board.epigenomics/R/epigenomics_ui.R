## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.

EpigenomicsInputs <- function(id) {

  ns <- shiny::NS(id)

  bigdash::tabSettings(

    shiny::conditionalPanel(
      condition = sprintf("input['%s'] === 'Methylation ideograms'", ns("tabs")),
      withTooltip(
        shiny::selectizeInput(
          inputId = ns("select_chromosome"),
          label = "Chromosome:",
          choices = NULL,
          multiple = TRUE
        ),
        "Select a chromosome of interest.",
        placement = "top"
      )
    ),

    ## withTooltip(
    ##   shiny::selectizeInput(
    ##     inputId = ns("search_gene"),
    ##     label = tspan("Gene:"),
    ##     choices = NULL,
    ##     options = list(maxOptions = 1001)
    ##   ),
    ##   "Type a gene of interest.",
    ##   placement = "top"
    ## ),

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
    tabs <- shiny::tabsetPanel(

      id = ns("tabs"),

      shiny::tabPanel(
        "Methylomics landscape",
        bs_alert("The Methylomics landscape panel provides an overview of the methylomics profiles across chromosomes and samples."),
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          bslib::layout_columns(
            col_widths = 12,
            row_heights = list(1.33, 1),            
            bslib::layout_columns(
              col_widths = c(6, 6),
              epigenomics_plot_beta_dist_ui(
                id = ns("betaDist"),
                title = "Global beta distribution",
                label = "",
                caption = "Density plot of beta values. ",
                info.text = "Density plot of beta values. Beta values from methylomics array cover the [0-1] range and should exhibit a bimodal distribution with enrichment at around 0.2 (hypomethylation) and 0.8 (hypermethylation).",
                info.methods = "Density plot is drawn from the whole dataset, including all chromosomes and all samples usng the ggplot2::geom_density R function. Dashed lines at beta values 0.2 and 0.8 are depicted.",
                info.references = NULL,
                info.extra_link = NULL,
                height = c("100%", "600px"),
                width = c("100%", "100%")
              ),
              dataview_table_beta_ui(
                id = ns("methyltable"),
                height = c("50%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%"),
                title = "Methylation table",
                info.text = "Table reporting providing information on beta methylation profiles stratified per chromosome. CpG probes are mapped to corresponding chromosome using annotation file. By default, if phenotype is 'ungrouped', the average methylation (beta) value for each chromosome is calculated across all samples ('Ave'), and for each available sample. If a phenotype is selected, the average methylation value is calculated for each group.",
                caption = "Table providing information on beta methylation profiles stratified per chromosome."
              )
            ),
            bslib::layout_columns(
              col_widths = 12,
              epigenomics_plot_boxplot_beta_ui(
                id = ns("boxplotBeta"),
                title = "Beta chromosomal profiles",
                label = "",
                caption = "Boxplots of beta values per chromosome.",
                info.text = "Boxplots of beta values per chromosome. Beta values from methylomics array cover the [0-1] range and should enrich at around 0.2 (hypomethylation) and 0.8 (hypermethylation).",
                info.methods = "Boxplots of beta values per chromosome. For each chromosome, the average beta value if first computed within each sample. For each chromosome, the boxplot show the distribution of average beta values per sample.",
                info.references = NULL,
                info.extra_link = NULL,
                height = c("100%", "600px"),
                width = c("100%", "100%")
              )
            )
          )
        )
      ),

      shiny::tabPanel(
        "Methylation ideograms",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          epigenomics_plot_methylIdeogram_ui(
            id = ns("methylIdeogram"),
            title = "Methylation Ideograms",
            label = "",
            caption = "Chromosome ideogram showing average beta values per genomic region.",
            info.text = "Displays methylation beta values across chromosomes' ideogram.",
            info.methods = "Genome-wide karyotype-style methylation visualization is generated using karyoploteR. It begins by harmonizing inputs: restricting beta matrix and annotation to their shared CpGs, auto-detecting chromosome and position columns by name pattern and removing features with missing coordinates. It then filters to the requested chromosomes (defaulting to chr1–22, X, Y) and computes per-CpG mean beta values across all samples; if a two-level phenotype vector is provided, it also computes separate per-group means. Each CpG is converted into a GRanges object and assigned a scatter colour based on its mean methylation: blue for hypomethylated (β < 0.2), yellow for intermediate, and red for hypermethylated (β > 0.8). Positions falling outside the reference genome's chromosome lengths (hg19 or mm10) are discarded. All CpGs are then binned along each chromosome in fixed-size windows (default 1 Mb), computing the mean beta and CpG count per bin, and applies an adaptive LOESS smooth weighted by bin density to produce a clean trend line; this is done once for the overall mean and, in two-group mode, once per group. If two groups are present, it further computes a per-bin delta-beta (group1 minus group2) by matching bin positions across groups. In the figure, panel 1 shows the scatter of individual CpG means coloured by methylation level (subsampled to 200,000 points for performance), overlaid with dashed reference lines at β = 0.2 and β = 0.8, and a smoothed mean trend line; in two-group mode with pheno_lines = TRUE this becomes two coloured lines (one per group) with a legend placed to the right of the shortest chromosome. Panel 2, if present, is either a delta-beta display with red bars rising above zero where group1 exceeds group2, and blue bars falling below where group2 exceeds group1, with an overlaid delta line and a symmetric axis scaled to the observed range. In single-group mode, log-scaled grey bars showing the probe density per bin as a coverage track. Finally, a black border rectangle is drawn around each chromosome's data area in both panels.",
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
