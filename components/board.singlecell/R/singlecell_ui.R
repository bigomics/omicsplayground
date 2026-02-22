##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

SingleCellInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    withTooltip(
      shiny::selectInput(ns("samplefilter"), "Filter samples:",
        choices = NULL, multiple = TRUE
      ),
      "Filter relevant samples (cells).",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(
      shiny::selectInput(ns("clustmethod"), "Layout", c("default", "pca"),
        selected = "default"
      ),
      "Specify a layout for the figures: t-SNE, UMAP or PCA-based layout.",
      placement = "top", options = list(container = "body")
    )
  )
}

SingleCellUI <- function(id) {
  fullH <- 750 ## full height of panel
  imgH <- 680 ## row height of panel
  tabH <- 200 ## row height of panel
  modH <- TABLE_HEIGHT_MODAL

  ns <- shiny::NS(id) ## namespace
  div(
    boardHeader(title = "Cell Profiling", info_link = ns("infotext")),
    shiny::tabsetPanel(
      id = ns("tabs"),
      shiny::tabPanel(
        "Cell type",
        bslib::layout_columns(
          col_widths = c(6, 6),
          height = "calc(100vh - 181px)",
          singlecell_plot_icpplot_ui(
            id = ns("icpplot"),
            title = "Cell type profiling",
            info.text = " Currently, we have implemented a total of 8 methods and 9 reference datasets to predict immune cell types (4 datasets), tissue types (2 datasets), cell lines (2 datasets) and cancer types (1 dataset). However, we plan to expand the collection of methods and databases and to infer other cell types.",
            caption = "Plot infering the type of cells using computational deconvolution methods and reference datasets from the literature.",
            label = "a",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%"),
            parent = ns
          ),
          singlecell_plot_phenoplot_ui(
            id = ns("phenoplot"),
            title = "Phenotypes",
            info.text = "The plots show the distribution of the phenotypes superposed on the t-SNE clustering. Often, we can expect the t-SNE distribution to be driven by the particular phenotype that is controlled by the experimental condition or unwanted batch effects.",
            caption = "t-SNE plot of the samples with superposed phenotypes.",
            label = "b",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      ),
      shiny::tabPanel(
        "Mapping",
        bslib::layout_columns(
          col_widths = c(6, 6),
          height = "calc(100vh - 181px)",
          singlecell_plot_mappingplot_ui(
            id = ns("mappingplot"),
            title = "Cell type mapping",
            info.text = "Cell type profiling infers the type of cells using computational deconvolution methods and reference datasets from the literature. Currently, we have implemented a total of 8 methods and 9 reference datasets to predict immune cell types (4 datasets), tissue types (2 datasets), cell lines (2 datasets) and cancer types (1 dataset). However, we plan to expand the collection of methods and databases and to infer other cell types.",
            caption = "Dot plot indicating the chosen cell type profile of each individual sample/cell or phenotypic group. Useful for determining a sample/cell origin.",
            label = "a",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%"),
            parent = ns
          ),
          singlecell_plot_crosstabPlot_ui(
            id = ns("crosstabPlot"),
            title = "Proportions",
            info.text = "The Proportions tab visualizes the interrelationships between two categorical variables (so-called cross tabulation). Although this feature is very suitable for a single-cell sequencing data, it provides useful information about the proportion of different cell types in samples obtained by the bulk sequencing method.",
            caption = "Proportion plot indicating the chosen cell-type composition of each experimental phenotype.",
            label = "b",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%"),
            parent = ns
          )
        )
      ),
      shiny::tabPanel(
        "Markers",
        bslib::layout_columns(
          col_widths = c(6, 6),
          height = "calc(100vh - 181px)",
          singlecell_plot_markersplot_ui(
            id = ns("markersplot"),
            title = "Expression of marker genes",
            info.text = "The Markers section produces for the top marker genes, a t-SNE with samples colored in red when the gene is overexpressed in corresponding samples. The top genes (N=36) with the highest standard deviation are plotted. In the plotting options, users can also restrict the marker analysis by selecting a particular functional group in which genes are divided into 89 groups, such as chemokines, transcription factors, genes involved in immune checkpoint inhibition, and so on.",
            caption = "Selected Marker gene expression across all samples.",
            label = "a",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%"),
            parent = ns
          ),
          singlecell_plot_cytoplot_ui(
            id = ns("cytoplot"),
            title = "Cytometry plot",
            info.text = "The aim of this feature is to observe the distribution of samples in relation to the selected gene pairs. For instance, when applied to single-cell sequencing data from immunological cells, it can mimic flow cytometry analysis and distinguish T helper cells from the other T cells by selecting the CD4 and CD8 gene combination.",
            caption = "Cytometry-like plot of samples based on gene pair combinations. ",
            label = "b",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%"),
            parent = ns
          )
        )
      ),
      shiny::tabPanel(
        "AI Summary",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          AiTextCardUI(
            ns("singlecellAISummary"),
            title = "AI Cell Profiling Summary",
            info.text = "AI-generated summary of the cell type deconvolution and marker gene analysis results.",
            caption = "AI-generated cell profiling summary.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    )
  )
}
