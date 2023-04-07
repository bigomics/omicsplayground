##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

SingleCellInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), br(),
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
          "Specify a layout for the figures: t-SNE or PCA-based layout.",
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
    boardHeader(title = "Single Cell Board", info_link = ns("infotext")),
    div(
      shiny::fillCol(
        flex = c(1),
        height = 780,
        shiny::tabsetPanel(
          id = ns("tabs"),
          shiny::tabPanel(
            "Cell type",
            bslib::layout_column_wrap(
              width = 1/2,
              singlecell_plot_icpplot_ui(
                  id = ns("icpplot"),
                  title = "Cell type profiling",
                  info.text = " Currently, we have implemented a total of 8 methods and 9 reference datasets to predict immune cell types (4 datasets), tissue types (2 datasets), cell lines (2 datasets) and cancer types (1 dataset). However, we plan to expand the collection of methods and databases and to infer other cell types.",
                  caption = "Plot infering the type of cells using computational deconvolution methods and reference datasets from the literature.",
                  label = "a",
                  height = c("70vh", modH),
                  width = c("auto", "100%"),
                  parent = ns
              ),
              singlecell_plot_phenoplot_ui(
                  id = ns("phenoplot"),
                  title = "Phenotypes",
                  info.text = "The plots show the distribution of the phenotypes superposed on the t-SNE clustering. Often, we can expect the t-SNE distribution to be driven by the particular phenotype that is controlled by the experimental condition or unwanted batch effects.",
                  caption = "t-SNE plot of the samples with superposed phenotypes.",
                  label = "b",
                  height = c("70vh", modH),
                  width = c("auto", "100%")
              )
            )
          ),
          shiny::tabPanel(
            "Mapping",
            bslib::layout_column_wrap(
              width = 1/2,
              singlecell_plot_mappingplot_ui(
                  id = ns("mappingplot"),
                  title = "Cell type mapping",
                  info.text = "Cell type profiling infers the type of cells using computational deconvolution methods and reference datasets from the literature. Currently, we have implemented a total of 8 methods and 9 reference datasets to predict immune cell types (4 datasets), tissue types (2 datasets), cell lines (2 datasets) and cancer types (1 dataset). However, we plan to expand the collection of methods and databases and to infer other cell types.",
                  caption = "Dot plot indicating the chosen cell type profile of each individual sample/cell or phenotypic group. Useful for determining a sample/cell origin.",
                  label = "a",
                  height = c("70vh", modH),
                  width = c("100%", "100%"),
                  parent = ns
              ),
              singlecell_plot_crosstabPlot_ui(
                  id = ns("crosstabPlot"),
                  title = "Proportions",
                  info.text = "The Proportions tab visualizes the interrelationships between two categorical variables (so-called cross tabulation). Although this feature is very suitable for a single-cell sequencing data, it provides useful information about the proportion of different cell types in samples obtained by the bulk sequencing method.",
                  caption = "Proportion plot indicating the chosen cell-type composition of each experimental phenotype.",
                  label = "b",
                  height = c("70vh", modH),
                  width = c("100%", "100%"),
                  parent = ns
              )
            )
          ),
          shiny::tabPanel(
            "Markers",
            bslib::layout_column_wrap(
              width = 1/2,
              singlecell_plot_markersplot_ui(
                  id = ns("markersplot"),
                  title = "Expression of marker genes",
                  info.text = "The Markers section produces for the top marker genes, a t-SNE with samples colored in red when the gene is overexpressed in corresponding samples. The top genes (N=36) with the highest standard deviation are plotted. In the plotting options, users can also restrict the marker analysis by selecting a particular functional group in which genes are divided into 89 groups, such as chemokines, transcription factors, genes involved in immune checkpoint inhibition, and so on.",
                  caption = "Selected Marker gene expression across all samples.",
                  label = "a",
                  height = c("70vh", modH),
                  width = c("100%", "100%"),
                  parent = ns
              ),
              singlecell_plot_cytoplot_ui(
                  id = ns("cytoplot"),
                  title = "Cytometry plot",
                  info.text = "The aim of this feature is to observe the distribution of samples in relation to the selected gene pairs. For instance, when applied to single-cell sequencing data from immunological cells, it can mimic flow cytometry analysis and distinguish T helper cells from the other T cells by selecting the CD4 and CD8 gene combination.",
                  caption = "Cytometry-like plot of samples based on gene pair combinations. ",
                  label = "b",
                  height = c("70vh", modH),
                  width = c("100%", "100%"),
                  parent = ns
              )
            )
          )
          # These tabPanels are not shown in the UI Stefan 22.03.2022
          # shiny::tabPanel("CNV",
          #     shiny::fillCol(
          #     height = 750,
          #     flex = c(NA,1),
          #     tags$div(
          #         HTML("<b>Inferred copy number variation (CNV)</b>. Genomic copy number is estimated from gene expression
          #         data by computing a moving average of the relative expression along the chromosomes.
          #         The heatmap shows the inferred copy number of the samples versus chromosomes.
          #         The samples are annotated further with phenotype information on the right side of the figure.")
          #     ),
          #     plotWidget(ns("cnaplot"))
          # )),
          # shiny::tabPanel("iTALK",
          #     shiny::fillCol(
          #     height = 750,
          #     flex = c(NA,NA,1),
          #     tags$div(
          #         HTML("<b>Visualization of ligand-receptor interaction.</b> The iTALK R package is designed to profile
          #         and visualize the ligand-receptor mediated intercellular cross-talk signals from single-cell
          #         RNA sequencing data (scRNA-seq). iTALK uses a manually curated list of ligand-receptor gene
          #         pairs further classified into 4 categories based on the primary function of the ligand:
          #         cytokines/chemokines, immune checkpoint genes, growth factors, and others. <b>(a)</b>
          #         The Ligand-Receptor plot visualizes the communication structure of ligand-receptor genes
          #         as a circle plot. <b>(b)</b> The heatmap visualizes the expression level/log fold change
          #         of the ligand/receptor genes.  <b>(c)</b> The NetView plot visualizes the communication structure
          #         of ligand-receptor genes as a graph..")
          #     ),
          #     shiny::inputPanel(
          #         withTooltip( shiny::selectInput(ns("italk_groups"),NULL, choices=""),
          #             "Select the phenotype parameter to divide samples into groups.",
          #             placement="right"),
          #         withTooltip( shiny::selectInput(ns("italk_category"),NULL,
          #                             choices=c('cytokine','growth factor','checkpoint','other')),
          #             "Select the gene category for finding L/R pairs.",
          #             placement="right")
          #     ),
          #     shiny::fillRow(
          #         flex = c(1,1),
          #         plotWidget(ns("italk_LRPlot")),
          #         plotWidget(ns("italk_heatmap")),
          #         plotWidget(ns("italk_netview"))
          #     )
          # )),
          # shiny::tabPanel("Monocle",
          #     shiny::fillCol(
          #     height = fullH,
          #     flex = c(NA,0.05,1),
          #     tags$div(
          #         HTML("<b>Single-cell trajectory analysis</b>.  <b>(a)</b> Heatmap visualizing the expression
          #         of the group-specific top markers. <b>(b)</b> Single-cell trajectory analysis how cells
          #         choose between one of several possible end states. Reconstruction algorithms can robustly
          #         reveal branching trajectories, along with the genes that cells use to navigate these decisions.
          #         <b>(c)</b> Gene expression distribution for selected marker gene.")
          #     ),
          #     shiny::br(),
          #     ##------- Page layout -------
          #     shiny::fillRow(
          #         flex = c(1.2,1),
          #         plotWidget(ns("monocle_plotTopMarkers")),
          #         shiny::fillCol(
          #             flex = c(1,1),
          #             plotWidget(ns("monocle_plotTrajectory")),
          #             plotWidget(ns("monocle_plotGene"))
          #         )
          #     )
          # ))
        )
      )
    )
  )
}
