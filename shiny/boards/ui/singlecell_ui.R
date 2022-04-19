##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

SingleCellInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    bigdash::tabSettings(
        withTooltip(shiny::actionLink(ns("info"), "Info", icon=icon("info-circle")),
                "Show more information about this module."),
        shiny::hr(),br(),
        withTooltip( shiny::actionLink(ns("options"), "Options",
                            icon=icon("cog", lib = "glyphicon")),
                "Toggle options", placement="top"),
        shiny::br(),br(),
        shiny::conditionalPanel(
            "input.options % 2 == 1", ns=ns, 
            shiny::tagList(
                withTooltip(shiny::selectInput(ns("samplefilter"),"Filter samples:",
                                    choices=NULL, multiple=TRUE),
                        "Filter relevant samples (cells).",
                        placement="top", options = list(container = "body")),
                
                withTooltip(shiny::selectInput(ns('clustmethod'),"Layout", c("default","pca"),
                                    selected="default"),
                        "Specify a layout for the figures: t-SNE or PCA-based layout.",
                        placement="top", options = list(container = "body"))
            )
        )
    )
}

SingleCellUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        flex = c(1),
        height = 780,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("Cell type",
                shiny::fillRow(
                flex = c(1.5,0.08,1),
                shiny::fillCol(
                    height = 750,
                    flex = c(NA,1),
                    tags$div(
                    HTML("<strong>Cell type profiling</strong> infers the type of cells using computational deconvolution
                        methods and reference datasets from the literature. Currently, we have implemented a total
                        of 8 methods and 9 reference datasets to predict immune cell types (4 datasets),
                        tissue types (2 datasets), cell lines (2 datasets) and cancer types (1 dataset).
                        However, we plan to expand the collection of methods and databases and to infer other cell types.")
                    ),
                    plotWidget(ns("icpplot"))
                ),
                shiny::br(),
                shiny::fillCol(
                    height = 750,
                    flex = c(NA,1), 
                    tags$div(
                    HTML("<b>Phenotype plots.</b> The plots show the distribution of the phenotypes superposed on the
                     t-SNE clustering. Often, we can expect the t-SNE distribution to be driven by the particular
                      phenotype that is controlled by the experimental condition or unwanted batch effects.")
                    ),
                    plotWidget(ns("phenoplot"))
                )
            )),
            shiny::tabPanel("Mapping",
                shiny::fillRow(
                flex = c(1.3,0.08,1),
                shiny::fillCol(
                    height = 750,
                    flex = c(NA,1),
                    tags$div(
                    HTML("<b>Cell type mapping.</b> The inferred cell types can be by matched to the phenotype variable
                    of the data set. The reference set can be a cell type reference database but also cancer types,
                    tissue types or cell lines.")
                    ),
                    plotWidget(ns("mappingplot"))
                ),
                shiny::br(),
                shiny::fillCol(
                    height = 750,
                    flex = c(NA,1),
                    tags$div(
                    HTML("<b>Proportion plot.</b> Plot visualizing the overlap between two categorical variables
                    (so-called cross tabulation). Although this feature is very suitable for a single-cell sequencing data,
                    it provides useful information about the proportion of different cell types in samples
                    obtained by the bulk sequencing method.")
                    ),
                    plotWidget(ns("crosstabPlot"))
                )
            )),
            shiny::tabPanel("Markers",
                shiny::fillRow(
                flex = c(1.3,0.08,1),
                shiny::fillCol(
                    height = 750,
                    flex = c(NA,1),
                    tags$div(
                    HTML("<b>T-SNE distribution of expression of marker genes.</b> Good biomarkers will show a
                    distribution pattern strongly correlated with some phenotype. The top genes with the highest
                    standard deviation are shown. The red color shading is proportional to the (absolute)
                    expression of the gene in corresponding samples.")
                    ),
                    plotWidget(ns("markersplot")) 
                ),
                shiny::br(),
                shiny::fillCol(
                    height = 750,
                    flex = c(NA,1),
                    tags$div(
                    HTML("<b>Cyto plot.</b> This plot shows the distribution of samples in relation to the expression
                    of selected gene pairs. It mimics the scatter plots used for gating in flow cytometry analysis.")
                    ),
                    plotWidget(ns("cytoplot"))
                )
            ))
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
}

