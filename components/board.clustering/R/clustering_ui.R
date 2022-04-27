##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

ClusteringInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    bigdash::tabSettings(
        withTooltip( shiny::actionLink(ns("clust_info"), "Tutorial", icon = shiny::icon("youtube")),
                "Show more information and video tutorial about this module."),
        shiny::hr(), shiny::br(),             
        withTooltip( shiny::selectInput(ns("hm_features"),"Features:", choices=NULL, multiple=FALSE),
                "Select a family of features.", placement="top"),
        shiny::conditionalPanel(
            "input.hm_features == '<custom>'", ns=ns,
            withTooltip( shiny::textAreaInput(ns("hm_customfeatures"), NULL, value = NULL,
                                    height = "150px", width = "100%", 
                                    rows=5, placeholder="Paste your custom gene list"),
                    "Paste a custom list of genes to be used as features.",
                    placement="bottom")
        ),
        shiny::conditionalPanel(
            "input.hm_features == '<contrast>'", ns=ns,
            tipifyR( shiny::selectInput(ns("hm_contrast"), NULL, choices=NULL),
                    "Select contrast to be used as signature.")
        ),
        withTooltip( shiny::selectInput(ns('hm_group'),'Group by:',choices=NULL),
                "Group the samples by condition.", 
                placement="top", options = list(container = "body")),
        withTooltip( shiny::selectInput(ns("hm_samplefilter"),"Filter samples:",
                            choices=NULL, multiple=TRUE),
                "Filter the relevant samples for the analysis.",
                placement="top", options = list(container = "body")),            
        withTooltip( shiny::actionLink(ns("hm_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                "Toggle advanced options.", placement="top"),
        shiny::br(),
        shiny::conditionalPanel(
            "input.hm_options % 2 == 1", ns=ns,
            shiny::tagList(
                        withTooltip( shiny::selectInput(ns("hm_level"),"Level:", choices=c("gene","geneset")),
                                        "Specify the level analysis: gene or geneset level.",
                                        placement="top", options = list(container = "body")),
                        withTooltip( shiny::checkboxInput(ns('hm_filterXY'),'exclude X/Y genes',FALSE),
                                        "Exclude genes on X/Y chromosomes.", 
                                        placement="top", options = list(container = "body")),
                        withTooltip( shiny::checkboxInput(ns('hm_filterMitoRibo'),
                                                                'exclude mito/ribo genes',FALSE),
                                        "Exclude mitochondrial (MT) and ribosomal protein (RPS/RPL) genes.", 
                                        placement="top", options = list(container = "body"))
                )
        )
    )
}

ClusteringUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillRow(
        flex = c(1.6,0.05,1),
        height = 780,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Heatmap", 
                plotWidget(ns("hm_splitmap")),
                tags$div( class="caption",
                    HTML("<b>Clustered heatmap.</b> Heatmap showing gene expression sorted by 2-way hierarchical
                    clustering. Red corresponds to overexpression, blue to underexpression of the gene.
                    At the same time, gene clusters are functionally annotated in the
                    'Annotate clusters' panel on the right."
                    )
                )
            ),
            shiny::tabPanel("PCA/tSNE",
                plotWidget(ns("hm_PCAplot")),
                tags$div( class="caption",
                    HTML("<b>PCA/tSNE plot.</b> The plot visualizes the similarity in expression of
                     samples as a scatterplot in reduced dimension (2D or 3D).
                     Samples that are similar are clustered near to each other, while samples with different
                     expression are positioned farther away. Groups of samples with similar profiles
                     will appear as <i>clusters</i> in the plot.")
                )
            ),
            shiny::tabPanel("Parallel",
                plotWidget(ns("hm_parcoord")),
                tableWidget(ns("hm_parcoord_table")),                            
                tags$div( class="caption",
                    HTML("<b>Parallel Coordinates plot.</b> <b>(a)</b>The Parallel Coordinates plot displays
                        the expression levels of selected genes across all conditions.
                        On the x-axis the experimental conditions are plotted. The y-axis shows the expression level
                        of the genes grouped by condition. The colors correspond to the gene groups as
                        defined by the hierarchical clustered heatmap. <b>(b)</b>
                        Average expression of selected genes across conditions."
                    )
                )
            )
        ),
        shiny::br(),
        shiny::tabsetPanel(
            id = ns("tabs2"),
            shiny::tabPanel("Annotate clusters",
                uiOutput(ns("hm_annotateUI"))),
            shiny::tabPanel("Phenotypes",
                plotWidget(ns("clust_phenoplot")),
                tags$div( class="caption",
                        HTML("<b>Phenotype distribution.</b> The plots show the distribution of the phenotypes
                        superposed on the t-SNE clustering. Often, we can expect the t-SNE distribution to be
                        driven by the particular phenotype that is controlled by the experimental condition
                        or unwanted batch effects."
                        )
                )
            ),
            shiny::tabPanel("Feature ranking",
                plotWidget(ns("clust_featureRank")),      
                tags$div( class="caption",
                        HTML("<b>Feature-set ranking.</b> Ranked discriminant score for top feature sets.
                         The plot ranks the discriminative power of feature sets (or gene sets) as the
                         cumulative discriminant score for all phenotype variables."
                        )
                )
            )
        )
    )
}
