##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

DataViewInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        tags$div(
         HTML("<b>DataView.</b> Information and descriptive statistics to quickly lookup a gene, check the total counts, or view the data tables.")
        ),
        shiny::tagList(
            shinyBS::tipify( shiny::actionLink(ns("data_info"), "Tutorial", icon = shiny::icon("youtube")),
                   "Show more information about this module."),
            shiny::hr(), shiny::br(),
            shinyBS::tipify( shiny::selectInput(ns("search_gene"),"Gene:", choices=NULL),
                   "Enter a gene of interest for the analysis.", placement="top"),
            shinyBS::tipify( shiny::selectInput(ns("data_samplefilter"),"Filter samples:",
                                choices=NULL, multiple=TRUE),
                   "Filter the relevant samples for the analysis.", placement="top"),
            shinyBS::tipify( shiny::selectInput(ns('data_groupby'),'Group by:', choices=NULL),
                   "Select phenotype for grouping the samples.", placement="top"),
            shiny::br(),
            shinyBS::tipify( shiny::actionLink(ns("data_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Toggle advanced options.", placement="top"),
            shiny::br(),br(),
            shiny::conditionalPanel(
                "input.data_options % 2 == 1", ns=ns,
                shinyBS::tipify( shiny::radioButtons(ns('data_type'),'Data type:',
                                     choices=c("counts","logCPM"), selected="logCPM", inline=TRUE),
                       "Choose an input data type for the analysis.", placement="bottom")
            )
        )
    )
}

DataViewUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tabsetPanel(
        id = ns("tabs"),
        shiny::tabPanel("Plots",
            shiny::fillCol(
                height = 750,
                flex = c(NA,0.03,1),
                tags$div(
                    HTML("<b>Gene plots.</b> <b>(a)</b> Further information about the selected gene from public databases.
                    <b>(b)</b> Abundance/expression of selected gene across groups. <b>(c)</b>
                    Average rank of the selected gene compared to other genes. <b>(d)</b> t-SNE of samples colored by
                    expression of selected gene. <b>(e)</b> Top correlated genes. Darker color corresponds to higher
                    expression of the gene. <b>(f)</b> Tissue expression of selected gene.")
                ),
                shiny::br(),
                shiny::fillRow(
                    flex = c(1,0.06,5),
                    plotWidget(ns("data_geneInfo")),
                    shiny::br(),
                    shiny::fillCol(
                        flex = c(1,0.2,1),
                        shiny::fillRow(
                            flex = c(1.5,1,1), id = "genePlots_row1",
                            height = 355, ## width=1600,
                            plotWidget(ns("genePlots_barplot")),
                            plotWidget(ns("genePlots_averageRankPlot")),
                            plotWidget(ns("genePlots_tsne"))
                        ),
                        shiny::br(),
                        shiny::fillRow(
                            flex = c(1.5,2), id = "genePlots_row2",
                            height = 355, ## width=1600,
                            plotWidget(ns("genePlots_correlationplot")),
                            plotWidget(ns("data_tissueplot"))
                        )
                    )
                )
            )),
        shiny::tabPanel("QC",
            shiny::fillCol(
                flex = c(NA,0.04,1,1),
                height = 750,
                tags$div(
                    HTML("<b>Counts distribution</b>. Plots associated with the counts, abundance or expression levels across
                    the samples/groups.  <b>(a)</b> Total counts per sample or average per group.
                    <b>(b)</b> Distribution of total counts per sample/group. The center horizontal bar correspond to
                    the median.  <b>(c)</b> Histograms of total counts distribution per sample/group. <b>(d)</b>
                    Abundance of major gene types per sample/group. <b>(e)</b> Average count by gene type per sample/group.")
                ),
                shiny::br(),
                shiny::fillRow(
                    flex = c(1,1,1), id = "counts_tab_row1", height=355,
                    plotWidget(ns("counts_tab_barplot")),
                    plotWidget(ns("counts_tab_boxplot")),
                    plotWidget(ns("counts_tab_histplot"))
                ),
                shiny::fillRow(
                    flex = c(1,1), id = "counts_tab_row2", height=355,
                    plotWidget(ns("counts_tab_abundanceplot")),
                    plotWidget(ns("counts_tab_average_countplot"))
                )
            )),
        shiny::tabPanel("Counts",
            shiny::fillCol(
                flex = c(NA,0.025,1),
                height = 750,
                tags$div(
                    HTML("<b>Gene table.</b> The table shows the gene expression values per sample, or average expression
                    values across the groups. The column 'rho' reports the correlation with the gene selected
                    in 'Search gene' in the left side bar.")
                ),
                shiny::br(),
                tableWidget(ns("data_rawdataTable"))
            )),
        shiny::tabPanel("Samples",
            shiny::fillCol(
                flex = c(NA,0.04,1.2,1),
                height = 750,
                tags$div(
                    HTML(paste(
                        "<b>(a)</b> <b>Phenotype clustering.</b> Clustered heatmap of sample information (i.e. phenotype data).",
                        "<b>(b)</b> <b>Phenotype association matrix.</b> Clustered heatmap of phenotype association. The values corresponds to the -log10(p) value of the corresponding statistical test between two phenotype variables. A higher value corresponds to stronger 'correlation'.",
                        "<b>(c)</b> <b>Sample information table.</b> Phenotype information about the samples. Phenotype variables starting with a 'dot' (e.g. '.cell cycle' and '.gender' ) have been estimated from the data."
                    ))
                ),
                shiny::br(),
                shiny::fillRow(
                    flex = c(2,0.07,1),
                    shiny::div(plotWidget(ns("data_phenoHeatmap")), style="overflow-y: auto;"),
                    shiny::br(),
                    plotWidget(ns("data_phenotypeAssociation"))
                ),
                tableWidget(ns("data_sampleTable"))
            )),
        shiny::tabPanel("Contrasts",
            shiny::fillCol(
                flex = c(NA,0.03,1), height = 750,
                tags$div(
                    HTML("<b>Contrast table.</b> summarizing the contrasts of all comparisons. Non-zero entries '+1' and '-1'
                    correspond to the group of interest and control group, respectively. Zero or empty entries
                    denote samples not use for that comparison.")
                ),
                shiny::br(),
                tableWidget(ns("data_contrastTable"))
            )),
        shiny::tabPanel("Resource info",
            shiny::fillCol(
                flex = c(NA,0.02,1),
                height = 750,
                tags$div(
                    HTML("<b>Resource information.</b> Details about the execution times of the methods,
                     dimensions and memory sizes of objects.")
                ),
                shiny::br(),
                shiny::fillRow(
                    flex = c(5,1, 2,1, 1.5, 2),
                    tableWidget(ns("datatable_timings")),
                    shiny::br(),
                    tableWidget(ns("datatable_objectdims")),
                    shiny::br(),
                    tableWidget(ns("datatable_objectsize")),
                    shiny::br()
                )
            ))
    )
}