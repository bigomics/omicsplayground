##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


#' DataView module UI Input function
#'
#' @description A shiny Module. Renders the input parts (sidebar contents) for the module.
#'
#' @param id Internal parameters for {shiny}.
#' #'
#' @export 
DataViewInputs <- function(id) {
  ns <- shiny::NS(id)  ## namespace

  bigdash::tabSettings(
    withTooltip( shiny::actionLink(ns("data_info"), "Tutorial", icon = shiny::icon("youtube")),
                "Show more information about this module."),
    shiny::hr(), shiny::br(),
    withTooltip( shiny::selectInput(ns("search_gene"),"Gene:", choices=NULL),
                "Enter a gene of interest for the analysis.", placement="top"),
    withTooltip( shiny::selectInput(ns("data_samplefilter"),"Filter samples:",
                                    choices=NULL, multiple=TRUE),
                "Filter the relevant samples for the analysis.", placement="top"),
    withTooltip( shiny::selectInput(ns('data_groupby'),'Group by:', choices=NULL),
                "Select phenotype for grouping the samples.", placement="top"),
    shiny::br(),
    withTooltip( shiny::actionLink(ns("data_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                "Toggle advanced options.", placement="top"),
    shiny::br(), shiny::br(),
    shiny::conditionalPanel(
      "input.data_options % 2 == 1", ns=ns,
      withTooltip( shiny::radioButtons(ns('data_type'),'Data type:',
                                       choices=c("counts","logCPM"), selected="logCPM", inline=TRUE),
                  "Choose an input data type for the analysis.", placement="bottom")
    )
  )
}

#' DataView module UI output function
#'
#' @description Renders the output part for the module as tabsetPanel object
#'
#' @param id Internal parameters for {shiny}.
#' #'
#' @export 
DataViewUI <- function(id) {
  ns <- shiny::NS(id)  ## namespace

  imgH <- c(330,600)   ## heights for small and fullscreen image

  shiny::tabsetPanel(

    id = ns("tabs"),

    ##----------------------------------------------------------------------------
    shiny::tabPanel(
      "Plots",
      shiny::fillRow(
        flex = c(1.2,5),
        dataview_module_geneinfo_ui(ns("geneinfo")),
        shiny::fillCol(
          flex = c(1,1,NA),
          height = 800, 
          shiny::fillRow(
            flex = c(1.5,1,1), id = "genePlots_row1",
            height = 375, ## width=1600,
            dataview_plot_expression_ui(ns("expressionplot"),height=imgH),
            dataview_plot_averagerank_ui(ns("averagerankplot"),height=imgH),
            dataview_plot_tsne_ui(ns("tsneplot"),height=imgH)       
          ),
          shiny::fillRow(
            flex = c(1.5,2), id = "genePlots_row2",
            height = 375, ## width=1600,
            dataview_plot_correlation_ui(ns("correlationplot"),height=imgH),
            dataview_plot_tissue_ui(ns("tissueplot"),height=imgH),
            ),
          tags$div(
            class="caption",
            HTML("<b>Gene plots.</b> <b>(a)</b> Further information about the selected gene from public databases.
            <b>(b)</b> Abundance/expression of selected gene across groups. <b>(c)</b>
            Average rank of the selected gene compared to other genes. <b>(d)</b> t-SNE of samples colored by
            expression of selected gene. <b>(e)</b> Top correlated genes. Darker color corresponds to higher
            expression of the gene. <b>(f)</b> Tissue expression of selected gene.")
          )
        )
      )
    ),

    ##----------------------------------------------------------------------------
    shiny::tabPanel(
      "QC",
      shiny::fillCol(
        flex = c(1,1,NA),
        height = 800,
        bigdash::swappable(
          ns("swap_qc"),
          shiny::fillRow(
            flex = c(1,1,1), id = "counts_tab_row1", height=375,
            bigdash::swappableItem(
              "totalcount_swap",
              dataview_plot_totalcounts_ui(ns("counts_total"), height=imgH),
              ),
            bigdash::swappableItem(
              "boxplot_swap",
              dataview_plot_boxplot_ui(ns("counts_boxplot"), height=imgH),
              ),
            bigdash::swappableItem(
              "histogram_swap",
              dataview_plot_histogram_ui(ns("counts_histplot"), height=imgH)                    
            )
          ),
          shiny::fillRow(
            flex = c(1,1), id = "counts_tab_row2", height=375,
            bigdash::swappableItem(
              "abundance_swap",
              dataview_plot_abundance_ui(ns("counts_abundance"), height=imgH),
              ),
            bigdash::swappableItem(
              "averagecounts_swap",
              dataview_plot_averagecounts_ui(ns("counts_averagecounts"), height=imgH)
            )
          ),
          tags$div(
            class="caption",
            HTML("<b>Counts distribution</b>. Plots associated with the counts, abundance or expression levels across
                the samples/groups.  <b>(a)</b> Total counts per sample or average per group.
                <b>(b)</b> Distribution of total counts per sample/group. The center horizontal bar correspond to
                the median.  <b>(c)</b> Histograms of total counts distribution per sample/group. <b>(d)</b>
                Abundance of major gene types per sample/group. <b>(e)</b> Average count by gene type per sample/group.")
          )
        )
      )
    ),
    ##----------------------------------------------------------------------------
    shiny::tabPanel(
      "Counts",
      shiny::fillCol(
        flex = c(NA,1),
        height = 750,
        tags$div( class="caption",
          HTML("<b>Gene table.</b> The table shows the gene expression values per sample, or average
                expression values across the groups. The column 'rho' reports the correlation with the
                gene selected in 'Search gene' in the left side bar.")
        ),
        dataview_table_rawdata_ui(ns("rawdatatable"))
      )
    ),
    ##----------------------------------------------------------------------------
    shiny::tabPanel(
    "Samples",
    shiny::fillCol(
      flex = c(1.35,1,NA),
      height = 800,
      shiny::fillRow(
        flex = c(2.3,1),
        shiny::div(
          dataview_plot_phenoheatmap_ui(ns("phenoheatmap"),height=imgH),
          style="overflow-y: auto;"),
        dataview_plot_phenoassociation_ui(ns("phenoassociation"),height=imgH)                    
      ),
      dataview_table_samples_ui(ns("sampletable")),
      tags$div(
        class="caption",
        HTML(paste(
          "<b>(a)</b> <b>Phenotype clustering.</b> Clustered heatmap of sample information
           (i.e. phenotype data).","<b>(b)</b> <b>Phenotype association matrix.</b> Clustered
           heatmap of phenotype association. The values corresponds to the -log10(p) value of
           the corresponding statistical test between two phenotype variables. A higher value
           corresponds to stronger 'correlation'.","<b>(c)</b> <b>Sample information table.
           </b> Phenotype information about the samples. Phenotype variables starting with
           a 'dot' (e.g. '.cell cycle' and '.gender' ) have been estimated from the data."
        ))
      ))      
    ),

    ##----------------------------------------------------------------------------
    shiny::tabPanel(
      "Contrasts",
      shiny::fillCol(
        flex = c(NA,1), height = 750,
        tags$div( class="caption",
          HTML("<b>Contrast table.</b> summarizing the contrasts of all comparisons. Non-zero entries
               '+1' and '-1' correspond to the group of interest and control group, respectively. Zero
                or empty entries denote samples not use for that comparison.")
        ),
        dataview_table_contrasts_ui(ns("contrastTable"))
      )),
    ##----------------------------------------------------------------------------
    shiny::tabPanel(
      "Resource info",
      dataview_table_rescources_ui(ns("resources"))
    )

  )
}
