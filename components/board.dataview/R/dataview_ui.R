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
  ns <- shiny::NS(id) ## namespace

  bigdash::tabSettings(
    #      withTooltip( shiny::actionLink(ns("data_info"), "Tutorial", icon = shiny::icon("youtube")),
    #          "Show more information about this module."),
    shiny::hr(), shiny::br(),
    withTooltip(shiny::selectInput(ns("search_gene"), "Gene:", choices = NULL),
      "Enter a gene of interest for the analysis.",
      placement = "top"
    ),
    withTooltip(
      shiny::selectInput(ns("data_samplefilter"), "Filter samples:",
        choices = NULL, multiple = TRUE
      ),
      "Filter the relevant samples for the analysis.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("data_groupby"), "Group by:", choices = NULL),
      "Select phenotype for grouping the samples.",
      placement = "top"
    ),
    shiny::br(),
    withTooltip(shiny::actionLink(ns("data_options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(), shiny::br(),
    shiny::conditionalPanel(
      "input.data_options % 2 == 1",
      ns = ns,
      withTooltip(
        shiny::radioButtons(ns("data_type"), "Data type:",
          choices = c("counts", "logCPM"), selected = "logCPM", inline = TRUE
        ),
        "Choose an input data type for the analysis.",
        placement = "bottom"
      )
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
  ns <- shiny::NS(id) ## namespace

  imgH <- c(330, 600) ## heights for small and fullscreen image
  imgH <- c("35vh", "70vh") ## heights for small and fullscreen image

  tabs <- shiny::tabsetPanel(
    id = ns("tabs"),


    # QC tab #####

    shiny::tabPanel(
      "Sample QC",
      div(
        class = "row",
        div(
          class = "col-lg-6 col-xxl-4 col-xxxl-3",
          dataview_plot_totalcounts_ui(
            ns("counts_total"),
            height = imgH,
            label = "a"
          )
        ),
        div(
          class = "col-lg-6 col-xxl-4 col-xxxl-3",
          dataview_plot_boxplot_ui(
            ns("counts_boxplot"),
            height = imgH,
            label = "b"
          )
        ),
        div(
          class = "col-lg-6 col-xxl-4 col-xxxl-3",
          dataview_plot_histogram_ui(
            ns("counts_histplot"),
            height = imgH,
            label = "c"
          )
        ),
        div(
          class = "col-lg-6 col-xxl-5 col-xxxl-3",
          dataview_plot_genetypes_ui(
            ns("counts_genetypes"),
            height = imgH, label = "d"
          )
        ),
        div(
          class = "col-lg-9 col-xxl-7 col-xxxl-5",
          dataview_plot_abundance_ui(
            ns("counts_abundance"),
            height = imgH, label = "e"
          )
        )
      ),
      tags$div(
        class = "caption",
        HTML("<b>Counts distribution</b>. Plots associated with the counts, abundance or expression levels across
            the samples/groups.  <b>(a)</b> Total counts per sample or average per group.
            <b>(b)</b> Distribution of total counts per sample/group. The center horizontal bar correspond to
            the median.  <b>(c)</b> Histograms of total counts distribution per sample/group. <b>(d)</b>
            Abundance of major gene types per sample/group. <b>(e)</b> Average count by gene type per sample/group.")
      )
    ),

    # Gene overview tab #####
    shiny::tabPanel(
      "Gene overview",
      div(
        class = "row",
        div(
          class = "col-md-2",
          dataview_module_geneinfo_ui(ns("geneinfo")),
        ),
        div(
          class = "col-md-10",
          div(
            class = "row",
            div(
              class = "col-lg-6 col-xxl-4 col-xxxl-3",
              dataview_plot_expression_ui(
                ns("expressionplot"),
                height = imgH,
                label = "a"
              )
            ),
            div(
              class = "col-lg-6 col-xxl-4 col-xxxl-3",
              dataview_plot_averagerank_ui(
                ns("averagerankplot"),
                height = imgH,
                label = "b"
              )
            ),
            div(
              class = "col-lg-6 col-xxl-4 col-xxxl-3",
              dataview_plot_tsne_ui(
                ns("tsneplot"),
                height = imgH,
                label = "c"
              )
            ),
#          ),
#          div(
#            class = "row",
            div(
              class = "col-lg-9 col-xxl-7 col-xxxl-5",
              dataview_plot_correlation_ui(
                ns("correlationplot"),
                height = imgH,
                label = "d"
              )
            ),
            div(
              class = "col-lg-9 col-xxl-5 col-xxxl-3",
              dataview_plot_tissue_ui(
                ns("tissueplot"),
                height = imgH,
                label = "e"
              )
            )
          ),
          tags$div(
            class = "caption",
            HTML("<b>Gene plots.</b> <b>(a)</b> Further information about the selected gene from public databases.
            <b>(b)</b> Abundance/expression of selected gene across groups. <b>(c)</b>
            Average rank of the selected gene compared to other genes. <b>(d)</b> t-SNE of samples colored by
            expression of selected gene. <b>(e)</b> Top correlated genes. Darker color corresponds to higher
            expression of the gene. <b>(f)</b> Tissue expression of selected gene.")
          )
        )
      )
    ),

    # counts table tab #####

    shiny::tabPanel(
      "Counts table",
      dataview_table_rawdata_ui(
        ns("rawdatatable"),
        height = c("75vh", TABLE_HEIGHT_MODAL),
        width = c("100%", "90%")
      ),
      tags$div(
        class = "caption",
        HTML("<b>Gene table.</b> The table shows the gene expression values per sample, or average
              expression values across the groups. The column 'rho' reports the correlation with the
              gene selected in 'Search gene' in the left side bar. If the data type selected is counts,
              the geometric mean is calculated.")
      )
    ),
    # Sample information #####
    shiny::tabPanel(
      "Sample information",
      div(
        class = "row",
        div(
          class = "col-md-6",
          shiny::div(
            dataview_plot_phenoheatmap_ui(
              ns("phenoheatmap"),
              height = imgH,
              label = "a"
            ),
            style = "overflow-y: auto;"
          )
        ),
        div(
          class = "col-md-6",
          dataview_plot_phenoassociation_ui(
            ns("phenoassociation"),
            height = imgH,
            label = "b"
          )
        )
      ),
      dataview_table_samples_ui(
        ns("sampletable"),
        height = c(280, TABLE_HEIGHT_MODAL),
        width = c("auto", "90%")
      ),
      tags$div(
        class = "caption",
        HTML(
          "<b>(a)</b> <b>Phenotype clustering.</b> Clustered heatmap of sample information
            (i.e. phenotype data).", "<b>(b)</b> <b>Phenotype association matrix.</b> Clustered
            heatmap of phenotype association. The values corresponds to the -log10(p) value of
            the corresponding statistical test between two phenotype variables. A higher value
            corresponds to stronger 'correlation'.", "<b>(c)</b> <b>Sample information table.
            </b> Phenotype information about the samples. Phenotype variables starting with
            a 'dot' (e.g. '.cell cycle' and '.gender' ) have been estimated from the data."
        )
      )
    ),


    # contrasts tab #####

    shiny::tabPanel(
      "Contrasts",
      dataview_table_contrasts_ui(
        ns("contrastTable"),
        height = c(500, TABLE_HEIGHT_MODAL),
        width = c("auto", "90%")
      ),
      tags$div(
        class = "caption",
        HTML(
          "<b>Contrast table.</b> summarizing the contrasts of all comparisons. Non-zero entries
          '+1' and '-1' correspond to the group of interest and control group, respectively. Zero
          or empty entries denote samples not use for that comparison."
        )
      )
    ),

    # Resource info #####

    shiny::tabPanel(
      "Resource info",
      dataview_table_rescources_ui(ns("resources"))
    )
  )

  div(
    boardHeader(title = "Data View", info_link = ns("data_info")),
    tabs
  )
}
