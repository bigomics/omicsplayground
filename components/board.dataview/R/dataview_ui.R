##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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


    # Gene overview tab #####
    shiny::tabPanel(
      "Gene overview",
      div(
        class = "row",
        div(
          class = "col-md-2",
          dataview_module_geneinfo_ui(
            id = ns("geneinfo")
          ),
        ),
        div(
          class = "col-md-10",
          div(
            class = "row",
            div(
              class = "col-lg-6 col-xxl-4 col-xxxl-3",
              dataview_plot_expression_ui(
                id = ns("expressionplot"),
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
          )
        )
      )
    ),
    
    # QC tab #####
    shiny::tabPanel(
      "Sample QC",
      shinyjqui::jqui_sortable(
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
      )),
      tags$div(
        class = "caption",
        HTML("<b>Counts distribution</b>. Plots associated with the counts, abundance or expression levels across
            the samples/groups.  <b>(a)</b> Total counts per sample or average per group.
            <b>(b)</b> Distribution of total counts per sample/group. The center horizontal bar correspond to
            the median.  <b>(c)</b> Histograms of total counts distribution per sample/group. <b>(d)</b>
            Abundance of major gene types per sample/group. <b>(e)</b> Average count by gene type per sample/group.")
      )
    ),

    # counts table tab #####
    shiny::tabPanel(
        "Counts table",
        bslib::layout_column_wrap(
            width = 1,
            dataview_table_rawdata_ui(
                ns("rawdatatable"),
                height = c("75vh", TABLE_HEIGHT_MODAL),
                width = c("100%", "100%")
            )
        )
    ),
    # Sample information #####
    shiny::tabPanel(
      "Sample information",
      bslib::layout_column_wrap(
        width = 1,
        bslib::layout_column_wrap(
          width = 1/2,
          heights_equal = "row",
          dataview_plot_phenoheatmap_ui(
            ns("phenoheatmap"),
            height = imgH,
            label = "a"
          ),
          dataview_plot_phenoassociation_ui(
            ns("phenoassociation"),
            height = imgH,
            label = "b"
          )
        ),
        bslib::layout_column_wrap(
          width = 1,
          dataview_table_samples_ui(
            ns("sampletable"),
            height = c(300, TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    ),


    # contrasts tab #####
    shiny::tabPanel(
      "Contrasts",
      dataview_table_contrasts_ui(
        ns("contrastTable"),
        height = c(500, TABLE_HEIGHT_MODAL),
        width = c("auto", "100%")
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
    div(boardHeader(title = "Data View", info_link = ns("board_info")),
    ),
    tabs
  )
}
