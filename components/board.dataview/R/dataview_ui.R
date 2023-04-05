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
            id = ns("geneinfo"),
            info.text = "Information about the selected gene and its function from public databases. For more information, follow the hyperlinks to public databases.",
            caption = "Information about the selected gene and its function from public databases."
            
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
                label = "a",
                info.text = "Samples (or cells) in the barplot can be ungrouped by setting the grouped under the main Options.",
                caption = "Barplot of abundance or expression of grouped samples (or cells) for the gene selected in the Search gene."
              )
            ),
            div(
              class = "col-lg-6 col-xxl-4 col-xxxl-3",
              dataview_plot_averagerank_ui(
                ns("averagerankplot"),
                height = imgH,
                label = "b",
                info.text = "Select the gene or feature of interest under the main Options.",
                caption = "Ranking of the selected gene by decreasing average expression."

              )
            ),
            div(
              class = "col-lg-6 col-xxl-4 col-xxxl-3",
              dataview_plot_tsne_ui(
                ns("tsneplot"),
                height = imgH,
                label = "c",
                info.text = "T-SNE clustering of samples (or cells) colored by an expression of the gene selected in the search_gene dropdown menu. The red color represents an over-expression of the selected gene across samples (or cells).",
                caption = "t-SNE of samples colored by expression of selected gene."
              )
            ),
            div(
              class = "col-lg-9 col-xxl-7 col-xxxl-5",
              dataview_plot_correlation_ui(
                ns("correlationplot"),
                height = imgH,
                label = "d",
                info.text = "Colors are from absolute expression levels of genes, where the low and high expressions range between the light and dark colors, respectively.",
                caption = "Barplot of the top positively and negatively correlated genes with the selected gene. Darker color corresponds to higher expression of the gene."
              )
            ),
            div(
              class = "col-lg-9 col-xxl-5 col-xxxl-3",
              dataview_plot_tissue_ui(
                ns("tissueplot"),
                height = imgH,
                label = "e",
                info.text = "Colors correspond to 'tissue clusters' as computed by unsupervised clustering. Select the gene or feature of interest under the main Options.",
                caption = paste("Top 15 expressing tissues for the selected gene in the tissue expression GTEx database. Colors represent tissue clusters.")
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
