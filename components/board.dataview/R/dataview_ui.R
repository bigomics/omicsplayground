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
            title = "Gene information",
            info.text = "Information about the selected gene and its function from public databases. For more information, follow the hyperlinks to public databases.",
            caption = "Information about the selected gene and its function from public databases.",
            height = c(600, 800),
            width = c("auto", "100%")
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
                title = "Gene expression",
                info.text = "Samples (or cells) in the barplot can be ungrouped by setting the grouped under the main Options.",
                caption = "Barplot of abundance or expression of grouped samples (or cells) for the gene selected in the Search gene.",
                height = imgH,
                label = "a"
                
              )
            ),
            div(
              class = "col-lg-6 col-xxl-4 col-xxxl-3",
              dataview_plot_averagerank_ui(
                ns("averagerankplot"),
                height = imgH,
                label = "b",
                title = "Average rank",
                info.text = "Select the gene or feature of interest under the main Options.",
                caption = "Ranking of the selected gene by decreasing average expression.",
                width = c("auto", "100%")
              )
            ),
            div(
              class = "col-lg-6 col-xxl-4 col-xxxl-3",
              dataview_plot_tsne_ui(
                ns("tsneplot"),
                height = imgH,
                label = "c",
                title = "t-SNE clustering",
                info.text = "T-SNE clustering of samples (or cells) colored by an expression of the gene selected in the search_gene dropdown menu. The red color represents an over-expression of the selected gene across samples (or cells).",
                caption = "t-SNE of samples colored by expression of selected gene.",
                width = c("auto", "100%")
              )
            ),
            div(
              class = "col-lg-9 col-xxl-7 col-xxxl-5",
              dataview_plot_correlation_ui(
                ns("correlationplot"),
                height = imgH,
                label = "d",
                title = "Top correlated genes",
                info.text = "Colors are from absolute expression levels of genes, where the low and high expressions range between the light and dark colors, respectively.",
                caption = "Barplot of the top positively and negatively correlated genes with the selected gene. Darker color corresponds to higher expression of the gene.",
                width = c("auto", "100%")
              )
            ),
            div(
              class = "col-lg-9 col-xxl-5 col-xxxl-3",
              dataview_plot_tissue_ui(
                ns("tissueplot"),
                height = imgH,
                width = c("auto", "100%"),
                label = "e",
                title = "Tissue expression (GTEX)",
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
            label = "a",
            title = "Total counts",
            info.text = "The samples (or cells) can be grouped/ungrouped in the grouped setting under the main Options.",
            caption = "Barplot of the average number of counts (abundance) for each group.",
            height = imgH,
            width =c("auto", "100%")
          )
        ),
        div(
          class = "col-lg-6 col-xxl-4 col-xxxl-3",
          dataview_plot_boxplot_ui(
            ns("counts_boxplot"),
            title = "Median counts distribution",
            info.text = "The samples (or cells) can be grouped/ungrouped in the grouped setting under the main Options.",
            caption = "Distribution of total counts per sample/group. The center horizontal bar correspond to the median.",
            height = imgH,
            label = "b"
          )
        ),
        div(
          class = "col-lg-6 col-xxl-4 col-xxxl-3",
          dataview_plot_histogram_ui(
            ns("counts_histplot"),
            title = "Density distribution of counts",
            info.text = "The samples (or cells) can be grouped/ungrouped in the grouped setting under the main Options.",
            caption = "Density distribution of total counts per sample/group",
            height = imgH,
            width = c("auto", "100%"),
            label = "c"
          )
        ),
        div(
          class = "col-lg-6 col-xxl-5 col-xxxl-3",
          dataview_plot_genetypes_ui(
            ns("counts_genetypes"),
            title = "Dataset abundance of major gene types",
            info.text = "The samples (or cells) can be grouped/ungrouped in the grouped setting under the main Options. Genetypes can be ribosomal protein genes, kinases or RNA binding motifs, etc..",
            caption = "Barplot showing the dataset relative abundance of counts in terms of major gene types.",
            height = imgH,
            width = c("auto", "100%"),
            label = "d"
          )
        ),
        div(
          class = "col-lg-9 col-xxl-7 col-xxxl-5",
          dataview_plot_abundance_ui(
            ns("counts_abundance"),
            title = "Abundance of major gene types per group",
            info.text = "The samples (or cells) can be grouped/ungrouped in the grouped setting under the main Options. Genetypes can be ribosomal protein genes, kinases or RNA binding motifs, etc..",
            caption = "Barplot showing the group or sample relative abundance of counts in terms of major gene types.",
            height = imgH,
            label = "e",
            width = c("auto", "100%")
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
                title = "Gene expression table",
                info.text = "The column 'rho' reports the correlation with the gene selected in 'Search gene' in the left side bar. If the data type selected is counts, the geometric mean is calculated. The SD column reports the standard deviation of expression across samples (or cells).",
                caption = "The table shows the gene expression values per sample, or average expression values across the groups.",
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
            title = "Phenotype clustering",
            info.text = "Column ordering has been performed using hierarchical clustering on a one-hot encoded matrix.",
            caption = "Clustered heatmap of sample information (i.e. phenotype data)",
            height = imgH,
            width = c("auto", "100%"),
            label = "a"
          ),
          dataview_plot_phenoassociation_ui(
            ns("phenoassociation"),
            height = imgH,
            width = c("auto", "100%"),
            label = "b",
            title = "Phenotype association",
            info.text = "Phenotype association matrix. Clustered heatmap of phenotype association. The values correspond to the -log10(p) value of the corresponding statistical test between two phenotype variables. A higher value corresponds to stronger 'correlation'.",
            caption = "Clustered heatmap of phenotype association."
          )
        ),
        bslib::layout_column_wrap(
          width = 1,
          dataview_table_samples_ui(
            ns("sampletable"),
            height = c(300, TABLE_HEIGHT_MODAL),
            width = c("auto", "100%"),
            title = "Sample information",
            info.text = "Phenotype variables starting with a 'dot' (e.g. '.cell cycle' and '.gender' ) have been estimated from the data.",
            caption = "Phenotype information about the samples." 
          )
        )
      )
    ),

    # contrasts tab #####
    shiny::tabPanel(
      "Contrasts",
      dataview_table_contrasts_ui(
        ns("contrastTable"),
        title = "Contrast table",
        info.text = "Here, you can check which samples belong to which groups for the different comparisons. Non-zero entries '+1' and '-1' correspond to the group of interest and control group, respectively. Zero or empty entries denote samples not use for that comparison.", 
        caption = "Table summarizing the contrasts of all comparisons.",
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
