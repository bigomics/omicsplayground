##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

TimeSeriesInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace

  settings_taglist <- tagList(
    withTooltip(shiny::selectInput(ns("module"), "Select module:", choices = NULL,
      multiple=FALSE),
      "Select module to show."
    ),
    withTooltip(shiny::selectInput(ns("contrast"), "Contrast:", choices = NULL),
      "Select contrast to show in table."
    ),
    shiny::br(),
    bslib::accordion(
      id = ns("options_accordion"),
      open = FALSE,
      bslib::accordion_panel(
        "Advanced options",
        icon = icon("cog", lib = "glyphicon"),
        withTooltip(shiny::selectInput(
          ns("timevar"), "Time variable:", choices = NULL),
          "Select phenotypes to show in heatmap and phenotype distribution plots."
        ),
        withTooltip(
          shiny::checkboxInput(ns("timefactor"),
            "Time as factor",
            TRUE
          ),
          "Treat time as factor"
        ),
        withTooltip(
          shiny::selectInput(
            ns("knn"),
            "Number of modules:",
            choices =  c(4,6,9,12),
            selected = 9
          ),
          "Choose number of KNN clusters."
        ),
        withTooltip(        
          shiny::selectInput(
            ns("maxfeatures"),
            "Max features per module:",
            choices=c(50,100,200,500),
            selected=100
          ),
          "Set maximum features to show"
        ),
        withTooltip(        
          shiny::checkboxInput(ns("filtermodules"), "Filter modules", FALSE),
          "Apply eigenmodule filter"
        )
      )
    )
  )

  bigdash::tabSettings(
    settings_taglist
  )
}


TimeSeriesUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  board_info <- "The TimeSeries Board performs unsupervised clustering analysis. After having done the QC, it is probably the first way to explore your data. The main purpose is to discover patterns and subgroups in the data, show correlation with known phenotypes, detect outliers, or investigate batch effects."

  parallel_info <- HTML("<b>Time clustering</b> groups features with similar variation in time together using KNN. We can see different trends in the expression of groups of genes, or so-called gene modules. The parallel plot is interactive so you can manually order the time points.")

  div(
    boardHeader(title = "Time Series", info_link = ns("board_info")),
    shiny::tabsetPanel(
      id = ns("tabs1"),

      shiny::tabPanel(
        "Statistics",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          row_heights = c("auto",1),
          bs_alert(parallel_info),
          bslib::layout_columns(
            height = "calc(100vh - 181px)",
            col_widths = c(6,6),
            TimeSeriesBoard.features_table(
              id = ns("features"),
              title = "Feature table"
            ),
            TimeSeriesBoard.features_plot(
              id = ns("features"),
              title = "Selected features"
            )
          )
        )
      ), ## end of tabpanel

      shiny::tabPanel(
        "Clustering",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          bs_alert(parallel_info),
          bslib::layout_columns(
            height = "calc(100vh - 181px)",
            col_widths = c(6, 6),
            bslib::layout_columns(
              col_widths = 12,
              TimeSeriesBoard.clustering_plot_ui(
                id = ns("clustering"),
                title = "Time series clustering"
              ),
              TimeSeriesBoard.parcoord_table_ui(
                id = ns("parcoord"),
                title = "Genes in module"
              )
            ),
            bslib::layout_columns(
              col_widths = 12,
              TimeSeriesBoard.enrichment_lolliplot_ui(
                id = ns("enrichment"),
                title = "Enrichment plot"
              ),
              TimeSeriesBoard.enrichment_table_ui(
                id = ns("enrichment"),
                title = "Enriched genesets"
              )
            )
          )
        ),
        bslib::layout_columns(
          col_widths = c(6,6),
          height = "45vh",
          TimeSeriesBoard.parcoord_plot_ui(
            id = ns("parcoord"),
            title = "Parallel coordinates"
          )
       )
      ) ## end of tabpanel
      
    )
  )
}
