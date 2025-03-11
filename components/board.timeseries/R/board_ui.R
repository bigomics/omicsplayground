##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

TimeSeriesInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace

  settings_taglist <- tagList(
    withTooltip(shiny::selectInput(ns("timevar"), "Time variable:", choices = NULL),
      "Select phenotypes to show in heatmap and phenotype distribution plots.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("module"), "Module:", choices = NULL,
                                   multiple=TRUE),
      "Select module(s) to show in parallel plot.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("contrast"), "Contrast:", choices = NULL),
      "Select contrast to show in table.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("groupvar"), "Color by:", choices = NULL),
      "Select phenotypes to show for creating groups.",
      placement = "top"
    ),
    shiny::br(),
    bslib::accordion(
      id = ns("options_accordion"),
      open = FALSE,
      bslib::accordion_panel(
        "Advanced options",
        icon = icon("cog", lib = "glyphicon"),
        withTooltip(
          shiny::selectInput(
            ns("knn"), "Number of clusters:",
            choices = c(2,3,4,5,7,10,15), selected=7
          ),
          "Choose number of KNN clusters."
        ),
        withTooltip(
          shiny::checkboxInput(ns("timefactor"), "Time as factor", TRUE),
          "Treat time as factor"
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
        "Time clustering",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          bs_alert(parallel_info),
          bslib::layout_columns(
            col_widths = c(6, 6),
            TimeSeriesBoard.clustering_plot_ui(
              id = ns("clustering"),
              title = "Time series clustering"
            ),
            bslib::layout_columns(
              col_widths = 12,
              TimeSeriesBoard.parcoord_plot_ui(
                id = ns("parcoord"),
                title = "Parallel coordinates"
              ),
              TimeSeriesBoard.parcoord_table_ui(
                id = ns("parcoord"),
                title = "Module genes"
              )
            )
          )
        )
      ), ## end of tabpanel
      shiny::tabPanel(
        "Features",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          row_heights = c("auto",1),
          bs_alert(parallel_info),
          bslib::layout_columns(
            col_widths = c(5, 7),
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
      ) ## end of tabpanel
      
    )
  )
}
