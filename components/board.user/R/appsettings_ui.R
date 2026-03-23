##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

AppSettingsInputs <- function(id) {
  ns <- shiny::NS(id)
  bigdash::tabSettings()
}

AppSettingsUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  div(
    boardHeader(title = "App Settings", info_link = ns("board_info")),
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "App Settings & News",
        bslib::layout_columns(
          height = "calc(100vh - 183px)",
          col_widths = c(6, 6),
          PlotModuleUI(
            ns("newfeatures"),
            outputFunc = htmlOutput,
            info.text = "New features and changes",
            title = "New features"
          ),
          TableModuleUI(
            ns("packages"),
            info.text = "Packages and versions used in Omics Playground.",
            title = "Package versions"
          )
        )
      ),
      # Resource info #####
      shiny::tabPanel(
        "Resource info",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 183px)",
          user_table_resources_ui(ns("resources"))
        )
      )
    )
  )
}

AppSettingsUI2 <- function(id) {
  ns <- shiny::NS(id) ## namespace

  div(
    boardHeader(title = "Settings", info_link = ns("board_info")),

    bslib::navset_pill_list(
      id = ns("tabs1"),
      widths = c(2,10),
      bslib::nav_panel(
        "New features",
        PlotModuleUI(
          ns("newfeatures"),
          outputFunc = htmlOutput,
          info.text = "New features and changes",
          title = "New features"
        )
      ),

      bslib::nav_panel(
        "Package versions",
        TableModuleUI(
          ns("packages"),
          info.text = "Packages and versions used in Omics Playground.",
          title = "Package versions"
        )
      ),

      bslib::nav_panel(
        "Resource info",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 183px)",
          user_table_resources_ui(ns("resources"))
        )
      )

    )
  )
}
