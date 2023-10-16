##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

UserSettingsInputs <- function(id) {
  ns <- shiny::NS(id)
  bigdash::tabSettings()
}

UserSettingsUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  div(
    boardHeader(title = "App Settings", info_link = ns("board_info")),
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "App Settings & News",
        bslib::layout_column_wrap(
          height = "calc(100vh - 183px)",
          width = 1,
          style = htmltools::css(grid_template_columns = "2fr 5fr 5fr"),
          CardUI(
            bslib::input_switch(ns("enable_beta"), "Enable beta features"),
            bslib::input_switch(ns("enable_info"), "Show info alerts", value = TRUE),
            selector_switch(
              class = "card-footer-checked",
              label = "show captions",
              is.checked = FALSE
            ),
            title = "Application options"
          ),
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
        bslib::layout_column_wrap(
          width = 1,
          height = "calc(100vh - 183px)",
          user_table_resources_ui(ns("resources"))
        )
      )
    )
  )
}
