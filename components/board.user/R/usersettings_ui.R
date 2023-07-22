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
          style = htmltools::css(grid_template_columns = "6fr 6fr"),          
          bslib::card(
            ## bslib::card_header(),
            bslib::card_body(
              shiny::h3("Application options"),
              shinyWidgets::prettySwitch(ns("enable_beta"), "Enable beta features"),
              shinyWidgets::prettySwitch(ns("enable_info"), "Show info alerts", value = TRUE),
              selector_switch(
                class = "card-footer-checked",
                label = "show captions",
                is.checked = FALSE
              )
            ),
            bslib::card_footer(NULL)
          ),
          bslib::card(
            ## bslib::card_header(),
            bslib::card_body(
              shiny::h3("New features"),
              shiny::htmlOutput(ns("newfeatures"))
            ),
            bslib::card_footer(NULL)
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
