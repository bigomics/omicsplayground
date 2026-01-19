##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

AdminPanelUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  div(
    class = "row",
    boardHeader(title = "Admin Panel", info_link = ns("board_info")),
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "Overview",
        bslib::layout_columns(
          height = "calc(100vh - 183px)",
          col_widths = c(4, 8),
          wellPanel(
            shiny::h4("Admin Information"),
            shiny::uiOutput(ns("admin_info")),
            shiny::hr(),
            shiny::h4("System Status"),
            shiny::uiOutput(ns("system_status"))
          ),
          bslib::layout_columns(
            col_widths = 12,
            wellPanel(
              shiny::h4("User Statistics"),
              shiny::tableOutput(ns("user_stats"))
            ),
            br()
          )
        )
      ),
      shiny::tabPanel(
        "User Management"
      ),
      shiny::tabPanel(
        "System Settings"
      )
    )
  )
}
