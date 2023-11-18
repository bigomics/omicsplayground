##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

UserProfileUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  div(
    class = "row",
    boardHeader(title = "User Profile", info_link = ns("board_info")),
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "User profile",
        bslib::layout_column_wrap(
          height = "calc(100vh - 183px)",
          width = 1,
          style = htmltools::css(grid_template_columns = "4fr 8fr"),
          wellPanel(
            shiny::h4("Subscription"),
            uiOutput(ns("plan")),
            shiny::tableOutput(ns("userdata"))
          ),
          bslib::layout_column_wrap(
            width = 1,
            ## shiny::plotOutput(ns("usage")),
            PlotModuleUI(
              ns("usage"),
              plotlib = "plotly",
              download.fmt = c("png", "pdf", "csv"),
              title = "Platform usage",
              height = c("100%", TABLE_HEIGHT_MODAL)
            ),
            br()
          )
        )
      )
    )
  )
}
