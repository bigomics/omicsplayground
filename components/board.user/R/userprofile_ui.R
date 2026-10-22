##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

UserProfileUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  OmicsBoardUI(
    ns = ns,
    title = "User Profile",
    header_margin = "0px",
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "User profile",
        bslib::layout_columns(
          height = "100%",
          col_widths = c(4, 8),
          wellPanel(
            shiny::h4("Subscription"),
            uiOutput(ns("plan")),
            shiny::tableOutput(ns("userdata"))
          ),
          bslib::layout_columns(
            col_widths = 12,
            PlotModuleUI(
              ns("usage"),
              plotlib = "plotly",
              download.fmt = c("png", "pdf", "csv", "svg"),
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
