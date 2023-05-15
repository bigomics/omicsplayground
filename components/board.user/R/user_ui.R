##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

UserInputs <- function(id) {
  ns <- shiny::NS(id)
  bigdash::tabSettings()
}

UserUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  div(
      class = "row",
      boardHeader(title = "Profile", info_link = ns("board_info")),
      div(
          class = "col-md-7",
          shiny::tabsetPanel(
              id = ns("tabs1"),
              shiny::tabPanel(
                  "App Settings",
                  bslib::layout_column_wrap(
                      height = "calc(100vh - 183px)",
                      width = 1,
                      tagList(
                          shinyWidgets::prettySwitch(ns("enable_beta"), "enable beta features"),
                          shinyWidgets::prettySwitch(ns("enable_tabinfo"), "enable tab info / alerts")
                      )
                  )
              ),
              shiny::tabPanel(
                  "Subscription",
                  bslib::layout_column_wrap(
                      height = "calc(100vh - 183px)",
                      width = 1,
                      tagList(
                          shiny::h4("Personal"),
                          uiOutput(ns("plan")),
                          shiny::tableOutput(ns("userdata"))
                      )
                  )
              )
          )
      ),
      div(
          class = "col-md-5",
          shiny::tabsetPanel(
              id = ns("tabs2"),
              shiny::tabPanel(
                  "News",
                  bslib::layout_column_wrap(
                      height = "calc(100vh - 183px)",
                      width = 1,
                      shiny::htmlOutput(ns("news"))
                  )
              )
          )
      )
  )

}
