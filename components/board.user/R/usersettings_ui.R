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
        class = "row",
        boardHeader(title = "Settings", info_link = ns("board_info")),
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
                            shinyWidgets::prettySwitch(ns("enable_beta"), "Enable beta features"),
                            shinyWidgets::prettySwitch(ns("enable_info"), "Show info alerts", value = TRUE)
                        )
                    )
                )
            )
        )
    )

}
