##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' AdminPanel module UI Input function
#'
#' @description A shiny Module. Renders the input parts (sidebar contents) for the module.
#'
#' @param id Internal parameters for {shiny}.
#'
#' @export
AdminPanelInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace

  bigdash::tabSettings(
    shiny::br()
  )
}

#' AdminPanel module UI output function
#'
#' @description Renders the output part for the module as tabsetPanel object
#'
#' @param id Internal parameters for {shiny}.
#'
#' @export
AdminPanelUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- "calc(100vh - 181px)"

  tabs <- shiny::tabsetPanel(
    id = ns("tabs1"),
    shiny::tabPanel(
      "Overview",
      bslib::layout_columns(
        height = fullH,
        col_widths = c(4, 8),
        row_heights = list("auto", 1),
        bslib::layout_columns(
          col_widths = 12,
          admin_module_info_ui(
            id = ns("admin_info"),
            title = "Admin Information",
            info.text = "Information about the currently logged-in administrator.",
            caption = "Details about the current admin user.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          ),
          admin_module_status_ui(
            id = ns("system_status"),
            title = "System Status",
            info.text = "Current system status information including R version and platform details.",
            caption = "System environment information.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        ),
        bslib::layout_columns(
          col_widths = 12,
          admin_table_users_ui(
            id = ns("user_stats"),
            title = "User Statistics",
            info.text = "Summary statistics about registered users on the platform.",
            caption = "Table of user statistics and platform metrics.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    ),
    shiny::tabPanel(
      "User Management",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        row_heights = list("auto", 1),
        bs_alert("Manage user credentials. Click on a cell to edit it. Use the buttons below to add or remove users, and save changes when done."),
        admin_table_credentials_ui(
          id = ns("credentials"),
          title = "User Credentials",
          info.text = "Editable table of user credentials from the CREDENTIALS file.",
          caption = "Edit user credentials and save changes.",
          height = c("100%", TABLE_HEIGHT_MODAL)
        )
      )
    ),
    shiny::tabPanel(
      "System Settings"
    )
  )

  div(
    boardHeader(title = "Admin Panel", info_link = ns("board_info")),
    tabs
  )
}
