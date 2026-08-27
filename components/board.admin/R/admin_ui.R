##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
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

  fullH <- "100%"

  ## Both choosers list every board of the full sidebar menu, labelled with its
  ## group so "Samples" and "Features" are not ambiguous.
  menu_tree <- opg_menu_tree()
  basic_choices <- unlist(lapply(names(menu_tree), function(g) {
    stats::setNames(names(menu_tree[[g]]), paste0(g, " - ", menu_tree[[g]]))
  }))

  ## same shell as the Settings board (appsettings_ui.R): nav list on the left
  tabs <- bslib::navset_pill_list(
    id = ns("tabs1"),
    widths = c(2, 10),
    bslib::nav_panel(
      "Overview",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        row_heights = list(1),
        admin_table_users_ui(
          id = ns("user_stats"),
          title = "Platform Statistics",
          info.text = "Summary statistics about registered users, datasets, and storage on the platform.",
          caption = "Platform overview metrics.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    ),
    bslib::nav_panel(
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
    bslib::nav_panel(
      "Basic Menu",
      bslib::layout_columns(
        col_widths = c(6, 6, 12),
        ## bound the height so the two lists scroll inside their cards and the
        ## Save row stays on screen instead of running off the bottom
        height = fullH,
        row_heights = list(1, "auto"),
        bslib::card(
          bslib::card_header("Basic menu items"),
          bslib::card_body(
            bs_alert("Choose which menu items are kept when a user switches on 'Basic menu' in Settings. The choice applies to the whole deployment; users see it after a page reload."),
            shiny::checkboxGroupInput(
              ns("basic_menu"),
              label = NULL,
              choices = basic_choices,
              selected = opt$BASIC_MENU
            )
          )
        ),
        bslib::card(
          bslib::card_header("Locked settings"),
          bslib::card_body(
            bs_alert("Choose which boards have their advanced settings greyed out in the settings sidebar. Unticked boards keep their options fully usable in basic mode. Boards with no advanced settings are unaffected either way."),
            shiny::checkboxGroupInput(
              ns("basic_locked"),
              label = NULL,
              choices = basic_choices,
              selected = opt$BASIC_LOCKED
            )
          )
        ),
        div(
          shiny::textOutput(ns("basic_menu_status")),
          shiny::actionButton(ns("save_basic_menu"), "Save", class = "btn-primary")
        )
      )
    ),
    bslib::nav_panel(
      "Data Management",
      bslib::layout_columns(
        col_widths = 12,
        height = fullH,
        row_heights = list("auto", 1),
        bs_alert("Browse, move or delete datasets across users, shared and public directories."),
        admin_table_datamanager_ui(
          id = ns("datamanager"),
          title = "Data Management",
          info.text = "Manage .pgx files across user folders, data_shared, and data_public directories.",
          caption = "Select files and use the action buttons to move, copy, or delete.",
          height = c("100%", TABLE_HEIGHT_MODAL)
        )
      )
    )
  )

  OmicsBoardUI(
    ns = ns,
    title = "Admin Panel",
    tabs
  )
}
