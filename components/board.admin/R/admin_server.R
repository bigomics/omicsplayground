##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' AdminPanel module server function
#'
#' @description A shiny Module (server code).
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @param auth Reactive list that provides authentication info (ADMIN, username, email)
#' @param credentials_file Path to the CREDENTIALS CSV file (optional)
#'
#' @export
AdminPanelBoard <- function(id, auth, credentials_file = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    dbg("[AdminPanelBoard] >>> initializing AdminBoard...")

    ## Check if user is admin - this is a critical security check
    is_admin <- reactive({
      isTRUE(auth$ADMIN)
    })

    ## ----------------------------------------------------------------------
    ## More Info (pop up window)
    ## ----------------------------------------------------------------------

    shiny::observeEvent(input$board_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Admin Panel</strong>"),
        shiny::HTML(
          "The Admin Panel provides administrative functions for managing
          the OmicsPlayground platform. This panel is only accessible to
          users with admin privileges."
        ),
        easyClose = TRUE, size = "l"
      ))
    })

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    admin_module_info_server(
      "admin_info",
      auth = auth
    )

    admin_module_status_server(
      "system_status",
      auth = auth
    )

    ## ================================================================================
    ## ===============================  TABLES ========================================
    ## ================================================================================

    admin_table_users_server(
      "user_stats",
      auth = auth
    )

    admin_table_credentials_server(
      "credentials",
      auth = auth,
      credentials_file = credentials_file
    )

    ## ================================================================================
    ## =================================== END ========================================
    ## ================================================================================
  })
}
