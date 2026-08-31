##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


## Append one row per admin action to the per-host audit log.
## Hostname is part of the filename because multiple deploys may share
## the same mounted etc/ folder; admin events must not be mixed across
## servers or the audit trail becomes ambiguous.
log_admin_action <- function(admin_email, action, subjects,
                             source_labels = "", destination = "") {
  if (length(subjects) == 0) {
    return(invisible())
  }
  log.file <- file.path(ETC, paste0("PGXADMIN-", opg_hostname(), ".log"))
  log.entry <- data.frame(
    date = format(Sys.time(), tz = "CET"),
    admin = admin_email,
    action = action,
    subject = subjects,
    source = source_labels,
    destination = destination,
    stringsAsFactors = FALSE
  )
  tryCatch(
    {
      if (file.exists(log.file)) {
        write.table(log.entry,
          file = log.file, col.names = FALSE,
          row.names = FALSE, sep = ",", append = TRUE
        )
      } else {
        write.table(log.entry,
          file = log.file, col.names = TRUE,
          row.names = FALSE, sep = ","
        )
      }
    },
    error = function(e) dbg("[admin_log] write error: ", e$message)
  )
}


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

    OmicsBoard(
      session, pgx = NULL, title = "Admin Panel",
      infotext = "The Admin Panel provides administrative functions for managing
          the OmicsPlayground platform. This panel is only accessible to
          users with admin privileges."
    )

    ## Check if user is admin - this is a critical security check
    is_admin <- reactive({
      isTRUE(auth$ADMIN)
    })

    ## ================================================================================
    ## ===============================  TABLES ========================================
    ## ================================================================================

    admin_table_users_server(
      "user_stats",
      auth = auth,
      credentials_file = credentials_file
    )

    admin_table_credentials_server(
      "credentials",
      auth = auth,
      credentials_file = credentials_file
    )

    admin_table_datamanager_server(
      "datamanager",
      auth = auth
    )

    ## ----------------------------------------------------------------------
    ## Basic menu chooser. Persisted per deploy in etc/BASIC_MENU-<HOSTNAME>
    ## and read back by global.R at startup. `opt` is per R process (ShinyProxy
    ## gives each user their own), so other sessions pick this up on reload.
    ## ----------------------------------------------------------------------

    basic_menu_status <- shiny::reactiveVal("")
    output$basic_menu_status <- shiny::renderText(basic_menu_status())

    ## Live preview before Save: as soon as a checkbox changes, apply the
    ## selection to the admin's own session (via the same reactiveVals the
    ## Save handler pushes into), so the admin can check how the basic-mode
    ## UI looks before persisting anything. Session-local only: it vanishes
    ## on reload, and only Save writes etc/ and opt for other sessions.
    ## Both effects are visible while the admin's own session has Basic menu
    ## on (menu filtering and the body.basic-mode greying both key off it).
    shiny::observe({
      if (is.null(input$basic_menu)) return() ## panel not rendered yet
      ## Unticking *all* menu items previews the runtime default fallback
      ## (opg_basic_menu_boards() does); Save itself still rejects it.
      boards_rv <- getUserOption(session, "basic_menu_boards")
      if (!is.null(boards_rv)) {
        boards_rv(opg_basic_menu_boards(boards = input$basic_menu))
      }
      locked_rv <- getUserOption(session, "locked_boards")
      if (!is.null(locked_rv)) {
        locked_rv(input$basic_locked %||% character(0))
      }
    })

    shiny::observeEvent(input$save_basic_menu, {
      shiny::req(isTRUE(auth$ADMIN))
      boards <- input$basic_menu
      locked <- input$basic_locked %||% character(0) ## unticking all is allowed
      if (!length(boards)) {
        basic_menu_status("Select at least one menu item.")
        return()
      }
      tryCatch(
        {
          writeLines(boards, etc_host_file("BASIC_MENU"))
          writeLines(locked, etc_host_file("BASIC_LOCKED"))
          opt$BASIC_MENU <<- boards
          opt$BASIC_LOCKED <<- locked
          ## Live-refilter the admin's own sidebar now, without a reload --
          ## other sessions still pick this up on their next page load, since
          ## opt is per-process (see note above).
          boards_rv <- getUserOption(session, "basic_menu_boards")
          if (!is.null(boards_rv)) {
            boards_rv(opg_basic_menu_boards(boards = boards))
          }
          ## Same live treatment for the locked-settings selection: the
          ## locked_boards observer in opg_server.R re-toggles the class in
          ## the admin's own settings sidebar.
          locked_rv <- getUserOption(session, "locked_boards")
          if (!is.null(locked_rv)) locked_rv(locked)
          log_admin_action(
            admin_email = auth$email,
            action = "basic_menu",
            subjects = paste(
              c(boards, if (length(locked)) paste0("locked:", locked) else "locked:none"),
              collapse = ";"
            )
          )
          basic_menu_status(paste0(
            "Saved at ", format(Sys.time(), "%H:%M:%S"),
            " - your menu and locked settings update now; other users see it after a page reload."
          ))
        },
        error = function(e) basic_menu_status(paste("Save failed:", e$message))
      )
    })

    ## ================================================================================
    ## =================================== END ========================================
    ## ================================================================================
  })
}
