##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## Append one row per admin action to the per-host audit log.
## Hostname is part of the filename because multiple deploys may share
## the same mounted etc/ folder; admin events must not be mixed across
## servers or the audit trail becomes ambiguous.
log_admin_action <- function(admin_email, action, subjects,
                             source_labels = "", destination = "") {
  if (length(subjects) == 0) return(invisible())
  host <- NULL
  if (exists("opt", inherits = TRUE)) host <- opt$HOSTNAME
  if (is.null(host) || !nzchar(host)) {
    host <- toupper(Sys.info()[["nodename"]])
  }
  host <- gsub("[^A-Za-z0-9._-]", "_", host)
  log.file <- file.path(ETC, paste0("PGXADMIN-", host, ".log"))
  log.entry <- data.frame(
    date = format(Sys.time(), tz = "CET"),
    admin = admin_email,
    action = action,
    subject = subjects,
    source = source_labels,
    destination = destination,
    stringsAsFactors = FALSE
  )
  tryCatch({
    if (file.exists(log.file)) {
      write.table(log.entry, file = log.file, col.names = FALSE,
                  row.names = FALSE, sep = ",", append = TRUE)
    } else {
      write.table(log.entry, file = log.file, col.names = TRUE,
                  row.names = FALSE, sep = ",")
    }
  }, error = function(e) dbg("[admin_log] write error: ", e$message))
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

    ## ================================================================================
    ## =================================== END ========================================
    ## ================================================================================
  })
}
