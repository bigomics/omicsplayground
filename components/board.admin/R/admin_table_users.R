##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


admin_table_users_ui <- function(
  id,
  title,
  height,
  width = c("auto", "100%"),
  caption,
  info.text
) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("tbl"),
    width = width,
    height = height,
    title = title,
    info.text = info.text,
    caption = caption
  )
}

admin_table_users_server <- function(id, auth,
                                     credentials_file = NULL,
                                     scrollY = "calc(100vh - (240px + 140px))") {
  moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      shiny::req(isTRUE(auth$ADMIN))
      dbg("[AdminBoard::user_stats] table_data")

      ## User directories
      user_dirs <- list.dirs(PGX.DIR, full.names = FALSE, recursive = FALSE)
      user_dirs <- grep("@", user_dirs, value = TRUE)
      n_user_dirs <- length(user_dirs)

      ## Count .pgx files across all user directories
      all_pgx <- list.files(PGX.DIR, pattern = "\\.pgx$", recursive = TRUE)
      n_pgx_total <- length(all_pgx)

      ## Shared datasets
      n_shared <- 0
      if (dir.exists(SHARE.DIR)) {
        n_shared <- length(list.files(SHARE.DIR, pattern = "\\.pgx$", recursive = TRUE))
      }

      ## Public datasets
      n_public <- 0
      if (dir.exists(PUBLIC.DIR)) {
        n_public <- length(list.files(PUBLIC.DIR, pattern = "\\.pgx$", recursive = TRUE))
      }

      ## Registered users from credentials file
      n_registered <- NA
      if (!is.null(credentials_file) && file.exists(credentials_file)) {
        cred <- tryCatch(
          read.csv(credentials_file, colClasses = "character", stringsAsFactors = FALSE),
          error = function(e) NULL
        )
        if (!is.null(cred)) n_registered <- nrow(cred)
      }

      ## Build stats table
      metrics <- c(
        "Registered users",
        "User directories",
        "User datasets (.pgx)",
        "Shared datasets",
        "Public datasets"
      )

      values <- c(
        ifelse(is.na(n_registered), "N/A", as.character(n_registered)),
        as.character(n_user_dirs),
        as.character(n_pgx_total),
        as.character(n_shared),
        as.character(n_public)
      )

      dt <- data.frame(
        Metric = metrics,
        Value = values,
        check.names = FALSE
      )
      dt
    })

    table.RENDER <- function() {
      dt <- table_data()
      req(dt)
      DT::datatable(dt,
        class = "compact hover",
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = 1),
        options = list(
          dom = "lfrtip",
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = scrollY,
          scrollResize = TRUE,
          deferRender = TRUE
        )
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    modal_table.RENDER <- function() {
      dt <- table_data()
      req(dt)
      DT::datatable(dt,
        class = "compact hover",
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        options = list(
          dom = "lfrtip",
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = SCROLLY_MODAL,
          deferRender = TRUE
        )
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "20px", lineHeight = "70%")
    }

    TableModuleServer(
      "tbl",
      func = table.RENDER,
      func2 = modal_table.RENDER,
      selector = "none"
    )
  }) ## end of moduleServer
} ## end of server
