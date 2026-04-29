##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


admin_table_datamanager_ui <- function(
  id,
  title = "Data Management",
  height = c("100%", 800),
  info.text = "",
  caption = ""
) {
  ns <- shiny::NS(id)

  if (length(height) == 1) height <- c(height, 800)

  ifnotchar.int <- function(s) {
    suppressWarnings(
      ifelse(!is.na(as.integer(s)), paste0(as.integer(s), "px"), s)
    )
  }
  height.1 <- ifnotchar.int(height[1])

  header <- shiny::fillRow(
    flex = c(NA, 1, NA),
    class = "tablemodule-header",
    shiny::div(class = "tablemodule-title", title = title, title),
    "",
    DropdownMenu(
      shiny::div(
        class = "tablemodule-info",
        shiny::HTML(paste0("<b>", as.character(title), ".", "</b>", "&nbsp;", as.character(info.text)))
      ),
      width = "250px",
      size = "xs",
      icon = shiny::icon("info"),
      status = "default"
    )
  )

  bslib::card(
    class = "tablemodule",
    full_screen = FALSE,
    style = paste0("height:", height.1, ";overflow:auto;"),
    bslib::as.card_item(div(header)),
    bslib::card_body(
      style = "flex: 1 1 auto; overflow: auto;",
      shiny::fluidRow(
        shiny::column(
          4,
          shiny::selectInput(
            ns("filter_folder"),
            label = "Filter by user",
            choices = NULL,
            selected = NULL,
            width = "100%"
          )
        ),
        shiny::column(
          4,
          shiny::selectInput(
            ns("dest_folder"),
            label = "Destination",
            choices = NULL,
            selected = NULL,
            width = "100%"
          )
        ),
        shiny::column(
          4,
          style = "padding-top: 25px;",
          shiny::actionButton(ns("move_btn"), "Move", icon = shiny::icon("arrows-alt"), class = "btn-sm btn-primary"),
          shiny::actionButton(ns("delete_btn"), "Delete", icon = shiny::icon("trash"), class = "btn-sm btn-danger")
        )
      ),
      DT::DTOutput(ns("files_tbl")) %>% bigLoaders::useSpinner()
    ),
    bslib::as.card_item(
      shiny::div(
        style = "padding: 10px 15px;",
        shiny::span(
          style = "color: #888;",
          shiny::textOutput(ns("status_msg"), inline = TRUE)
        )
      )
    ),
    bslib::card_body(
      class = "card-footer",
      div(class = "caption", shiny::HTML(paste0(
        "<b>", as.character(title), ".</b>",
        "&nbsp;", as.character(caption)
      )))
    )
  )
}

admin_table_datamanager_server <- function(id, auth) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    status <- shiny::reactiveVal("")
    refresh_trigger <- shiny::reactiveVal(0)
    busy <- shiny::reactiveVal(FALSE)

    ## Helper: disable/enable action buttons during operations
    action_btns <- c("move_btn", "delete_btn")
    disable_actions <- function() {
      busy(TRUE)
      for (btn in action_btns) {
        shinyjs::disable(btn)
      }
    }
    enable_actions <- function() {
      for (btn in action_btns) {
        shinyjs::enable(btn)
      }
      busy(FALSE)
    }

    ## -----------------------------------------------------------
    ## Scan all folders for .pgx files
    ## -----------------------------------------------------------
    ## Mapping from display label to internal folder path
    folder_label <- function(path) {
      if (path == "data/") {
        return("(root)")
      }
      if (path == "data_shared/") {
        return("shared")
      }
      if (path == "data_public/") {
        pub_label <- auth$options$PUBLIC_DATASETS_LABEL
        return(if (!is.null(pub_label) && nchar(pub_label) > 0) pub_label else "public")
      }
      ## User directories: "data/user@example.com/" -> "user@example.com"
      lbl <- sub("^data/", "", path)
      lbl <- sub("/$", "", lbl)
      lbl
    }

    folder_list <- shiny::reactive({
      refresh_trigger()
      shiny::req(isTRUE(auth$ADMIN))

      folders <- character(0)

      ## User directories inside PGX.DIR
      user_dirs <- list.dirs(PGX.DIR, full.names = FALSE, recursive = FALSE)
      user_dirs <- grep("@", user_dirs, value = TRUE)
      user_dirs <- sort(user_dirs)
      if (length(user_dirs) > 0) {
        folders <- c(folders, paste0("data/", user_dirs, "/"))
      }

      ## Public dir (root "data/" and "data_shared/" are intentionally omitted
      ## from the dropdowns — files in those locations still appear under the
      ## "(all)" filter via file_data()).
      if (dir.exists(PUBLIC.DIR)) folders <- c(folders, "data_public/")

      folders
    })

    ## Update filter and destination dropdowns
    shiny::observe({
      fl <- folder_list()
      shiny::req(fl)
      labels <- sapply(fl, folder_label)
      ## Preserve current filter selection across refreshes
      current <- input$filter_folder
      sel <- if (!is.null(current) && current %in% c("__all__", fl)) current else "__all__"
      shiny::updateSelectInput(session, "filter_folder",
        choices = c("(all)" = "__all__", stats::setNames(fl, labels)),
        selected = sel
      )
      shiny::updateSelectInput(session, "dest_folder",
        choices = stats::setNames(fl, labels)
      )
    })

    ## -----------------------------------------------------------
    ## Build file table
    ## -----------------------------------------------------------
    file_data <- shiny::reactive({
      refresh_trigger()
      shiny::req(isTRUE(auth$ADMIN))

      opg <- dirname(PGX.DIR) ## project root

      ## Helper: read dataset info from datasets-info.csv in a directory
      read_dataset_info <- function(dir) {
        info_file <- file.path(dir, "datasets-info.csv")
        if (!file.exists(info_file)) {
          return(NULL)
        }
        tryCatch(
          {
            df <- read.csv(info_file, colClasses = "character", stringsAsFactors = FALSE)
            if ("dataset" %in% names(df)) df else NULL
          },
          error = function(e) NULL
        )
      }

      ## Helper: scan a directory for .pgx (non-recursive)
      scan_pgx <- function(dir, label) {
        if (!dir.exists(dir)) {
          return(NULL)
        }
        ff <- list.files(dir, pattern = "\\.pgx$", full.names = TRUE)
        if (length(ff) == 0) {
          return(NULL)
        }
        ## Lookup info from datasets-info.csv
        info_df <- read_dataset_info(dir)
        datasets <- sub("\\.pgx$", "", basename(ff))

        creators <- rep("", length(ff))
        descriptions <- rep("", length(ff))
        if (!is.null(info_df)) {
          idx <- match(datasets, info_df$dataset)
          if ("creator" %in% names(info_df)) {
            creators <- ifelse(!is.na(idx), info_df$creator[idx], "")
          }
          if ("description" %in% names(info_df)) {
            descriptions <- ifelse(!is.na(idx), info_df$description[idx], "")
          }
        }

        ## Fallback: infer creator from user folder name
        folder_email <- sub("^data/", "", label)
        folder_email <- sub("/$", "", folder_email)
        if (grepl("@", folder_email)) {
          creators <- ifelse(nchar(creators) == 0, folder_email, creators)
        }

        data.frame(
          folder = label,
          user = folder_label(label),
          file = datasets,
          description = descriptions,
          creator = unname(creators),
          full_path = ff,
          stringsAsFactors = FALSE
        )
      }

      parts <- list()

      ## Root data/ (non-recursive, top-level pgx only)
      parts[[1]] <- scan_pgx(PGX.DIR, "data/")

      ## User directories
      user_dirs <- list.dirs(PGX.DIR, full.names = FALSE, recursive = FALSE)
      user_dirs <- grep("@", user_dirs, value = TRUE)
      for (ud in user_dirs) {
        parts[[length(parts) + 1]] <- scan_pgx(file.path(PGX.DIR, ud), paste0("data/", ud, "/"))
      }

      ## Shared
      if (dir.exists(SHARE.DIR)) {
        parts[[length(parts) + 1]] <- scan_pgx(SHARE.DIR, "data_shared/")
      }

      ## Public
      if (dir.exists(PUBLIC.DIR)) {
        parts[[length(parts) + 1]] <- scan_pgx(PUBLIC.DIR, "data_public/")
      }

      df <- do.call(rbind, parts)
      if (is.null(df) || nrow(df) == 0) {
        df <- data.frame(
          folder = character(0),
          user = character(0),
          file = character(0),
          description = character(0),
          creator = character(0),
          full_path = character(0),
          stringsAsFactors = FALSE
        )
      }
      df
    })

    ## Filtered data
    display_data <- shiny::reactive({
      df <- file_data()
      shiny::req(df)
      flt <- input$filter_folder
      if (!is.null(flt) && flt != "__all__") {
        df <- df[df$folder == flt, , drop = FALSE]
      }
      df
    })

    ## -----------------------------------------------------------
    ## Render table
    ## -----------------------------------------------------------
    output$files_tbl <- DT::renderDT(
      {
        df <- display_data()
        shiny::req(df)
        show <- df[, c("user", "file", "description", "creator"), drop = FALSE]
        names(show) <- c("User", "File", "Description", "Creator")
        DT::datatable(
          show,
          class = "compact hover",
          rownames = FALSE,
          selection = list(mode = "multiple", target = "row"),
          options = list(
            dom = "tip",
            scrollX = TRUE,
            scrollY = "45vh",
            deferRender = TRUE,
            pageLength = 50,
            order = list(list(0, "asc"), list(1, "asc"))
          )
        ) %>%
          DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "1.3")
      },
      server = FALSE
    )

    ## -----------------------------------------------------------
    ## Helper: resolve destination directory from label
    ## -----------------------------------------------------------
    resolve_dir <- function(label) {
      opg <- dirname(PGX.DIR)
      if (label == "data/") {
        return(PGX.DIR)
      }
      if (label == "data_shared/") {
        return(SHARE.DIR)
      }
      if (label == "data_public/") {
        return(PUBLIC.DIR)
      }
      ## User directory e.g. "data/user@example.com/"
      sub_dir <- sub("^data/", "", label)
      sub_dir <- sub("/$", "", sub_dir)
      file.path(PGX.DIR, sub_dir)
    }

    ## -----------------------------------------------------------
    ## Get selected rows / file paths
    ## -----------------------------------------------------------
    selected_rows <- reactive({
      sel <- input$files_tbl_rows_selected
      if (is.null(sel) || length(sel) == 0) {
        return(NULL)
      }
      df <- display_data()
      df[sel, , drop = FALSE]
    })

    selected_files <- reactive({
      rows <- selected_rows()
      if (is.null(rows)) return(NULL)
      rows$full_path
    })

    ## -----------------------------------------------------------
    ## Move files (with confirmation)
    ## -----------------------------------------------------------
    shiny::observeEvent(input$move_btn, {
      shiny::req(isTRUE(auth$ADMIN))
      files <- selected_files()
      if (is.null(files)) {
        status("No files selected.")
        return()
      }
      dest <- input$dest_folder
      shiny::req(dest)
      dest_label <- folder_label(dest)
      shiny::showModal(shiny::modalDialog(
        title = "Confirm Move",
        shiny::HTML(paste0(
          "<p>Move <b>", length(files), "</b> file(s) to <b>", dest_label, "</b>?</p>",
          "<ul>", paste0("<li>", basename(files), "</li>", collapse = ""), "</ul>"
        )),
        footer = shiny::tagList(
          shiny::actionButton(ns("confirm_move"), "Move", class = "btn-warning"),
          shiny::modalButton("Cancel")
        ),
        easyClose = TRUE
      ))
    })

    shiny::observeEvent(input$confirm_move, {
      shiny::req(isTRUE(auth$ADMIN))
      shiny::req(!busy())
      rows <- selected_rows()
      shiny::req(rows)
      files <- rows$full_path
      dest <- input$dest_folder
      shiny::req(dest)
      dest_dir <- resolve_dir(dest)
      shiny::removeModal()

      disable_actions()
      on.exit(enable_actions())
      shiny::withProgress(message = "Moving files...", value = 0, {
        if (!dir.exists(dest_dir)) {
          dir.create(dest_dir, recursive = TRUE)
        }

        success <- logical(length(files))
        for (i in seq_along(files)) {
          f <- files[i]
          shiny::incProgress(0.5 / length(files), detail = basename(f))
          target <- file.path(dest_dir, basename(f))
          if (normalizePath(f, mustWork = FALSE) == normalizePath(target, mustWork = FALSE)) next
          ok <- file.rename(f, target)
          if (!ok) {
            ok <- file.copy(f, target, overwrite = FALSE)
            if (ok) file.remove(f)
          }
          if (ok) success[i] <- TRUE
        }
        n_ok <- sum(success)

        log_admin_action(
          admin_email = auth$email,
          action = "move",
          subjects = basename(files[success]),
          source_labels = rows$folder[success],
          destination = dest
        )

        ## Update datasets-info.csv in source and destination folders
        source_dirs <- unique(dirname(files))
        dirs_to_update <- c(source_dirs, dest_dir)
        for (i in seq_along(dirs_to_update)) {
          shiny::incProgress(0.5 / length(dirs_to_update), detail = paste("Updating index:", basename(dirs_to_update[i])))
          tryCatch(
            playbase::pgxinfo.updateDatasetFolder(dirs_to_update[i], update.sigdb = FALSE),
            error = function(e) dbg("[admin_datamanager] pgxinfo update error: ", e$message)
          )
        }

        status(paste0("Moved ", n_ok, " of ", length(files), " file(s) to ", dest))
      })
      refresh_trigger(refresh_trigger() + 1)
    })

    ## -----------------------------------------------------------
    ## Delete files (with confirmation)
    ## -----------------------------------------------------------
    shiny::observeEvent(input$delete_btn, {
      shiny::req(isTRUE(auth$ADMIN))
      files <- selected_files()
      if (is.null(files)) {
        status("No files selected.")
        return()
      }
      shiny::showModal(shiny::modalDialog(
        title = "Confirm Delete",
        shiny::HTML(paste0(
          "<p>Are you sure you want to delete <b>", length(files), "</b> file(s)?</p>",
          "<ul>", paste0("<li>", basename(files), "</li>", collapse = ""), "</ul>"
        )),
        footer = shiny::tagList(
          shiny::actionButton(ns("confirm_delete"), "Delete", class = "btn-danger"),
          shiny::modalButton("Cancel")
        ),
        easyClose = TRUE
      ))
    })

    shiny::observeEvent(input$confirm_delete, {
      shiny::req(isTRUE(auth$ADMIN))
      shiny::req(!busy())
      rows <- selected_rows()
      shiny::req(rows)
      files <- rows$full_path
      shiny::removeModal()

      disable_actions()
      on.exit(enable_actions())
      shiny::withProgress(message = "Deleting files...", value = 0, {
        source_dirs <- unique(dirname(files))

        shiny::incProgress(0.3, detail = paste("Removing", length(files), "file(s)"))
        ## Soft-delete: rename ".pgx" -> ".pgx_" so files are hidden from
        ## the app but remain recoverable on disk (same pattern as the
        ## per-user delete in board.loading).
        removed <- vapply(files, function(f) {
          isTRUE(file.rename(f, paste0(f, "_")))
        }, logical(1))
        n_ok <- sum(removed)

        log_admin_action(
          admin_email = auth$email,
          action = "delete",
          subjects = basename(files[removed]),
          source_labels = rows$folder[removed]
        )

        ## Update datasets-info.csv in affected folders
        for (i in seq_along(source_dirs)) {
          shiny::incProgress(0.7 / length(source_dirs), detail = paste("Updating index:", basename(source_dirs[i])))
          tryCatch(
            playbase::pgxinfo.updateDatasetFolder(source_dirs[i], update.sigdb = FALSE),
            error = function(e) dbg("[admin_datamanager] pgxinfo update error: ", e$message)
          )
        }

        status(paste0("Deleted ", n_ok, " of ", length(files), " file(s)."))
      })
      refresh_trigger(refresh_trigger() + 1)
    })

    ## -----------------------------------------------------------
    ## Refresh
    ## -----------------------------------------------------------
    shiny::observeEvent(input$refresh_btn, {
      refresh_trigger(refresh_trigger() + 1)
      status("Refreshed.")
    })

    ## -----------------------------------------------------------
    ## Status message
    ## -----------------------------------------------------------
    output$status_msg <- shiny::renderText({
      status()
    })
  }) ## end of moduleServer
}
