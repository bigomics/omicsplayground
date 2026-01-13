##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

loading_table_datasets_public_ui <- function(
  id,
  title,
  info.text,
  caption,
  height,
  width,
  delete_button = FALSE,
  load_button = FALSE
) {
  ns <- shiny::NS(id)

  tagList(
    TableModuleUI(
      ns("datasets"),
      info.text = info.text,
      width = width,
      caption = caption,
      height = height,
      title = title
    ),
    div(
      if (load_button) {
        shiny::actionButton(
          ns("loadbutton"),
          label = "Load selected",
          icon = icon("file-import"),
          class = "btn btn-primary"
        )
      },
      shiny::actionButton(
        ns("importbutton"),
        label = "Import dataset",
        icon = icon("file-import"),
        class = "btn btn-primary"
      ),
      if (delete_button) {
        shiny::actionButton(
          ns("deletebutton"),
          label = "Delete dataset",
          icon = icon("trash"),
          class = "btn btn-danger"
        )
      }
    )
  )
}

loading_table_datasets_public_server <- function(id,
                                                 pgx_public_dir,
                                                 reload_pgxdir_public,
                                                 auth,
                                                 reload_pgxdir,
                                                 loadAndActivatePGX = NULL) {
  moduleServer(id, function(input, output, session) {
    getPGXINFO_PUBLIC <- shiny::reactive({
      shiny::req(auth$logged)
      if (is.null(auth$logged) || !auth$logged) {
        warning("[LoadingBoard:getPGXINFO_PUBLIC] user not logged in!")
        return(NULL)
      }

      reload_pgxdir_public()

      ## update meta files
      info <- NULL
      shiny::withProgress(message = "Checking datasets library...", value = 0.33, {
        need_update <- playbase::pgxinfo.needUpdate(pgx_public_dir,
          check.sigdb = FALSE, verbose = FALSE
        )
      })

      if (need_update) {
        dbg("[loading_server:getPGXINFO_PUBLIC] updating public pgxdir =", pgx_public_dir)
        pgx.showSmallModal("Updating datasets library<br>Please wait...")
        shiny::withProgress(message = "Updating datasets library...", value = 0.33, {
          ## before reading the info file, we need to update for new files
          playbase::pgxinfo.updateDatasetFolder(pgx_public_dir, update.sigdb = FALSE)
        })
        shiny::removeModal(session)
      }

      info <- playbase::pgxinfo.read(pgx_public_dir, file = "datasets-info.csv")
      return(info)
    })

    getFilteredPGXINFO_PUBLIC <- shiny::reactive({
      ## get the filtered table of pgx datasets
      req(auth$logged)
      if (!auth$logged) {
        return(NULL)
      }

      df <- getPGXINFO_PUBLIC()
      shiny::req(df)

      pgxfiles <- dir(pgx_public_dir, pattern = ".pgx$")
      sel <- sub("[.]pgx$", "", df$dataset) %in% sub("[.]pgx$", "", pgxfiles)
      df <- df[sel, , drop = FALSE]

      ## Apply filters
      if (nrow(df) > 0) {
        f1 <- f2 <- f3 <- rep(TRUE, nrow(df))
        notnull <- function(x) !is.null(x) && length(x) > 0 && x[1] != "" && !is.na(x[1])
        if (notnull(input$flt_datatype)) f2 <- (df$datatype %in% input$flt_datatype)
        if (notnull(input$flt_organism)) f3 <- (df$organism %in% input$flt_organism)
        df <- df[which(f1 & f2 & f3), , drop = FALSE]
        df$date <- as.Date(df$date, format = "%Y-%m-%d")
        df <- df[order(df$date, decreasing = TRUE), ]
        if (nrow(df) > 0) rownames(df) <- nrow(df):1
      }

      kk <- unique(c(
        "dataset", "description", "datatype", "nsamples",
        "ngenes", "nsets", "conditions", "date", "organism",
        "creator"
      ))
      kk <- intersect(kk, colnames(df))
      df <- df[, kk, drop = FALSE]

      df$dataset <- sub("[.]pgx$", "", df$dataset)
      df$conditions <- gsub("[,]", " ", df$conditions)
      df$conditions <- sapply(as.character(df$conditions), andothers, split = " ", n = 5)
      df$description <- gsub("[_]", " ", df$description) ## replace underscore
      ##      df$description <- playbase::shortstring(as.character(df$description), 200)
      df$nsets <- NULL
      df$organism <- NULL

      df
    })

    # disable buttons when no row is selected; enable when one is selected
    observeEvent(pgxtable_public$rows_selected(),
      {
        if (is.null(pgxtable_public$rows_selected())) {
          shinyjs::disable(id = "loadbutton")
          shinyjs::disable(id = "importbutton")
          shinyjs::disable(id = "deletebutton")
        } else {
          shinyjs::enable(id = "loadbutton")
          shinyjs::enable(id = "importbutton")
          shinyjs::enable(id = "deletebutton")
        }
      },
      ignoreNULL = FALSE
    )

    observeEvent(input$loadbutton, {
      ## Load button: loads the dataset directly from public directory without importing
      if (is.null(loadAndActivatePGX)) {
        warning("[loading_table_datasets_public] loadAndActivatePGX function not provided")
        return(NULL)
      }

      selected_row <- pgxtable_public$rows_selected()
      pgx_name <- pgxtable_public$data()[selected_row, "dataset"]

      ## Load directly from the public directory without copying
      loadAndActivatePGX(pgx_name, pgxdir = pgx_public_dir)
    })

    observeEvent(input$importbutton, {
      selected_row <- pgxtable_public$rows_selected()
      pgx_name <- pgxtable_public$data()[selected_row, "dataset"]
      pgx_file <- file.path(pgx_public_dir, paste0(pgx_name, ".pgx"))
      pgx_path <- auth$user_dir
      new_pgx_file <- file.path(pgx_path, paste0(pgx_name, ".pgx"))

      ## check number of datasets. If deletion is disabled, we count also .pgx_ files... :)
      numpgx <- length(dir(pgx_path, pattern = "*.pgx$"))
      if (!auth$options$ENABLE_DELETE) numpgx <- length(dir(pgx_path, pattern = "*.pgx$|*.pgx_$"))
      maxpgx <- as.integer(auth$options$MAX_DATASETS)
      if (numpgx >= maxpgx) {
        shinyalert_storage_full(numpgx, maxpgx) ## ui-alerts.R
        return(NULL)
      }

      if (file.exists(new_pgx_file)) {
        shinyalert::shinyalert(
          title = "Oops! File exists...",
          paste(
            "There is already a dataset called", pgx_name,
            "in your dataset folder. Please delete your file first."
          )
        )
        return()
      }

      ## Copy the file from Public folder to user folder
      shiny::withProgress(message = "Importing dataset...", value = 0.33, {
        file.copy(from = pgx_file, to = new_pgx_file)
        playbase::pgxinfo.updateDatasetFolder(pgx_path, update.sigdb = FALSE)
      })

      shinyalert::shinyalert(
        "Dataset imported",
        paste(
          "The dataset", pgx_name, "from", tolower(auth$options$PUBLIC_DATASETS_LABEL), "has now been successfully imported",
          "to your library. Feel free to load it as usual!"
        )
      )

      reload_pgxdir(reload_pgxdir() + 1)
    })

    observeEvent(input$deletebutton, {
      selected_row <- pgxtable_public$rows_selected()
      pgx_name <- pgxtable_public$data()[selected_row, "dataset"]
      pgx_file <- file.path(pgx_public_dir, paste0(pgx_name, ".pgx"))
      file.rename(pgx_file, paste0(pgx_file, "_"))
      reload_pgxdir_public(reload_pgxdir_public() + 1)
    })

    pgxTable_DT <- reactive({
      df <- getFilteredPGXINFO_PUBLIC()
      ## shiny::req(df)

      # need this, otherwise there is an error on user logout
      validate(need(nrow(df) > 0, paste("No", tolower(auth$options$PUBLIC_DATASETS_LABEL), "!")))

      df$creator <- sub("@.*", "", df$creator) ## hide full email

      target1 <- grep("date", colnames(df))
      target2 <- grep("description", colnames(df))
      target3 <- grep("conditions", colnames(df))
      target4 <- grep("dataset", colnames(df))

      DT::datatable(
        df,
        class = "compact hover",
        rownames = FALSE,
        editable = FALSE,
        extensions = c("Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "ft",
          pageLength = 9999,
          scrollX = FALSE,
          scrollY = "55vh",
          scrollResize = TRUE,
          deferRender = TRUE,
          autoWidth = TRUE,
          columnDefs = list(
            list(width = "60px", targets = target1),
            ##            list(width = "30vw", targets = target2),
            list(
              targets = target2, ## with no rownames column 1 is column 2
              render = DT::JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 150 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 150) + '...</span>' : data;",
                "}"
              )
            )
          )
        ) ## end of options.list
      )
    })

    pgxTable.RENDER <- function() {
      pgxTable_DT() %>%
        DT::formatStyle(0, target = "row", fontSize = "12px", lineHeight = "95%")
    }

    pgxTable_modal.RENDER <- function() {
      pgxTable_DT() %>%
        DT::formatStyle(0, target = "row", fontSize = "16px", lineHeight = "95%")
    }


    pgxtable_public <- TableModuleServer(
      "datasets",
      func = pgxTable.RENDER,
      func2 = pgxTable_modal.RENDER,
      selector = "single"
    )

    ## please refer to TableModule for return values

    return(pgxtable_public)
  })
}
