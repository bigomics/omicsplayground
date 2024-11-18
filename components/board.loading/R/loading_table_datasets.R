##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

loading_table_datasets_ui <- function(
    id,
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  ## not sure if this should be here or in settings (IK)
  options <- tagList(
    shiny::checkboxGroupInput(ns("flt_datatype"), "Datatype", choices = "")
  )

  tagList(
    TableModuleUI(
      ns("datasets"),
      info.text = info.text,
      caption = caption,
      width = width,
      height = height,
      title = title,
      options = options
    ),
    div(
      id = "load-action-buttons",
      # this button is needed to trigger download but should be hidden
      shiny::downloadLink(
        ns("download_pgx_btn"),
        label = "",
        icon = NULL,
        width = "0%"
      ),
      # this button is needed to trigger download but should be hidden
      shiny::downloadLink(
        ns("download_zip_btn"),
        label = "",
        icon = NULL,
        width = "0%"
      ),
      shiny::downloadLink(
        ns("recompute_pgx_btn"),
        label = "",
        icon = NULL,
        width = "0%"
      )
    ) ## end of buttons div
  )
}

loading_table_datasets_server <- function(id,
                                          pgx_topdir,
                                          pgx_shared_dir,
                                          pgx_public_dir,
                                          auth,
                                          loadAndActivatePGX,
                                          loadPGX,
                                          refresh_shared,
                                          reload_pgxdir_public,
                                          reload_pgxdir,
                                          recompute_pgx,
                                          loadbutton,
                                          newuploadbutton,
                                          new_upload) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    share_pgx <- reactiveVal(NULL)

    getPGXINFO <- shiny::reactive({
      shiny::req(auth$logged)
      if (is.null(auth$logged) || !auth$logged) {
        return(NULL)
      }

      ## upstream trigger
      reload_pgxdir()
      pgxdir <- auth$user_dir

      shiny::withProgress(message = "Checking datasets library...", value = 0.33, {
        need_update <- playbase::pgxinfo.needUpdate(pgxdir,
          check.sigdb = FALSE,
          verbose = FALSE
        )
      })

      ## before reading the info file, we need to update for new files
      if (need_update) {
        pgx.showSmallModal("Updating your library<br>Please wait...")
        shiny::withProgress(message = "Updating your library...", value = 0.33, {
          playbase::pgxinfo.updateDatasetFolder(
            pgxdir,
            force = FALSE,
            delete.old = TRUE,
            update.sigdb = FALSE
          )
        })
        shiny::removeModal(session)
      }

      info <- playbase::pgxinfo.read(pgxdir, file = "datasets-info.csv")
      shiny::removeModal(session)
      return(info)
    })

    getFilteredPGXINFO <- shiny::reactive({
      ## get the filtered table of pgx datasets
      df <- getPGXINFO()
      if (is.null(df)) {
        return(NULL)
      }

      pgxdir <- auth$user_dir
      pgxfiles <- dir(pgxdir, pattern = ".pgx$")
      sel <- sub("[.]pgx$", "", df$dataset) %in% sub("[.]pgx$", "", pgxfiles)
      df <- df[sel, , drop = FALSE]

      ## Apply filters
      if (nrow(df) > 0) {
        f1 <- f2 <- f3 <- rep(TRUE, nrow(df))
        notnull <- function(x) !is.null(x) && length(x) > 0 && x[1] != "" && !is.na(x[1])
        if (notnull(input$flt_datatype)) f2 <- (df$datatype %in% input$flt_datatype)
        df <- df[which(f1 & f2 & f3), , drop = FALSE]
        df$date <- as.Date(df$date, format = "%Y-%m-%d")
        df <- df[order(df$date, decreasing = TRUE), ]
        if (nrow(df) > 0) rownames(df) <- nrow(df):1
      }

      kk <- unique(c(
        "dataset", "description", "organism", "datatype", "nsamples",
        "ngenes", "nsets", "conditions", "date",
        "creator"
      ))
      kk <- intersect(kk, colnames(df))
      df <- df[, kk, drop = FALSE]
      df
    })

    observeEvent(getPGXINFO(), {
      df <- getPGXINFO()
      if (is.null(df)) {
        return()
      }
      datatypes <- sort(setdiff(df$datatype, c(NA, "")))
      shiny::updateCheckboxGroupInput(session, "flt_datatype", choices = datatypes)
    })


    ## -------------------------------------------------------------------
    ## make a pgx public (i.e. share publicly)
    ## -------------------------------------------------------------------
    observeEvent(
      input$share_public_pgx,
      {
        selected_row <- as.numeric(stringr::str_split(input$share_public_pgx, "_row_")[[1]][2])
        pgx_name <- table_data()[selected_row, "dataset"]

        alert_val <- shinyalert::shinyalert(
          inputId = "share_public_confirm",
          title = paste("Share this dataset?"),
          paste(
            "Your dataset", pgx_name, "will be copied",
            "to the public folder. Other users will be able import and explore",
            "this dataset. Are you sure?"
          ),
          html = TRUE,
          showCancelButton = TRUE,
          showConfirmButton = TRUE
        )
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    )


    observeEvent(input$share_public_confirm, {
      if (input$share_public_confirm) {
        selected_row <- as.numeric(stringr::str_split(input$share_public_pgx, "_row_")[[1]][2])
        pgx_name <- table_data()[selected_row, "dataset"]
        pgx_name <- sub("[.]pgx$", "", pgx_name)
        pgx_file <- file.path(auth$user_dir, paste0(pgx_name, ".pgx"))
        new_pgx_file <- file.path(pgx_public_dir, paste0(pgx_name, ".pgx"))

        ## abort if file exists
        if (file.exists(new_pgx_file)) {
          shinyalert::shinyalert(
            title = "Oops! File exists...",
            paste(
              "There is already a dataset called", pgx_name,
              "in the Public folder. Sorry about that! Please rename your file",
              "if you still want to share it."
            )
          )
          return()
        }

        ## file.copy(from = pgx_file, to = new_pgx_file)
        shiny::withProgress(message = "Copying file to public folder...", value = 0.33, {
          pgx0 <- playbase::pgx.load(pgx_file)
          unknown.creator <- pgx0$creator %in% c(NA, "", "user", "anonymous", "unknown")
          if ("creator" %in% names(pgx0) && !unknown.creator) {
            file.copy(from = pgx_file, to = new_pgx_file)
          } else {
            pgx0$creator <- auth$email
            if (pgx0$creator %in% c(NA, "", "user", "anonymous", "unknown")) pgx0$creator <- "unknown"
            playbase::pgx.save(pgx0, file = new_pgx_file)
          }
        })

        reload_pgxdir_public(reload_pgxdir_public() + 1)

        shinyalert::shinyalert(
          title = "Successfully shared!",
          paste(
            "Your dataset", pgx_name, "has now been successfully",
            "been publicly shared. Thank you!"
          )
        )
      }
    })

    shiny::observeEvent(loadbutton(), {
      pgxfile <- table_selected_pgx()
      # Make sure there is a row selected
      if (is.null(pgxfile)) {
        return(NULL)
      }
      pgxfilename <- file.path(auth$user_dir, pgxfile)
      if (!file.exists(pgxfilename)) {
        message("[LoadingBoard@load_react] ERROR pgxfile not found : ", pgxfilename, "\n")
        return(NULL)
      }

      loadAndActivatePGX(pgxfile)
    })

    ## ---------------------------- create table module -----------------------------------

    table_data <- shiny::reactive({
      df <- getFilteredPGXINFO()
      if (is.null(df)) {
        return(NULL)
      }
      df$dataset <- sub("[.]pgx$", "", df$dataset)
      df$conditions <- gsub("[,]", " ", df$conditions)
      df$conditions <- sapply(as.character(df$conditions), andothers, split = " ", n = 5)
      df$description <- gsub("[_]", " ", df$description) ## replace underscore
      ##  df$description <- playbase::shortstring(as.character(df$description), 200)
      df$nsets <- NULL
      return(df)
    })

    pgxTable_DT <- reactive({
      df <- table_data()
      is.dt <- is.data.frame(df)
      if (is.null(df) || !is.dt || nrow(df) == 0) {
        shinyalert::shinyalert(
          title = "Empty?",
          text = paste(
            "Your dataset library seems empty. Please upload new data or import",
            "a dataset from the public datasets folder."
          )
        )
      }
      validate(need(nrow(df) > 0, "Need at least one dataset!"))

      ## df$creator <- NULL
      target1 <- grep("date", colnames(df))
      target2 <- grep("description", colnames(df))
      target3 <- grep("conditions", colnames(df))
      target4 <- grep("dataset", colnames(df))

      # create action menu for each row
      menus <- c()
      for (i in 1:nrow(df)) {
        download_pgx_menuitem <- NULL
        share_public_menuitem <- NULL
        share_dataset_menuitem <- NULL
        delete_pgx_menuitem <- NULL

        if (auth$options$ENABLE_PGX_DOWNLOAD) {
          download_pgx_menuitem <- shiny::actionButton(
            ns(paste0("download_pgx_row_", i)),
            label = "Download PGX",
            icon = shiny::icon("download"),
            class = "btn btn-outline-dark",
            style = "border: none;",
            width = "100%",
            onclick = paste0('Shiny.onInputChange(\"', ns("download_pgx"), '\",this.id,{priority: "event"})')
          )
        }
        if (auth$options$ENABLE_PUBLIC_SHARE && dir.exists(pgx_public_dir)) {
          share_public_menuitem <- shiny::actionButton(
            ns(paste0("share_public_row_", i)),
            label = "Share public",
            icon = shiny::icon("share-nodes"),
            class = "btn btn-outline-info",
            style = "border: none;",
            width = "100%",
            onclick = paste0('Shiny.onInputChange(\"', ns("share_public_pgx"), '\",this.id,{priority: "event"})')
          )
        }
        if (auth$options$ENABLE_USER_SHARE && dir.exists(pgx_shared_dir)) {
          share_dataset_menuitem <- shiny::actionButton(
            ns(paste0("share_dataset_row_", i)),
            label = "Share with user",
            icon = shiny::icon("share-nodes"),
            class = "btn btn-outline-info",
            style = "border: none;",
            width = "100%",
            onclick = paste0('Shiny.onInputChange(\"', ns("share_pgx"), '\",this.id,{priority: "event"})')
          )
        }

        ## instead of disabling we grey out and have a popup message
        delete_pgx_menuitem <- shiny::actionButton(
          ns(paste0("delete_dataset_row_", i)),
          label = "Delete dataset",
          icon = shiny::icon("trash"),
          class = "btn btn-outline-danger",
          style = "border: none;",
          width = "100%",
          onclick = paste0('Shiny.onInputChange(\"', ns("delete_pgx"), '\",this.id,{priority: "event"});')
        )
        recompute_pgx_menuitem <- shiny::actionButton(
          ns(paste0("recompute_pgx_row_", i)),
          label = "Reanalyse",
          icon = shiny::icon("gears"),
          class = "btn btn-outline-dark",
          style = "border: none;",
          width = "100%",
          onclick = paste0('Shiny.onInputChange(\"', ns("recompute_pgx"), '\",this.id,{priority: "event"});')
        )

        new_menu <- DropdownMenu(
          div(
            style = "width: 160px;",
            div(
              download_pgx_menuitem,
              shiny::actionButton(
                ns(paste0("download_zip_row_", i)),
                label = "Download ZIP",
                icon = shiny::icon("file-archive"),
                class = "btn btn-outline-dark",
                style = "border: none;",
                width = "100%",
                onclick = paste0('Shiny.onInputChange(\"', ns("download_zip"), '\",this.id,{priority: "event"})')
              ),
              recompute_pgx_menuitem,
              share_public_menuitem,
              share_dataset_menuitem,
              delete_pgx_menuitem
            )
          ),
          size = "sm",
          icon = shiny::icon("ellipsis-vertical"),
          status = "outline-dark",
          circle = FALSE,
          border = 0
        )
        menus <- c(menus, as.character(new_menu))
      }


      observeEvent(input$share_pgx,
        {
          share_pgx(input$share_pgx)
        },
        ignoreInit = TRUE
      )


      DT::datatable(
        df,
        class = "compact hover",
        rownames = menus,
        escape = FALSE,
        editable = list(
          target = "cell",
          disable = list(columns = c(1, 3:ncol(df)))
        ),
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
            ),
            list(sortable = FALSE, targets = ncol(df))
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
        DT::formatStyle(0, target = "row", fontSize = "20px", lineHeight = "95%")
    }

    table_module <- TableModuleServer(
      "datasets",
      func = pgxTable.RENDER,
      func2 = pgxTable_modal.RENDER,
      selector = "single"
    )

    ## --------------- edit rows of  pgxtable ---------------------
    observeEvent(
      input[["datasets-datatable_cell_edit"]],
      {
        row <- input[["datasets-datatable_cell_edit"]]$row
        col <- input[["datasets-datatable_cell_edit"]]$col
        val <- input[["datasets-datatable_cell_edit"]]$value

        df <- table_data()
        col_edited <- colnames(df)[col]

        dataset_edited <- df$dataset[row]
        pgxinfo <- getPGXINFO()

        row_edited <- match(dataset_edited, pgxinfo$dataset)
        pgxinfo[row_edited, col_edited] <- val
        fname <- file.path(auth$user_dir, "datasets-info.csv")
        write.csv(pgxinfo, fname)

        ## also rewrite description in actual pgx file
        pgx_name <- dataset_edited
        pgx_file <- file.path(auth$user_dir, paste0(pgx_name, ".pgx"))
        pgx <- playbase::pgx.load(pgx_file, verbose = FALSE) ## override any name

        row_edited <- match(dataset_edited, pgxinfo$dataset)
        new_val <- pgxinfo[row_edited, col_edited]
        pgx[[col_edited]] <- new_val

        dbg("[datasets-datatable_cell_edit] updating", col_edited, " -> ", new_val)
        dbg("[datasets-datatable_cell_edit] saving changes to", pgx_file)
        playbase::pgx.save(pgx, file = pgx_file)
        remove(pgx)
        dbg("[datasets-datatable_cell_edit] done!")
      },
      ignoreInit = TRUE
    )


    table_selected_pgx <- shiny::reactive({
      req(table_module)

      sel <- table_module$rows_selected()

      if (is.null(sel) || length(sel) == 0) {
        return(NULL)
      }
      df <- table_data()
      pgxfile <- as.character(df$dataset[sel])
      pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx") ## add/replace .pgx
      pgxfile
    })

    selected_sharePGX <- reactive({
      selected_row <- as.numeric(stringr::str_split(input$share_pgx, "_row_")[[1]][2])
      pgx_name <- table_data()[selected_row, "dataset"]
      pgx_name
    })

    ## ----------------- Download ZIP -----------------------
    observeEvent(input$download_zip,
      {
        shinyjs::click(id = "download_zip_btn")
      },
      ignoreNULL = TRUE
    )

    output$download_zip_btn <- shiny::downloadHandler(
      filename = function() {
        sel <- row_idx <- as.numeric(stringr::str_split(input$download_zip, "_row_")[[1]][2])
        df <- getFilteredPGXINFO()
        pgxfile <- as.character(df$dataset[sel])
        pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx") ## add/replace .pgx
        newfile <- sub("pgx$", "zip", pgxfile)
        newfile
      },
      content = function(file) {
        sel <- row_idx <- as.numeric(stringr::str_split(input$download_zip, "_row_")[[1]][2])
        df <- getFilteredPGXINFO()
        pgxfile <- as.character(df$dataset[sel])
        pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx") ## add/replace .pgx
        pgxname <- sub("[.]pgx$", "", pgxfile)

        pgx <- loadPGX(pgxfile)
        dir.create(tmp <- tempfile())
        tmp2 <- file.path(tmp, pgxname)
        dir.create(tmp2)

        exp.matrix <- sign(pgx$model.parameters$exp.matrix)
        exp.matrix <- playbase::contrastAsLabels(exp.matrix) ## new recommended style
        exp.matrix[is.na(exp.matrix)] <- ""

        write.csv(round(pgx$counts, digits = 2), file = file.path(tmp2, "counts.csv"))
        write.csv(pgx$samples, file = file.path(tmp2, "samples.csv"))
        write.csv(exp.matrix, file = file.path(tmp2, "contrasts.csv"))
        write.csv(round(pgx$X, digits = 4), file = file.path(tmp2, "normalized.csv"))

        zipfile <- tempfile(fileext = ".zip")
        zip::zip(zipfile,
          files = paste0(pgxname, "/", c(
            "counts.csv", "samples.csv",
            "contrasts.csv", "normalized.csv"
          )),
          root = tmp
        )
        file.copy(zipfile, file)
        remove(pgx)
        gc()
      }
    )
    ## ---------------- RECOMPUTE PGX -------------------
    shiny::observeEvent(input$recompute_pgx, {
      shinyalert::shinyalert(
        title = "Reanalyze",
        text = "Are you sure you want to reanalyze your data? All current contrasts will be kept.",
        showCancelButton = TRUE,
        cancelButtonText = "Cancel",
        confirmButtonText = "OK",
        callbackR = function(x) {
          if (x) {
            # Load PGX
            sel <- row_idx <- as.numeric(stringr::str_split(input$recompute_pgx, "_row_")[[1]][2])
            df <- getFilteredPGXINFO()
            pgxfile <- as.character(df$dataset[sel])
            pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx")
            pgx <- loadPGX(pgxfile)
            load_uploaded_data <- shiny::reactiveVal(NULL)
            reload_pgxdir <- shiny::reactiveVal(0)
            recompute_pgx(pgx)

            # bigdash.selectTab(session, "upload-tab")
            # shinyjs::runjs('$("[data-value=\'Upload\']").click();') # Should be Comparisons?

            return(0)
          } else {
            return(0)
          }
        }
      )
    })

    ## ---------------- DOWNLOAD PGX FILE ----------------
    observeEvent(input$download_pgx,
      {
        shinyjs::click(id = "download_pgx_btn")
      },
      ignoreNULL = TRUE
    )

    output$download_pgx_btn <- shiny::downloadHandler(

      ## filename = "userdata.pgx",
      filename = function() {
        sel <- row_idx <- as.numeric(stringr::str_split(input$download_pgx, "_row_")[[1]][2])
        df <- getFilteredPGXINFO()
        pgxfile <- as.character(df$dataset[sel])
        pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx")
        pgxfile
      },
      content = function(file) {
        sel <- row_idx <- as.numeric(stringr::str_split(input$download_pgx, "_row_")[[1]][2])
        df <- getFilteredPGXINFO()
        pgxfile <- as.character(df$dataset[sel])
        pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx")
        pgx <- loadPGX(pgxfile)
        temp <- tempfile()
        save(pgx, file = temp)
        file.copy(temp, file)
      }
    )

    ## --------------- DELETE PGX ------------------

    shiny::observeEvent(input$delete_pgx,
      {
        row_idx <- as.numeric(stringr::str_split(input$delete_pgx, "_row_")[[1]][2])
        df <- getFilteredPGXINFO()
        pgxfile <- as.character(df$dataset[row_idx])
        pgxname <- sub("[.]pgx$", "", pgxfile)
        pgxfile <- paste0(pgxname, ".pgx") ## add/replace .pgx
        pgxfile1 <- file.path(auth$user_dir, pgxfile)
        sel <- NULL

        deletePGX <- function(x) {
          if (input$confirmdelete) {
            pgxfile2 <- paste0(pgxfile1, "_") ## mark as deleted
            file.rename(pgxfile1, pgxfile2)

            ## also delete entry in PGXINFO and allFC (bit slôw)
            pgx.dir <- auth$user_dir
            playbase::pgxinfo.delete(pgx.dir, pgxname)

            ## signal upstream for update
            reload_pgxdir(reload_pgxdir() + 1)
          }
        }

        if (auth$options$ENABLE_DELETE) {
          shinyalert::shinyalert(
            "Delete this dataset?",
            paste("Are you sure you want\nto delete '", pgxfile, "'?"),
            confirmButtonText = "Delete",
            showCancelButton = TRUE,
            callbackR = deletePGX,
            inputId = "confirmdelete"
          )
        } else {
          msg <- paste(
            "Deleting is disabled.",
            "Please <a href='https://events.bigomics.ch/upgrade' target='_blank'>",
            "<b><u>upgrade</u></b></a> your account to enable it."
          )
          shinyalert::shinyalert(
            title = "Oops!",
            text = HTML(msg),
            showCancelButton = TRUE,
            showConfirmButton = FALSE,
            html = TRUE
          )
        }
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    )

    ## -------------------------------------------------------------------
    ## share dataset with specific user
    ## -------------------------------------------------------------------

    #' will either take free-form or selected email
    input_share_user <- reactive({
      share_user <- input$share_user
      if (!is.null(input$share_user2) && input$share_user2 != "") {
        share_user <- input$share_user2
      }
      share_user <- trimws(tolower(share_user)) ## force lower case email
      share_user
    })

    observeEvent(input$share_user2, {
      if (input$share_user2 != "") {
        updateTextInput(session, "share_user", value = "")
      }
    })

    share_dialog <- function(pgxname, choices = NULL) {
      select_user <- NULL
      if (!is.null(choices) && length(choices) > 0) {
        select_user <- tagList(
          shiny::selectizeInput(ns("share_user2"),
            "Or select a coworker from the list:",
            choices = c("Choose coworker..." = "", choices),
            options = list(create = TRUE)
          )
        )
      }

      shiny::modalDialog(
        tagList(
          paste("Your dataset", pgxname, "will be shared with the user below."),
          br(), br(),
          shiny::textInput(ns("share_user"),
            label = "Enter email who will receive the dataset:",
            placeholder = "Type email..."
          ),
          select_user,

          ## shiny::textInput(ns("share_user2"), "Re-enter use email:")

          shiny::textOutput(ns("error_alert")) %>%
            tagAppendAttributes(style = "color: red;")
        ),
        title = "Share this dataset?",
        footer = tagList(
          actionButton(ns("share_dialog_cancel"), "Cancel"),
          actionButton(ns("share_dialog_confirm"), "Share")
        )
      )
    }


    # put user dataset into shared folder
    observeEvent(input$share_pgx,
      {
        ## sharing folder has to exists
        if (!dir.exists(pgx_shared_dir)) {
          shinyalert::shinyalert(
            title = "Oops! Cannot share...",
            text = paste(
              "This server does not support sharing.",
              "Please contact your administrator."
            )
          )
          share_pgx(NULL)
          return()
        }

        ## user has to be logged in and have email for them to share with other users
        if (auth$email == "") {
          shinyalert::shinyalert(
            title = "Oops! Cannot share...",
            text = paste(
              "You need to be logged in with a valid email",
              "address to share pgx files with other users."
            )
          )
          share_pgx(NULL)
          return()
        }

        ## check how many are already in queue
        pp <- paste0("__from__", auth$email, "__$")
        num_shared_queue <- length(dir(pgx_shared_dir, pattern = pp))

        if (num_shared_queue >= opt$MAX_SHARED_QUEUE) { ## NB opt is global...
          shinyalert::shinyalert(
            title = "Oops! Too many shared...",
            text = paste(
              "You have already too many shared datasets in the waiting queue.",
              "Please contact your administrator."
            )
          )
          share_pgx(NULL)
          return()
        }

        ## show share dialog
        coworkers <- get_coworkers(pgx_topdir, auth$email)
        pgxname <- sub("[.]pgx$", "", selected_sharePGX())
        if (length(coworkers) == 0) coworkers <- NULL
        shiny::showModal(share_dialog(pgxname, choices = coworkers))
      },
      ignoreNULL = TRUE
    )

    observeEvent(input$share_dialog_cancel, {
      share_pgx(NULL)
      output$error_alert <- renderText({
        ""
      })
      shiny::removeModal()
    })

    observeEvent(input$share_dialog_confirm, {
      share_user <- input_share_user()
      if (share_user == "") {
        output$error_alert <- renderText({
          "Please enter an email."
        })
        return()
      }

      if (!is_valid_email(share_user)) {
        output$error_alert <- renderText({
          "Email is not valid. Please use only work or business emails."
        })
        return()
      }
      if (share_user == auth$email) {
        output$error_alert <- renderText({
          "Error. You cannot share with yourself..."
        })
        return()
      }

      output$error_alert <- renderText({
        ""
      })

      shiny::removeModal()


      pgx_name <- selected_sharePGX()

      alert_val <- shinyalert::shinyalert(
        inputId = "share_confirm",
        title = "Are you sure?",
        tagList(
          paste("Your dataset", pgx_name, "will be shared with", share_user)
        ),
        html = TRUE,
        showCancelButton = TRUE,
        showConfirmButton = TRUE
      )
    })


    observeEvent(input$share_confirm, {
      # if confirmed, then share the data
      if (input$share_confirm) {
        pgx_name <- selected_sharePGX()
        pgx_name <- sub("[.]pgx$", "", pgx_name)
        pgx_file <- file.path(auth$user_dir, paste0(pgx_name, ".pgx"))

        ## The shared file will be copied to the data_shared
        ## folder with the name of the sender and receiver in the
        ## file name.
        share_user <- input_share_user()
        current_user <- auth$email
        new_pgx_file <- file.path(
          pgx_shared_dir,
          paste0(
            pgx_name, ".pgx", "__to__", share_user,
            "__from__", current_user, "__"
          )
        )

        if (file.exists(new_pgx_file)) {
          shinyalert::shinyalert(
            title = "Oops! File exists...",
            paste(
              "There is already a dataset called", pgx_name,
              "being shared. Please rename your file."
            )
          )
          return()
        }

        ## load and save the pgx file to new directory
        shiny::withProgress(message = "Preparing to share...", value = 0.33, {
          pgx0 <- playbase::pgx.load(pgx_file)
          unknown.creator <- pgx0$creator %in% c(NA, "", "user", "anonymous", "unknown")
          if ("creator" %in% names(pgx0) && !unknown.creator) {
            file.copy(from = pgx_file, to = new_pgx_file)
          } else {
            pgx0$creator <- auth$email
            if (pgx0$creator %in% c(NA, "", "user", "anonymous", "unknown")) pgx0$creator <- "unknown"
            playbase::pgx.save(pgx0, file = new_pgx_file)
          }

          ## write transaction to log file
          log.entry <- data.frame(date = date(), from = "jane@demo.com", to = "tarzan@demo.com", file = "example-data.pgx")
          log.file <- file.path(pgx_shared_dir, "PGXSHARE-TRANSACTIONS.log")
          log.entry <- data.frame(date = date(), from = auth$email, to = share_user, file = paste0(pgx_name, ".pgx"))
          if (file.exists(log.file)) {
            write.table(log.entry, file = log.file, col.names = FALSE, row.names = FALSE, sep = ",", append = TRUE)
          } else {
            write.table(log.entry, file = log.file, col.names = TRUE, row.names = FALSE, sep = ",")
          }
        })

        share_user <- input_share_user()
        shinyalert::shinyalert(
          title = "Successfully shared!",
          paste(
            "Your dataset", pgx_name, "has now been successfully",
            "been shared with", share_user
          )
        )

        # send email to user
        sender <- trimws(auth$email)
        gmail_creds <- file.path(ETC, "gmail_creds")
        sendShareMessage(pgx_name, sender, share_user, path_to_creds = gmail_creds)

        refresh_shared(refresh_shared() + 1)
      }

      share_pgx(NULL)
    })




    ## please refer to TableModule for return values

    return(table_module)
  })
}
