##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

LoadingBoard <- function(id,
                         pgx_dir,
                         pgx,
                         auth,
                         limits = c(
                           "samples" = 1000, "comparisons" = 20,
                           "genes" = 20000, "genesets" = 10000,
                           "datasets" = 10
                         ),
                         enable_upload = TRUE,
                         enable_delete = TRUE,
                         enable_save = TRUE,
                         enable_userdir = TRUE,
                         r_global) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    loadedDataset <- shiny::reactiveVal(0) ## counts/trigger dataset upload

    ## info that's needed
    r_local <- reactiveValues(
      selected_row = NULL,
      found_example_trigger = NULL
    )

    observeEvent(pgxtable$rows_selected(), {
      r_local$selected_row <- pgxtable$rows_selected()
    })

    observeEvent(r_global$load_example_trigger, {

      data_names <- as.character(pgxtable$data()$dataset)

      # get the row which corresponds to "example-data"
      example_row <- which(data_names == "example-data")[1]

      # if not found, throw error modal that example-data doesnt exist
      if (is.na(example_row)) {
        shiny::showModal(modalDialog(
          title = "No example data found",
          size = "l",
          'Sorry, it appears that the example dataset cannot be found. You may
          have deleted it in a previous session. You can still load a copy of the
          example dataset by clicking the "Upload New Data" button.'
        ))
        r_global$load_example_trigger <- NULL
        return(NULL)
      } else {
        shinyjs::runjs("$('.tab-sidebar:eq(1)').trigger('click');")
        shinyjs::runjs("$('.settings-lock').trigger('mouse');")
        shinyjs::runjs("$('#settings-container').trigger('mouseover');") # initialize css
        shinyjs::runjs("$('.settings-lock').trigger('click');") # lock settings tab

        r_local$selected_row <- example_row
        r_local$found_example_trigger <- TRUE
      }
    })


    ## ================================================================================
    ## Modules
    ## ================================================================================
    loading_tsne_server("tsne", watermark = FALSE)

    pgxtable <- loading_table_datasets_server(
      "pgxtable",
      pgxTable_data = pgxTable_data
    )

    ## -----------------------------------------------------------------------------
    ## Description
    ## -----------------------------------------------------------------------------

    shiny::observeEvent(input$module_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Load Dataset</strong>"),
        shiny::HTML(module_infotext),
        easyClose = TRUE, size = "xl"
      ))
    })

    module_infotext <- paste0(
      "This panel shows the available datasets within the platform. The table
        reports a brief description as well as the total number of samples,
        genes, gene sets (or pathways), corresponding phenotypes and the creation
        date.<br><br><b>Selecting the dataset:</b> Users can select a dataset in
        the table. The Dataset info shows the information of the dataset of
        interest and users can load the data by clicking the 'Load dataset'
        button.<br><br><br><center><iframe width='560' height='315'
        src='https://www.youtube.com/embed/elwT6ztt3Fo'
        title='YouTube video player' frameborder='0'
        allow='accelerometer; autoplay; clipboard-write; encrypted-media;
        gyroscope; picture-in-picture' allowfullscreen></iframe><center>"
    )

    ## -----------------------------------------------------------------------------
    ## User interface
    ## -----------------------------------------------------------------------------

    currentSection <- shiny::reactive({
      cdata <- session$clientData
      sub("section-", "", cdata[["url_hash"]])
    })

    output$rowselected <- shiny::reactive({
      !is.null(selectedPGX()) && length(selectedPGX()) > 0
    })
    shiny::outputOptions(output, "rowselected", suspendWhenHidden = FALSE)

    observe({
      df <- getPGXINFO()
      datatypes <- sort(setdiff(df$datatype, c(NA, "")))
      organisms <- sort(setdiff(df$organism, c(NA, "")))
      shiny::updateCheckboxGroupInput(session, "flt_datatype", choices = datatypes)
      shiny::updateCheckboxGroupInput(session, "flt_organism", choices = organisms)
    })

    ## -----------------------------------------------------------------------------
    ## READ initial PGX file info
    ## -----------------------------------------------------------------------------

    ## reactive value for updating table
    reload_pgxdir <- shiny::reactiveVal(0)

    getPGXDIR <- shiny::reactive({
      reload_pgxdir() ## force reload

      email <- auth$email()
      email <- gsub(".*\\/", "", email)
      pdir <- pgx_dir ## from module input

      ## USERDIR=FALSE
      if (enable_userdir) {
        pdir <- paste0(pdir, "/", email)
        if (!is.null(email) && !is.na(email) && email != "") pdir <- paste0(pdir, "/")
        if (!dir.exists(pdir)) {
          dir.create(pdir)
          file.copy(file.path(pgx_dir, "example-data.pgx"), pdir)
        }
      }
      pdir
    })

    getPGXINFO <- shiny::reactive({
      req(auth)
      if (!auth$logged()) {
        warning("[LoadingBoard:getPGXINFO] user not logged in!")
        return(NULL)
      }
      info <- NULL
      pdir <- getPGXDIR()
      info <- pgx.scanInfoFile(pdir, file = "datasets-info.csv", verbose = TRUE)
      if (is.null(info)) {
        aa <- rep(NA, 9)
        names(aa) <- c(
          "dataset", "datatype", "description", "nsamples",
          "ngenes", "nsets", "conditions", "organism", "date"
        )
        info <- data.frame(rbind(aa))[0, ]
      }
      info
    })

    getFilteredPGXINFO <- shiny::reactive({
      ## get the filtered table of pgx datasets
      req(auth)
      if (!auth$logged()) {
        warning("[LoadingBoard:getFilteredPGXINFO] user not logged in!
                    not showing table!")
        return(NULL)
      }
      df <- getPGXINFO()
      if (is.null(df)) {
        return(NULL)
      }

      pgxdir <- getPGXDIR()
      pgxfiles <- dir(pgxdir, pattern = ".pgx$")
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
        rownames(df) <- nrow(df):1
      }

      kk <- unique(c(
        "dataset", "description", "datatype", "nsamples",
        "ngenes", "nsets", "conditions", "date", "organism"
      ))
      kk <- intersect(kk, colnames(df))
      df <- df[, kk, drop = FALSE]
      df
    })

    selectedPGX <- shiny::reactive({
      req(pgxtable)
      sel <- r_local$selected_row
      if (is.null(sel) || length(sel) == 0) {
        return(NULL)
      }
      df <- getFilteredPGXINFO()
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      pgxfile <- as.character(df$dataset[sel])
      pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx") ## add/replace .pgx
      pgxfile
    })

    ## =============================================================================
    ## ========================== OBSERVE/REACT ====================================
    ## =============================================================================

    loadPGX <- function(pgxfile) {
      req(auth$logged())
      if (!auth$logged()) {
        return(NULL)
      }

      pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx") ## add/replace .pgx
      pgxdir <- getPGXDIR()

      pgx.path <- pgxdir[file.exists(file.path(pgxdir, pgxfile))][1]
      pgxfile1 <- file.path(pgx.path, pgxfile)

      pgx <- NULL
      if (file.exists(pgxfile1)) {
        shiny::withProgress(message = "Loading data...", value = 0.33, {
          pgx <- local(get(load(pgxfile1, verbose = 0))) ## override any name
        })
      } else {
        warning("[LoadingBoard::loadPGX] ***ERROR*** file not found : ", pgxfile)
        return(NULL)
      }
      if (!is.null(pgx)) {
        pgx$name <- pgxfile
        return(pgx)
      } else {
        warning("[LoadingBoard::loadPGX] ERROR loading pgx object")
        return(NULL)
      }
    }

    output$downloadpgx <- shiny::downloadHandler(
      ## filename = "userdata.pgx",
      filename = function() {
        selectedPGX()
      },
      content = function(file) {
        pgxfile <- selectedPGX()
        if (is.null(pgxfile) || pgxfile == "" || length(pgxfile) == 0) {
          return(NULL)
        }
        pgx <- loadPGX(pgxfile)
        temp <- tempfile()
        save(pgx, file = temp)
        file.copy(temp, file)
      }
    )

    output$downloadzip <- shiny::downloadHandler(
      ## filename = "userdata.zip",
      filename = function() {
        sub("pgx$", "zip", selectedPGX())
      },
      content = function(file) {
        pgxfile <- selectedPGX()
        if (is.null(pgxfile) || pgxfile == "" || length(pgxfile) == 0) {
          return(NULL)
        }
        pgxname <- sub("[.]pgx$", "", pgxfile)
        pgx <- loadPGX(pgxfile)
        dir.create(tmp <- tempfile())
        tmp2 <- file.path(tmp, pgxname)
        dir.create(tmp2)

        exp.matrix <- sign(pgx$model.parameters$exp.matrix)
        exp.matrix <- contrastAsLabels(exp.matrix) ## new recommended style
        exp.matrix[is.na(exp.matrix)] <- ""

        write.csv(round(pgx$counts, digits = 2), file = file.path(tmp2, "counts.csv"))
        write.csv(pgx$samples, file = file.path(tmp2, "samples.csv"))
        write.csv(exp.matrix, file = file.path(tmp2, "contrasts.csv"))
        write.csv(round(pgx$X, digits = 4), file = file.path(tmp2, "normalized.csv"))

        zipfile <- tempfile(fileext = ".zip")
        zip::zip(zipfile,
          files = paste0(pgxname, "/", c("counts.csv", "samples.csv", "contrasts.csv", "normalized.csv")),
          root = tmp
        )
        file.copy(zipfile, file)
        remove(pgx)
        gc()
      }
    )

    shiny::observeEvent(input$deletebutton, {
      pgxfile <- selectedPGX()
      if (is.null(pgxfile) || pgxfile == "" || length(pgxfile) == 0) {
        return(NULL)
      }

      pgx.path <- getPGXDIR()
      pgxfile1 <- file.path(pgx.path, pgxfile)
      pgxfile1
      sel <- NULL

      deletePGX <- function() {
        if (input$confirmdelete) {
          cat(">>> deleting", pgxfile, "\n")
          pgxfile2 <- paste0(pgxfile1, "_") ## mark as deleted
          file.rename(pgxfile1, pgxfile2)
          reload_pgxdir(reload_pgxdir() + 1)
        } else {
          cat(">>> deletion cancelled\n")
        }
      }

      not.anonymous <- !is.na(auth$name()) && auth$name() != ""
      allow.delete <- !not.anonymous

      allow.delete <- TRUE
      if (!allow.delete) {
        warning(
          "[LoadingBoard::@deletebutton] WARNING:: ", pgxfile,
          " not owned by ", auth$name(), " \n"
        )
        shinyalert::shinyalert(
          title = "Error!",
          text = "You do not have permission to delete this dataset",
          type = "error"
        )
      } else {
        shinyalert::shinyalert(
          "Delete this dataset?",
          paste("Are you sure you want\nto delete '", pgxfile, "'?"),
          confirmButtonText = "Delete",
          showCancelButton = TRUE,
          callbackR = deletePGX,
          inputId = "confirmdelete"
        )
      }
    })


    ## ================================================================================
    ## ========================== LOAD DATA FROM LIST =================================
    ## ================================================================================


    load_react <- reactive({
      btn <- input$loadbutton
      btn2 <- r_local$found_example_trigger
      query <- parseQueryString(session$clientData$url_search)
      logged <- isolate(auth$logged()) ## avoid reloading when logout/login
      (!is.null(btn) || !is.null(query[["pgx"]])) && logged
    })

    shiny::observeEvent(load_react(), {
      if (!load_react()) {
        return(NULL)
      }

      on.exit({
        session$sendCustomMessage(
          "show-tabs",
          list()
        )
      })

      pgxfile <- NULL

      ## Observe URL query
      query <- parseQueryString(session$clientData$url_search)
      if (!is.null(query[["pgx"]])) {
        pgxfile <- query[["pgx"]]
        pgxfile <- basename(pgxfile) ## for security
        pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx") ## add/replace .pgx
      }

      ## Observe button press (over-rides URL query)
      btn <- input$loadbutton
      if (!is.null(btn) && btn != 0) {
        pgxfile <- selectedPGX()
      }
      ## Observe "try example dataset" press
      if (!is.null(r_local$found_example_trigger)) {
        pgxfile <- selectedPGX()
      }

      ## check if file is there
      if (is.na(pgxfile) || is.null(pgxfile) || pgxfile == "" || length(pgxfile) == 0) {
        message("[LoadingBoard@load_react] ERROR file not found : ", pgxfile, "\n")
        return(NULL)
      }

      ## During loading show loading pop-up modal
      pgx.showCartoonModal()

      ## ---------------------------------------------------------------------
      ## ----------------- Loaded PGX object ---------------------------------
      ## ---------------------------------------------------------------------

      loaded_pgx <- loadPGX(pgxfile)
      if (is.null(loaded_pgx)) {
        warning("[LoadingBoard@load_react] ERROR loading PGX file ", pgxfile, "\n")
        beepr::beep(10)
        shiny::removeModal()
        return(NULL)
      }

      ## ----------------- update PGX object ---------------------------------
      loaded_pgx <- pgx.initialize(loaded_pgx)

      if (is.null(loaded_pgx)) {
        warning("[LoadingBoard@load_react] ERROR in object initialization\n")
        beepr::beep(10)
        shiny::showNotification("ERROR in object initialization!\n")
        shiny::removeModal()
        return(NULL)
      }
      loaded_pgx$name <- sub("[.]pgx$", "", pgxfile) ## always use filename

      ## ----------------- update input --------------------------------------
      loadedDataset(loadedDataset() + 1) ## notify new data uploaded

      ## ***NEW*** update PGX from session
      if (1) {
        ## *** EXPERIMENTAL ***. Copying to pgx list to reactiveValues in
        ## session environment.
        for (i in 1:length(loaded_pgx)) {
          pgx[[names(loaded_pgx)[i]]] <- loaded_pgx[[i]]
        }
      }

      ## ----------------- remove modal on exit?? -------------------------
      remove(loaded_pgx)
      gc()
    })


    ## ================================================================================
    ## Header
    ## ================================================================================

    pgx_stats <- reactive({
      pgx <- getFilteredPGXINFO()
      shiny::req(pgx)
      ndatasets <- nrow(pgx)
      nsamples <- sum(as.integer(pgx$nsamples), na.rm = TRUE)
      paste(ndatasets, "Data sets &nbsp;&nbsp;&nbsp;", nsamples, "Samples")
    })

    output$navheader <- shiny::renderUI({
      fillRow(
        flex = c(NA, NA, 1),
        shiny::div(
          id = "navheader-current-section",
          HTML("Load dataset &nbsp;"),
          shiny::actionLink(
            ns("module_info"), "",
            icon = shiny::icon("info-circle"),
            style = "color: #ccc;"
          )
        ),
        shiny::div(HTML(pgx_stats()), id = "navheader-dataset-stats"),
        shiny::br()
      )
    })

    ## ================================================================================
    ## Data sets table
    ## ================================================================================

    ## reactive value for updating table
    touchtable <- shiny::reactiveVal(0)

    andothers <- function(s, split = " ", n = 8) {
      if (is.na(s)) {
        return("")
      }
      s <- sub("^[ ]*", "", s)
      s <- sub("[ ]+", " ", s)
      s1 <- strsplit(s, split = split)[[1]]
      if (length(s1) <= n) {
        return(s)
      }
      n2 <- setdiff(length(s1), n)
      paste(paste(head(s1, n), collapse = " "), "(+", n2, "others)")
    }

    pgxTable_data <- shiny::reactive({
      reload_pgxdir()

      df <- getFilteredPGXINFO()
      shiny::req(df)

      df$dataset <- gsub("[.]pgx$", " ", df$dataset)
      df$conditions <- gsub("[,]", " ", df$conditions)
      df$conditions <- sapply(as.character(df$conditions), andothers, split = " ", n = 5)
      df$description <- shortstring(as.character(df$description), 200)
      df$nsets <- NULL
      df$organism <- NULL
      df
    })

    ## ------------------------------------------------
    ## Board return object
    ## ------------------------------------------------
    res <- list(
      loaded = loadedDataset,
      auth = auth
    )
    return(res)
  })
}
