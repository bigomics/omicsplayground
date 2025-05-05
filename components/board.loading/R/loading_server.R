##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

LoadingBoard <- function(id,
                         pgx,
                         auth,
                         limits = c(
                           "samples" = 1000, "comparisons" = 20,
                           "genes" = 20000, "genesets" = 10000,
                           "datasets" = 10
                         ),
                         pgx_topdir,
                         load_example,
                         reload_pgxdir,
                         current_page,
                         load_uploaded_data,
                         recompute_pgx,
                         new_upload) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE



    reload_pgxdir_public <- reactiveVal(0)
    refresh_shared <- reactiveVal(0)
    is_data_loaded <- reactiveVal(NULL)

    ## static, not changing
    pgx_shared_dir <- stringr::str_replace_all(pgx_topdir, c("data" = "data_shared"))
    pgx_public_dir <- stringr::str_replace_all(pgx_topdir, c("data" = "data_public"))
    enable_public_tabpanel <- dir.exists(pgx_public_dir)

    ## -------------------------------------------------------------------
    ## Received/shared UI
    ## -------------------------------------------------------------------

    pgxreceived <- upload_module_received_server(
      id = "received",
      auth = auth,
      pgx_shared_dir = pgx_shared_dir,
      ##      max_datasets = auth$options$MAX_DATASETS,  ## wrong and not needed...
      reload_pgxdir = reload_pgxdir,
      current_page = current_page
    )

    pgxshared <- upload_module_shared_server(
      id = "shared",
      auth = auth,
      pgx_shared_dir = pgx_shared_dir,
      sendShareMessage = sendShareMessage,
      current_page = current_page,
      refresh = refresh_shared
    )

    output$sharing_alert <- renderUI({
      received_files <- pgxreceived$getReceivedFiles()
      shared_files <- pgxshared$getSharedFiles()
      num_received <- length(received_files)
      num_shared <- length(shared_files)

      no_sharing1 <- !auth$options$ENABLE_USER_SHARE
      no_sharing2 <- (num_received == 0 && num_shared == 0)
      no_sharing <- no_sharing1 || no_sharing2

      if (no_sharing) {
        tag <- bs_alert(HTML("This table shows the <b>available datasets</b> in your library. The <b>Signature t-SNE</b> shows similarity clustering of signatures using t-SNE. Select a dataset in the table and click the <b>Load selected</b> button below."))
        return(tag)
      }

      ## If not show alerts for sharing
      msg <- c()
      if (num_received > 0) {
        msg <- paste("You have received <strong>", num_received, "datasets</strong> that you need to accept.")
      }
      if (num_shared > 0) {
        msg1 <- paste("You have still <strong>", num_shared, "shared datasets</strong> waiting in the queue.")
        msg <- c(msg, msg1)
      }
      bs_alert(
        style = "warning",
        conditional = FALSE,
        shiny::HTML(paste(msg, "Please check the Sharing panel."))
      )
    })

    output$sharing_panel_ui <- renderUI({
      if (!auth$options$ENABLE_USER_SHARE) {
        return(
          "Your version does not allow sharing of datasets."
        )
      }
      received_files <- pgxreceived$getReceivedFiles()
      shared_files <- pgxshared$getSharedFiles()
      num_received <- length(received_files)
      num_shared <- length(shared_files)

      ## if (num_received == 0 && num_shared == 0) {
      ##   dbg("[sharing_panel_ui] no shared datasets!")
      ##   return(paste("No shared datasets in queue."))
      ## }

      out1 <- shiny::wellPanel(
        shiny::HTML("<b>Received datasets.</b> Accept or refuse the received dataset using the action buttons on the right."),
        br(), br(),
        pgxreceived$receivedPGXtable(),
        br()
      )

      out2 <- shiny::wellPanel(
        shiny::HTML("<b>Shared datasets.</b> Resend a message to the receiver or cancel sharing using the action buttons on the right."),
        br(), br(),
        pgxshared$sharedPGXtable(),
        br()
      )

      out <- shiny::tagList(out1, out2)
      return(out)
    })

    ## ======================================================================
    ## LOAD EXAMPLE TRIGGER
    ## ======================================================================
    observeEvent(load_example(),
      {
        # get the row which corresponds to "example-data"
        data_names <- as.character(pgxtable$data()$dataset)
        example_row <- which(data_names == "example-data")[1]
        has.exampledata <- ("example-data" %in% data_names)

        # if not found, throw error modal that example-data doesnt exist
        ## if (is.na(example_row)) {
        if (!has.exampledata) {
          shinyalert::shinyalert(
            title = "No example data",
            text = "Sorry, the example dataset could not be found.",
            type = "warning",
            closeOnClickOutside = FALSE
          )
          return(NULL)
        } else {
          loadAndActivatePGX("example-data")

          # open the left & right sidebar
          bigdash.openSettings(lock = TRUE)
          bigdash.openSidebar()
          bigdash.selectTab(session, selected = "dataview-tab")
          ## shiny::removeModal()
        }
      },
      ignoreInit = TRUE
    )

    ## ================================================================================
    ## Modules
    ## ================================================================================

    pgxtable <- loading_table_datasets_server(
      id = "pgxtable",
      pgx_shared_dir = pgx_shared_dir,
      pgx_public_dir = pgx_public_dir,
      pgx_topdir = pgx_topdir,
      auth = auth,
      loadAndActivatePGX = loadAndActivatePGX,
      loadPGX = loadPGX,
      refresh_shared = refresh_shared,
      reload_pgxdir_public = reload_pgxdir_public,
      reload_pgxdir = reload_pgxdir,
      recompute_pgx = recompute_pgx,
      loadbutton = reactive(input$loadbutton),
      new_upload = new_upload
    )

    loading_tsne_server(
      id = "tsne",
      pgx.dirRT = reactive(auth$user_dir),
      info.table = reactive(pgxtable$data()),
      r_selected = reactive(pgxtable$rows_all()),
      watermark = WATERMARK
    )


    if (enable_public_tabpanel) {
      pgxtable_public <- loading_table_datasets_public_server(
        id = "pgxtable_public",
        pgx_public_dir = pgx_public_dir,
        reload_pgxdir_public = reload_pgxdir_public,
        auth = auth,
        reload_pgxdir = reload_pgxdir
      )

      loading_tsne_server(
        id = "tsne_public",
        pgx.dir = reactive(pgx_public_dir),
        info.table = reactive(pgxtable_public$data()),
        r_selected = reactive(pgxtable_public$rows_all()),
        watermark = WATERMARK
      )
    }

    ## -----------------------------------------------------------------------------
    ## Description
    ## -----------------------------------------------------------------------------

    shiny::observeEvent(input$module_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Loading a dataset from your library</strong>"),
        shiny::HTML(module_infotext),
        easyClose = TRUE, size = "xl"
      ))
    })

    module_infotext <- tspan(paste0(
      "This panel shows the available datasets within the platform. The table
        reports a brief description as well as the total number of samples,
        genes, gene sets (or pathways), corresponding phenotypes and the creation
        date.<br><br><b>Selecting the dataset:</b> Users can select a dataset in
        the table. The Dataset info shows the information of the dataset of
        interest and users can analyze the data by clicking the 'Analyze dataset'
        button.<br><br><br><center><iframe width='560' height='315'
        src='https://www.youtube.com/embed/elwT6ztt3Fo'
        title='YouTube video player' frameborder='0'
        allow='accelerometer; autoplay; clipboard-write; encrypted-media;
        gyroscope; picture-in-picture' allowfullscreen></iframe><center>"
    ), js = FALSE)
    module_infotext <- paste0(
      "<center><iframe width='1120' height='630'
        src='https://www.youtube.com/embed/elwT6ztt3Fo'
        title='YouTube video player' frameborder='0'
        allow='accelerometer; autoplay; clipboard-write; encrypted-media;
        gyroscope; picture-in-picture' allowfullscreen></iframe><center>"
    )


    ## =============================================================================
    ## ========================== OBSERVE/REACT ====================================
    ## =============================================================================

    ## =========================== BUTTON ACTIONS =============================
    ## disable button if no row is selected
    observeEvent(pgxtable$rows_selected(),
      {
        shiny::req(pgxtable)
        if (is.null(pgxtable$rows_selected())) {
          shinyjs::disable(id = "loadbutton")
        } else {
          shinyjs::enable(id = "loadbutton")
        }
      },
      ignoreNULL = FALSE
    )

    loadPGX <- function(pgxfile) {
      req(auth$logged)
      if (!auth$logged) {
        return(NULL)
      }

      pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx") ## add/replace .pgx
      pgxfile1 <- file.path(auth$user_dir, pgxfile)

      pgx <- NULL
      if (file.exists(pgxfile1)) {
        pgx <- playbase::pgx.load(pgxfile1)
      } else {
        warning("[LoadingBoard::loadPGX] ***ERROR*** file not found : ", pgxfile)
        return(NULL)
      }
      if (!is.null(pgx)) {
        pgx$name <- pgxfile
        return(pgx)
      } else {
        warning("[LoadingBoard::loadPGX] ***ERROR*** loading pgx object")
        return(NULL)
      }
    }

    savePGX <- function(pgx, file) {
      req(auth$logged)
      if (!auth$logged) {
        warning("[LoadingBoard::savePGX] ***ERROR*** not logged in or authorized")
        return(NULL)
      }
      file <- paste0(sub("[.]pgx$", "", file), ".pgx") ## add/replace .pgx
      pgxdir <- auth$user_dir
      if (dir.exists(pgxdir)) {
        file1 <- file.path(pgxdir, file)
        playbase::pgx.save(pgx, file = file1)
      } else {
        warning("[LoadingBoard::savePGX] ***ERROR*** pgxdir not found : ", pgxdir)
      }
      return(NULL)
    }

    loadAndActivatePGX <- function(pgxfile) {
      ## During loading show loading pop-up modal
      pgx.showCartoonModal()

      loaded_pgx <- loadPGX(pgxfile)
      if (is.null(loaded_pgx)) {
        warning("[LoadingBoard@load_react] ERROR loading PGX file ", pgxfile, "\n")
        beepr::beep(10)
        shiny::removeModal()
        return(NULL)
      }

      ## ----------------- update PGX object ---------------------------------
      slots0 <- names(loaded_pgx)
      shiny::withProgress(message = "Initializing. Please wait...", value = 0.33, {
        loaded_pgx <- playbase::pgx.initialize(loaded_pgx)

        if (is.null(loaded_pgx)) {
          warning("[loading_server.R@load_react] ERROR in object initialization\n")
          beepr::beep(10)
          shiny::showNotification("ERROR in object initialization!\n")
          shiny::removeModal()
          return(NULL)
        }
        loaded_pgx$name <- sub("[.]pgx$", "", pgxfile) ## always use filename

        ## if PGX object has been updated with pgx.initialize, we save
        ## the updated object.
        slots1 <- names(loaded_pgx)
        if (length(slots1) != length(slots0)) {
          info("[loading_server.R] saving updated PGX")
          new_slots <- setdiff(slots1, slots0)
          savePGX(loaded_pgx, file = pgxfile)
        }

        ## Copying to pgx list to reactiveValues in
        ## session environment.
        info("[loading_server.R] copying pgx object to global environment")
        empty.slots <- setdiff(names(pgx), names(loaded_pgx))
        isolate({
          for (e in empty.slots) {
            pgx[[e]] <- NULL
          }
          for (i in 1:length(loaded_pgx)) {
            pgx[[names(loaded_pgx)[i]]] <- loaded_pgx[[i]]
          }
        })
      }) ## end of withProgress

      info("[loading_server.R] copying pgx done!")
      gc()
      remove(loaded_pgx)

      ## remove modal on exit??

      ## shiny::removeModal()
      bigdash.showTabsGoToDataView(session) ## in ui-bigdashplus.R

      ## notify new data uploaded
      if (is.null(is_data_loaded())) {
        is_data_loaded(1)
      } else {
        is_data_loaded(is_data_loaded() + 1)
      }
    }
    observeEvent(input$newuploadbutton, {
      new_upload(new_upload() + 1)
    })

    observeEvent(load_uploaded_data(), {
      upload_pgx <- sub("[.]pgx$", "", load_uploaded_data())
      loadAndActivatePGX(upload_pgx)
      load_uploaded_data(NULL)
    })

    # Generate report server module

    DatasetReportServer(id = "generate_report", auth = auth, pgxtable = pgxtable)

    ## ================================================================================
    ## Header
    ## ================================================================================

    pgx_stats <- reactive({
      pgx_info <- pgxtable$data()
      shiny::req(pgx_info)
      ndatasets <- nrow(pgx_info)
      nsamples <- sum(as.integer(pgx_info$nsamples), na.rm = TRUE)
      FC.file <- file.path(auth$user_dir, "datasets-allFC.csv")
      if (file.exists(FC.file)) {
        contrasts <- get_contrasts_from_user(auth)
        ncontrasts <- sum(contrasts, na.rm = TRUE)
        return(
          paste(ndatasets, "Data sets &nbsp;&nbsp;&nbsp;", nsamples, "Samples &nbsp;&nbsp;&nbsp;", ncontrasts, "Comparisons")
        )
      } else {
        return(
          paste(ndatasets, "Data sets &nbsp;&nbsp;&nbsp;", nsamples, "Samples")
        )
      }
    })

    output$pgx_stats_ui <- shiny::renderUI(HTML(pgx_stats()))

    ## ================================================================================
    ## Data sets table
    ## ================================================================================

    ## reactive value for updating table
    touchtable <- shiny::reactiveVal(0)


    ## ------------------------------------------------
    ## Board return object
    ## ------------------------------------------------
    res <- list(
      is_data_loaded = is_data_loaded
    )
    return(res)
  })
}
