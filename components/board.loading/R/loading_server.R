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
                         enable_userdir = TRUE,
                         enable_pgxdownload = TRUE,
                         enable_delete = TRUE,
                         enable_user_share = TRUE,
                         enable_public_share = TRUE,
                         r_global) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    reload_pgxdir_public <- reactiveVal(0)
    refresh_shared <- reactiveVal(0)
    update_received <- reactiveVal(0)

    ## static, not changing
    pgx_shared_dir = stringr::str_replace_all(pgx_topdir, c('data'='data_shared'))
    pgx_public_dir = stringr::str_replace_all(pgx_topdir, c('data'='data_public'))
    enable_public_tabpanel <- dir.exists(pgx_public_dir)

    ## only enable if folder exists
    if(enable_user_share && !dir.exists(pgx_shared_dir)) {
      info("[loading_server.R] share folder does not exist. disabling user sharing.")
      enable_user_share <- FALSE
    }
    if(enable_public_share && !dir.exists(pgx_public_dir)) {
      info("[loading_server.R] public folder does not exist. disabling public sharing.")
      enable_public_share <- FALSE
    }

    ##-------------------------------------------------------------------
    ## Received/shared UI
    ##-------------------------------------------------------------------

    pgxreceived <- upload_module_received_server(
      id = "received",
      auth = auth,
      pgx_shared_dir= pgx_shared_dir,
      getPGXDIR = getPGXDIR,
      max_datasets = limits['datasets'],
      enable_delete = enable_delete,
      r_global = r_global
    )

    pgxshared <- upload_module_shared_server(
      id = "shared",
      auth = auth,
      pgx_shared_dir= pgx_shared_dir,
      sendShareMessage = sendShareMessage,
      r_global = r_global,
      refresh = refresh_shared
    )

    output$sharing_alert <- renderUI({
      received_files <- pgxreceived$getReceivedFiles()
      shared_files   <- pgxshared$getSharedFiles()
      num_received <- length(received_files)
      num_shared <- length(shared_files)

      if(num_received==0 && num_shared==0 ) {
        tag <- bs_alert(HTML("This table shows the <b>available datasets</b> in your library. The table reports a brief description of each dataset. The <b>Signature t-SNE</b> shows similarity clustering of fold-change signatures using t-SNE. Select a dataset in the table and load the data by clicking the <b>Load Dataset</b> button below."))
        return(tag)
      }

      ## If not show alerts for sharing
      msg <- c()
      if(num_received>0) {
        msg <- paste("You have received <strong>",num_received,"datasets</strong> that you need to accept.")
      }
      if(num_shared>0) {
        msg1 <- paste("You have still <strong>",num_shared,"shared datasets</strong> waiting in the queue.")
        msg <- c(msg, msg1)
      }
      bs_alert(
        style = "warning",
        conditional = FALSE,
        shiny::HTML(paste(msg,"Please check the Sharing panel."))
      )
    })

    output$sharing_panel_ui <- renderUI({
      received_files <- pgxreceived$getReceivedFiles()
      shared_files   <- pgxshared$getSharedFiles()
      num_received <- length(received_files)
      num_shared <- length(shared_files)

      if( num_received==0 && num_shared==0 ) {
        return(paste("No datasets being shared."))
      }

      out <- tagList()
      if(length(out)==0) {
        if(num_received>0) {
          out1 <- shiny::wellPanel(
            shiny::HTML("<b>Received datasets.</b> Accept or refuse the received dataset using the action buttons on the right."),
            br(),br(),
            pgxreceived$receivedPGXtable(),
            br()
          )
          out <- tagList(out, out1)
        }
        if(num_shared>0) {
          out2 <- shiny::wellPanel(
            shiny::HTML("<b>Shared datasets.</b> Resend a message to the receiver or cancel sharing using the action buttons on the right."),
            br(),br(),
            pgxshared$sharedPGXtable(),
            br()
          )
          if(length(out)==0) {
            out <- out2
          } else {
            out <- tagList(out, br(), out2)
          }
        }
      }

      return(out)
    })


    ##======================================================================
    ## LOAD EXAMPLE TRIGGER
    ##======================================================================
    observeEvent(r_global$load_example_trigger, {

      # get the row which corresponds to "example-data"
      data_names <- as.character(pgxtable_data()$dataset)
      example_row <- which(data_names == "example-data")[1]
      has.exampledata <- ("example-data" %in% data_names)

      # if not found, throw error modal that example-data doesnt exist
      ##if (is.na(example_row)) {
      if (!has.exampledata) {
        shinyalert::shinyalert(
          title = "No example data",
          text ='Sorry, the example dataset could not be found.',
          type = "warning",
          closeOnClickOutside = FALSE
        )
        return(NULL)
      } else {

        loadAndActivatePGX("example-data")

        # open the left & right sidebar
        bigdash.openSettings(lock=TRUE)
        bigdash.openSidebar()
        bigdash.selectTab(session, selected = 'dataview-tab')
        ##shiny::removeModal()
      }
    }, ignoreInit = TRUE)

    ## ================================================================================
    ## Modules
    ## ================================================================================

    pgxtable <- loading_table_datasets_server(
      id = "pgxtable",
      ## getPGXINFO = getPGXINFO,
      getPGXDIR = getPGXDIR,
      pgx_shared_dir = pgx_shared_dir,
      pgx_topdir = pgx_topdir,
      auth = auth,
      r_global = r_global,
      loadAndActivatePGX = loadAndActivatePGX,
      loadPGX = loadPGX,
      refresh_shared = refresh_shared,
      enable_pgxdownload = enable_pgxdownload,
      enable_delete = enable_delete,
      enable_public_share = enable_public_share,
      enable_user_share = enable_user_share
    )

    pgxtable_data <- reactive({
      shiny::req(pgxtable)
      pgxtable$data()
    })

    loading_tsne_server(
      id = "tsne",
      pgx.dir = getPGXDIR,
      info.table = reactive(pgxtable_data()),
      r_selected = reactive(pgxtable$rows_all()),
      watermark = WATERMARK
    )


    if(enable_public_tabpanel) {

      pgxtable_public <- loading_table_datasets_public_server(
        id = "pgxtable_public",
        getPGXDIR = getPGXDIR,
        pgx_public_dir = pgx_public_dir,
        reload_pgxdir_public = reload_pgxdir_public,
        auth = auth,
        enable_delete = enable_delete,
        limits = limits,
        r_global = r_global
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
    module_infotext <- paste0(
      "<center><iframe width='1120' height='630'
        src='https://www.youtube.com/embed/elwT6ztt3Fo'
        title='YouTube video player' frameborder='0'
        allow='accelerometer; autoplay; clipboard-write; encrypted-media;
        gyroscope; picture-in-picture' allowfullscreen></iframe><center>"
    )

    ## -----------------------------------------------------------------------------
    ## READ initial PGX file info
    ## -----------------------------------------------------------------------------

    ## Get the pgx folder. If user folders are enabled, the user email
    ## is appended to the pgx dirname.
    getPGXDIR <- shiny::reactive({

      ## react on change of auth user
      user <- auth$email()
      if(auth$method == "password") user <- auth$name()
      user <- gsub(".*\\/", "", user)  ## get rid of dangerous characters that can skip folders...
      pdir <- pgx_topdir  ## from module input

      dbg("[LoadingBoard::getPGXDIR] authentication = ",auth$method)
      dbg("[LoadingBoard::getPGXDIR] pgx_topdir = ",pgx_topdir)
      dbg("[LoadingBoard::getPGXDIR] user = ",user)
      valid.user <-(!is.null(user) && !is.na(user) && user != "") 
      
      ## Append email to the pgx path.
      if (valid.user && enable_userdir) {
        pdir <- file.path(pdir, user)
        #If dir not exists, create and copy example pgx file
        example.file <- file.path(pgx_topdir, "example-data.pgx")
        if (!dir.exists(pdir) && file.exists(example.file)) {
          dir.create(pdir)
          file.copy(example.file, pdir)
        }
      }
      dbg("[LoadingBoard::getPGXDIR] user.pgxdir = ",pdir)
      pdir
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
        shiny::withProgress(message = "Loading PGX data...", value = 0.33, {
          pgx <- playbase::pgx.load(pgxfile1)
        })
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
      req(auth$logged())
      if (!auth$logged()) {
        warning("[LoadingBoard::savePGX] ***ERROR*** not logged in or authorized")
        return(NULL)
      }
      file <- paste0(sub("[.]pgx$", "", file), ".pgx") ## add/replace .pgx
      pgxdir <- getPGXDIR()[1]
      if (dir.exists(pgxdir)) {
        file1 <- file.path(pgxdir, file)
        shiny::withProgress(message = "Saving PGX data...", value = 0.33, {
          playbase::pgx.save(pgx, file=file1)
        })
      } else {
        warning("[LoadingBoard::savePGX] ***ERROR*** dir not found : ", pgxdir)
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
      if(length(slots1) != length(slots0)) {
        info("[loading_server.R] saving updated PGX")
        new_slots <- setdiff(slots1, slots0)
        savePGX(loaded_pgx, file=pgxfile)
      }

      ## Copying to pgx list to reactiveValues in
      ## session environment.
      info("[loading_server.R] copying pgx object to global environment")
      empty.slots <- setdiff(names(pgx),names(loaded_pgx))
      isolate({
        for (e in empty.slots) {
          pgx[[e]] <- NULL
        }
        for (i in 1:length(loaded_pgx)) {
          pgx[[names(loaded_pgx)[i]]] <- loaded_pgx[[i]]
        }
      })
      info("[loading_server.R] copying pgx done!")
      gc()
      remove(loaded_pgx)

      ## remove modal on exit??
      ## shiny::removeModal()
      bigdash.showTabsGoToDataView(session)  ## in ui-bigdashplus.R

      ## notify new data uploaded
      r_global$loadedDataset <- r_global$loadedDataset + 1
    }


    observeEvent(r_global$load_data_from_upload, {
        data_names <- as.character(pgxtable_data()$dataset)
        data_names <- sub("[.]pgx$","",data_names)
        upload_pgx <- sub("[.]pgx$","",r_global$load_data_from_upload)
        dbg("[load_data_from_upload] upload_pgx = ",upload_pgx)
        loadAndActivatePGX(upload_pgx)
        ##loadRowManual(loadRowManual() + 1)
        r_global$load_data_from_upload <- NULL
    }, ignoreNULL = TRUE)


    ## ================================================================================
    ## Header
    ## ================================================================================

    pgx_stats <- reactive({
      pgx <- pgxtable_data()
      shiny::req(pgx)
      ndatasets <- nrow(pgx)
      nsamples <- sum(as.integer(pgx$nsamples), na.rm = TRUE)
      paste(ndatasets, "Data sets &nbsp;&nbsp;&nbsp;", nsamples, "Samples")
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
      loaded = reactive(r_global$loadedDataset)
    )
    return(res)
  })
}

# util function
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

