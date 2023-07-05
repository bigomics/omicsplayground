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
    
    ## reactive variables used only within this module
    rl <- reactiveValues(
      reload_pgxdir_public = 0,
      pgxTablePublic_data = NULL,
      selected_row_public = NULL,
      share_pgx = NULL,
      share_public_pgx = NULL,
      refresh_shared = 0
    )

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
      refresh = reactive(rl$refresh_shared)
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

    ##======================================================================
    ## PUBLIC DATASET FOLDER
    ##======================================================================
    if(enable_public_tabpanel) {

      ##-------------------------------------------------------------------
      ## observers
      ##-------------------------------------------------------------------

      observeEvent( pgxtable_public$rows_selected(), {
        rl$selected_row_public <- pgxtable_public$rows_selected()
      }, ignoreNULL = FALSE)

      # disable buttons when no row is selected; enable when one is selected
      observeEvent( rl$selected_row_public, {
        if (is.null(rl$selected_row_public)) {
          shinyjs::disable(id = 'importbutton')
        } else {
          shinyjs::enable(id = 'importbutton')
        }
      }, ignoreNULL = FALSE)

      observeEvent( input$importbutton, {
          selected_row <- rl$selected_row_public
          pgx_name <- rl$pgxTablePublic_data[selected_row, 'dataset']
          pgx_file <- file.path(pgx_public_dir, paste0(pgx_name, '.pgx'))
          pgx_path <- getPGXDIR()
          new_pgx_file <- file.path(pgx_path, paste0(pgx_name, '.pgx'))

          ## check number of datasets. If deletion is disabled, we count also .pgx_ files... :)
          numpgx <- length(dir(pgx_path, pattern="*.pgx$"))
          if(!enable_delete) numpgx <- length(dir(pgx_path, pattern="*.pgx$|*.pgx_$"))
          maxpgx <- as.integer(limits['datasets'])
          if(numpgx >= maxpgx) {
            ## should use sprintf or glue here...
            msg = "You have reached your datasets limit. Please delete some datasets, or <a href='https://events.bigomics.ch/upgrade' target='_blank'><b><u>UPGRADE</u></b></a> your account."
            shinyalert::shinyalert(
                          title = "Your storage is full",
                          text = HTML(msg),
                          html = TRUE,                          
                          type = "warning"
                        )
            return(NULL)
          }

          if(file.exists(new_pgx_file)) {
            shinyalert::shinyalert(
              title = "Oops! File exists...",
              paste('There is already a dataset called', pgx_name,
                'in your dataset folder. Please delete your file first.')
            )
            return()
          }
          
          ## Copy the file from Public folder to user folder
          shiny::withProgress(message = "Importing dataset...", value = 0.33, {
            file.copy(from = pgx_file, to = new_pgx_file)
            playbase::pgxinfo.updateDatasetFolder(pgx_path, update.sigdb=FALSE)
          })
                    
          shinyalert::shinyalert(
            "Dataset imported",
            paste('The public dataset', pgx_name, 'has now been successfully imported',
              'to your library. Feel free to load it as usual!')
          )

          r_global$reload_pgxdir <- r_global$reload_pgxdir + 1
        }
      )

      ##-------------------------------------------------------------------
      ## make a pgx public (i.e. share publicly)
      ##-------------------------------------------------------------------
      observeEvent(
        rl$share_public_pgx,
        {
          selected_row <- as.numeric(stringr::str_split(rl$share_public_pgx, "_row_")[[1]][2])
          pgx_name <- pgxtable_data()[selected_row, "dataset"]

          alert_val <- shinyalert::shinyalert(
            inputId = 'share_public_confirm',
            title = paste("Share this dataset?"),
            paste('Your dataset', pgx_name, 'will be copied',
              'to the public folder. Other users will be able import and explore',
              'this dataset. Are you sure?'),
            html = TRUE,
            showCancelButton = TRUE,
            showConfirmButton = TRUE
          )
        },
        ignoreNULL = TRUE
      )

      observeEvent(input$share_public_confirm, {

        if (input$share_public_confirm) {

          selected_row <- as.numeric(stringr::str_split(rl$share_public_pgx, "_row_")[[1]][2])
          pgx_name <- pgxtable_data()[selected_row, "dataset"]
          pgx_name <- sub("[.]pgx$", "", pgx_name)
          pgx_path <- getPGXDIR()
          pgx_file <- file.path(pgx_path, paste0(pgx_name, '.pgx'))
          new_pgx_file <- file.path(pgx_public_dir, paste0(pgx_name, '.pgx'))
          pgx_file <- file.path(pgx_path, paste0(pgx_name, ".pgx"))

          new_pgx_file <- file.path(
            pgx_public_dir,
            paste0(pgx_name, ".pgx")
          )

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
            pgx0  <- playbase::pgx.load(pgx_file)
            unknown.creator <- pgx0$creator %in% c(NA,"","user","anonymous","unknown")
            if("creator" %in% names(pgx0) && !unknown.creator) {
              file.copy(from = pgx_file, to = new_pgx_file)
            } else {
              pgx0$creator <- session$user  ## really?
              if(pgx0$creator %in% c(NA,"","user","anonymous","unknown"))  pgx0$creator <- "unknown"
              playbase::pgx.save(pgx0, file = new_pgx_file)
            }
          })

          rl$reload_pgxdir_public <- rl$reload_pgxdir_public + 1
          ## r_global$reload_pgxdir <- r_global$reload_pgxdir + 1

          shinyalert::shinyalert(
            title = "Successfully shared!",
            paste(
              "Your dataset", pgx_name, "has now been successfully",
              "been publicly shared. Thank you!"
            )
          )
        }

        rl$share_public_pgx <- NULL

      })
    }


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
      rl = rl,
      r_global = r_global,
      loadAndActivatePGX = loadAndActivatePGX,
      loadPGX = loadPGX,
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
        table = reactive(rl$pgxTablePublic_data)
      )

      loading_tsne_server(
        id = "tsne_public",
        pgx.dir = reactive(pgx_public_dir),
        info.table = getFilteredPGXINFO_PUBLIC,
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
      email <- auth$email()
      email <- gsub(".*\\/", "", email)  ##??
      pdir  <- pgx_topdir ## from module input

      ## Append email to the pgx path.
      if (enable_userdir) {
        pdir <- paste0(pdir, "/", email)
        if (!is.null(email) && !is.na(email) && email != "") pdir <- paste0(pdir, "/")

        #If dir not exists, create and copy example pgx file
        if (!dir.exists(pdir)) {
          dir.create(pdir)
          file.copy(file.path(pgx_topdir, "example-data.pgx"), pdir)
        }
      }
      pdir
    })

    
    if(enable_public_tabpanel) {

      getPGXINFO_PUBLIC <- shiny::eventReactive({
        rl$reload_pgxdir_public
      }, {
        req(auth)
        if (!auth$logged()) {
          warning("[LoadingBoard:getPGXINFO_PUBLIC] user not logged in!")
          return(NULL)
        }

        ## update meta files
        info <- NULL
        shiny::withProgress(message = "Checking datasets library...", value = 0.33, {
          need_update <- playbase::pgxinfo.needUpdate(pgx_public_dir, check.sigdb=FALSE)
        })

      if(need_update) {
        pgx.showSmallModal("Updating datasets library<br>Please wait...")
        shiny::withProgress(message = "Updating datasets library...", value = 0.33, {
          ## before reading the info file, we need to update for new files
          playbase::pgxinfo.updateDatasetFolder(pgx_public_dir, update.sigdb=FALSE)

        })
        shiny::removeModal(session)
      }

      info <- playbase::pgxinfo.read(pgx_public_dir, file = "datasets-info.csv")
      return(info)
      })

      getFilteredPGXINFO_PUBLIC <- shiny::reactive({
        ## get the filtered table of pgx datasets
        req(auth)
        if (!auth$logged()) {
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
        df
      })

    }  ## end of if public


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

    if(enable_public_tabpanel) {
      observeEvent(
        ## c(getFilteredPGXINFO_PUBLIC(), rl$reload_pgxdir_public), {
        c(getFilteredPGXINFO_PUBLIC()), {
          df <- getFilteredPGXINFO_PUBLIC()
          df$dataset <- sub("[.]pgx$", "", df$dataset)
          df$conditions <- gsub("[,]", " ", df$conditions)
          df$conditions <- sapply(as.character(df$conditions), andothers, split = " ", n = 5)
          df$description <- playbase::shortstring(as.character(df$description), 200)
          df$nsets <- NULL
          df$organism <- NULL
          rl$pgxTablePublic_data <- df
        }
      )
    }

    ## ------------------------------------------------
    ## Board return object
    ## ------------------------------------------------
    res <- list(
      loaded = reactive(r_global$loadedDataset)
##      auth = auth
    )
    return(res)
  })
}
