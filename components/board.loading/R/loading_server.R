##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

LoadingBoard <- function(id,
                         pgx_topdir,
                         pgx,
                         auth,
                         limits = c(
                           "samples" = 1000, "comparisons" = 20,
                           "genes" = 20000, "genesets" = 10000,
                           "datasets" = 10
                         ),
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
      selected_row = NULL,
      found_example_trigger = 0,
      pgxTable_data = NULL,
      pgxTable_edited = 0,
      pgxTable_edited_row = NULL,
      pgxTablePublic_data = NULL,
      selected_row_public = NULL,
      download_pgx = NULL,
      download_zip = NULL,
      share_pgx = NULL,
      share_public_pgx = NULL,
      delete_pgx = NULL,
      delete_trigger = 0,
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

    # this allows for deselection (selected_row -> NULL)
    #observeEvent(pgxtable$rows_selected(), {
    #  rl$selected_row <- pgxtable$rows_selected()
    #}, ignoreNULL = FALSE)

    # disable buttons when no row is selected; enable when one is selected
    observeEvent( pgxtable$rows_selected(), {
      shiny::req(pgxtable)
      if (is.null(pgxtable$rows_selected())) {
        shinyjs::disable(id = 'loadbutton')
      } else {
        shinyjs::enable(id = 'loadbutton')
      }
    }, ignoreNULL = FALSE)


    ##-------------------------------------------------------------------
    ## util functions
    ##-------------------------------------------------------------------

    is_valid_email <- function(email) {
        is_personal <- grepl("gmail|ymail|outlook|yahoo|hotmail|mail.com$|icloud|msn.com$",email)
        valid_email <- grepl(".*@.*[.].*",email)
        valid_email <- valid_email && !grepl("[*/\\}{]",email) ## no special chars
        return(!is_personal && valid_email)
    }

    get_coworkers <- function(pgxdir, email) {
        domain <- sub(".*@","",email)
        if(email=="" || domain=="") return(NULL)
        cow <- dir(pgxdir, pattern=paste0(domain,"$"))
        cow <- setdiff(cow, email)
        cow
    }

    sendShareMessage <- function(pgxname, sender, share_user, path_to_creds='gmail_creds') {

      dbg("[sendShareMessage] pgxname = ", pgxname)
      dbg("[sendShareMessage] sender = ", sender)
      dbg("[sendShareMessage] share_user = ", share_user)
      
      if(!file.exists(path_to_creds)) return(NULL)
      
      blastula::smtp_send(
        blastula::compose_email(
          body = blastula::md(
            glue::glue(
              "Hello, {sender} shared a dataset with you on OmicsPlayground! ",
              "Login to accept the new dataset.")
          ),
          footer = blastula::md(
            glue::glue("Email sent on {blastula::add_readable_time()}.")
          )
        ),
        from = "app@bigomics.ch",
        to = share_user,
        subject = paste("Your friend",sender,"shared data on OmicsPlayground"),
        credentials = blastula::creds_file(path_to_creds)
      )
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
            shiny::HTML("<b>Send datasets.</b> Resend a message to the receiver or cancel sharing using the action buttons on the right."),
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

    ##-------------------------------------------------------------------
    ## share dataset with specific user
    ##-------------------------------------------------------------------
    
    share_dialog <- function(pgxname, choices=NULL) {

      select_user <- NULL
      if(!is.null(choices) && length(choices)>0) {
        select_user <- tagList(
          shiny::selectizeInput(ns("share_user2"),
            "Or select a coworker from the list:",
            choices = c("Choose coworker..."="",choices),
            options=list(create=TRUE))
        )
      }

      shiny::modalDialog(
        tagList(
            paste("Your dataset",pgxname,"will be shared with the user below."),
            br(), br(),
            shiny::textInput(ns("share_user"),
              label = "Enter email who will receive the dataset:",
              placeholder = "Type email..."
            ),
            select_user,
            ##shiny::textInput(ns("share_user2"), "Re-enter use email:")
            shiny::textOutput(ns('error_alert')) %>%
              tagAppendAttributes(style = 'color: red;')
        ),
        title = "Share this dataset?",
        footer = tagList(
            actionButton(ns("share_dialog_cancel"), "Cancel"),
            actionButton(ns("share_dialog_confirm"), "Share")
        )
      )
    }

    #' will either take free-form or selected email
    input_share_user <- reactive({
      share_user <- input$share_user
      if(!is.null(input$share_user2) && input$share_user2 != '') {
        share_user <- input$share_user2
      }
      share_user
    })

    observeEvent(input$share_user2,{
        if(input$share_user2!="") {
            updateTextInput(session, "share_user", value='')
        }
    })
    
    selectedPGX2 <- reactive({
      selected_row <- as.numeric(stringr::str_split(rl$share_pgx, "_row_")[[1]][2])
      pgx_name <- pgxtable_data()[selected_row, "dataset"]
      pgx_name
    })

    # put user dataset into shared folder
    observeEvent( rl$share_pgx, {

        ## sharing folder has to exists
        if(!dir.exists(pgx_shared_dir)) {
          shinyalert::shinyalert(
            title = "Oops! Cannot share...",
            text = paste('This server does not support sharing.',
                         'Please contact your administrator.')
          )
          rl$share_pgx <- NULL
          return()
        }

        ## user has to be logged in and have email for them to share with other users
        if (auth$email() == "") {
          shinyalert::shinyalert(
             title = "Oops! You're not logged in...",
             text = paste("You need to be logged in with a valid email",
                          "address to share pgx files with other users.")
          )
          rl$share_pgx <- NULL
          return()
        }

        ## check how many are already in queue
        pp <- paste0("__from__",auth$email(),"__$")
        num_shared_queue <- length(dir(pgx_shared_dir, pattern=pp))

        if(num_shared_queue >= opt$MAX_SHARED_QUEUE) {   ## NB opt is global...
          shinyalert::shinyalert(
            title = "Oops! Too many shared...",
            text = paste("You have already too many shared datasets in the waiting queue.",
                         "Please contact your administrator.")
          )
          rl$share_pgx <- NULL
          return()
        }

        ## show share dialog
        coworkers <- get_coworkers( pgx_topdir, auth$email())
        pgxname <- sub("[.]pgx$","",selectedPGX2())
        if(length(coworkers)==0) coworkers <- NULL
        shiny::showModal(share_dialog(pgxname, choices=coworkers))

    }, ignoreNULL = TRUE)

    observeEvent(input$share_dialog_cancel, {
        rl$share_pgx <- NULL
        output$error_alert <- renderText({''})
        shiny::removeModal()
    })

    observeEvent(input$share_dialog_confirm, {

        share_user <- input_share_user()
        if (share_user == '') {
            output$error_alert <- renderText({'Please enter an email.'})
            return()
        }

        if (!is_valid_email(share_user)) {
            output$error_alert <- renderText({'Email is not valid. Please use only work or business emails.'})
            return()
        }
        if (share_user==auth$email()) {
            output$error_alert <- renderText({'Error. You cannot share with yourself...'})
            return()
        }

        output$error_alert <- renderText({''})

        shiny::removeModal()

        pgx_name <- selectedPGX2()  ## need rl$share_pgx

        alert_val <- shinyalert::shinyalert(
            inputId = "share_confirm",
            title = "Are you sure?",
            tagList(
                paste("Your dataset", pgx_name, "will be shared with",share_user)
            ),
            html = TRUE,
            showCancelButton = TRUE,
            showConfirmButton = TRUE
        )
    })

    
    observeEvent(input$share_confirm, {

        # if confirmed, then share the data
        if (input$share_confirm) {

            pgx_name <- selectedPGX2()
            pgx_name <- sub("[.]pgx$", "", pgx_name)
            pgx_path <- getPGXDIR()
            pgx_file <- file.path(pgx_path, paste0(pgx_name, ".pgx"))

            dbg("[LoadingBoard:observeEvent(input$share_confirm)] pgx_name = ",pgx_name)
            dbg("[LoadingBoard:observeEvent(input$share_confirm)] pgx_path = ",pgx_path)                    
            dbg("[LoadingBoard:observeEvent(input$share_confirm)] pgx_file = ",pgx_file)

            ## The shared file will be copied to the data_shared
            ## folder with the name of the sender and receiver in the
            ## file name.
            share_user <- input_share_user()
            new_pgx_file <- file.path(
                pgx_shared_dir,
                paste0(pgx_name, ".pgx", '__to__', share_user,
                     '__from__', auth$email(), '__')
            )

            if(file.exists(new_pgx_file)) {
                shinyalert::shinyalert(
                    title = "Oops! File exists...",
                    paste('There is already a dataset called', pgx_name,
                          'being shared. Please rename your file.')
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
                pgx0$creator <- session$user ## really?
                if (pgx0$creator %in% c(NA, "", "user", "anonymous", "unknown")) pgx0$creator <- "unknown"
                playbase::pgx.save(pgx0, file = new_pgx_file)
              }

              ## write transaction to log file
              log.entry <- data.frame( date=date(), from="jane@demo.com", to="tarzan@demo.com",file="example-data.pgx")
              log.file <- file.path(pgx_shared_dir, "PGXSHARE-TRANSACTIONS.log")
              log.entry <- data.frame( date=date(), from=auth$email(), to=share_user, file=paste0(pgx_name,".pgx"))
              if(file.exists(log.file)) {
                write.table(log.entry, file=log.file, col.names=FALSE, row.names=FALSE, sep=',', append=TRUE)
              } else {
                write.table(log.entry, file=log.file, col.names=TRUE, row.names=FALSE, sep=',')
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
            sender <- auth$email()
            sendShareMessage(pgx_name, sender, share_user, path_to_creds='gmail_creds')

            rl$refresh_shared <- rl$refresh_shared + 1
          
        }

        rl$share_pgx <- NULL
    })

    ##======================================================================
    ## LOAD EXAMPLE TRIGGER
    ##======================================================================
    observeEvent(r_global$load_example_trigger, {

      data_names <- as.character(pgxtable_data()$dataset)

      # get the row which corresponds to "example-data"
      example_row <- which(data_names == "example-data")[1]

      # if not found, throw error modal that example-data doesnt exist
      if (is.na(example_row)) {
        shinyalert::shinyalert(
          title = "No example data found",
          text ='Sorry, the example dataset cannot be found. You may have deleted
            it in a previous session.',
          type = "warning",
          closeOnClickOutside = FALSE
        )
        return(NULL)
      } else {
        # open the left & right sidebar
        bigdash.openSettings(lock=TRUE)
        bigdash.openSidebar()

        # go to dataview
        bigdash.selectTab(session, selected = 'dataview-tab')

        rl$selected_row <- example_row
        rl$found_example_trigger <- rl$found_example_trigger+1
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

          ## check number of datasets
          numpgx <- length( dir(pgx_path, pattern="*.pgx$") )
          maxpgx <- as.integer(limits['datasets'])
          if(numpgx >= maxpgx) {
            dbg("[observeEvent:importbutton] numpgx = ",numpgx)
            dbg("[observeEvent:importbutton] maxpgx = ",maxpgx)
            ## should use sprintf or glue here...
            msg = "You have NUMPGX datasets in your library and your quota is MAXPGX datasets. Please delete some datasets, or <a href='https://events.bigomics.ch/upgrade'><b><u>UPGRADE</u></b></a> your account."
            msg <- sub("NUMPGX", numpgx, msg)
            msg <- sub("MAXPGX", maxpgx, msg)
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
          file.copy(from = pgx_file, to = new_pgx_file)
          playbase::pgxinfo.updateDatasetFolder(pgx_path, update.sigdb=FALSE)
          r_global$reload_pgxdir <- r_global$reload_pgxdir + 1
          
          shinyalert::shinyalert(
            "Dataset imported",
            paste('The public dataset', pgx_name, 'has now been successfully imported',
              'to your data files. Feel free to load it as usual!'
            )
          )
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
      getPGXINFO = getPGXINFO,
      getPGXDIR = getPGXDIR,
      auth = auth,
      rl = rl,
      r_global = r_global,
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

      loading_tsne_server(
        id = "tsne_public",
        pgx.dir = reactive(pgx_public_dir),
        info.table = getFilteredPGXINFO_PUBLIC,
        r_selected = reactive(pgxtable_public$rows_all()),
        watermark = WATERMARK
      )

      pgxtable_public <- loading_table_datasets_public_server(
        id = "pgxtable_public",
        table = reactive(rl$pgxTablePublic_data)
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
    ## User interface
    ## -----------------------------------------------------------------------------

    output$rowselected <- shiny::reactive({
      !is.null(selectedPGX()) && length(selectedPGX()) > 0
    })
    shiny::outputOptions(output, "rowselected", suspendWhenHidden = FALSE)

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

    getPGXINFO <- shiny::reactive({
      req(auth)
      if (!auth$logged()) {
        return(NULL)
      }

      ## upstream trigger
      r_global$reload_pgxdir  
      pgxdir <- getPGXDIR()
      
      shiny::withProgress(message = "Checking datasets library...", value = 0.33, {
        need_update <- playbase::pgxinfo.needUpdate(pgxdir, check.sigdb=FALSE)
      })
      
      ## before reading the info file, we need to update for new files
      if(need_update) {
        dbg("[loading_server:getPGXINFO] updating library...")
        pgx.showSmallModal("Updating your library<br>Please wait...")
        shiny::withProgress(message = "Updating your library...", value = 0.33, {
          playbase::pgxinfo.updateDatasetFolder(pgxdir, update.sigdb=FALSE)
        })
        shiny::removeModal(session)
      }

      info <- playbase::pgxinfo.read(pgxdir, file = "datasets-info.csv")
      shiny::removeModal(session)
      return(info)
    })

    getFilteredPGXINFO <- shiny::reactive({
      ## get the filtered table of pgx datasets
      req(auth)
      if (!auth$logged()) {
        return(NULL)
      }
      df <- getPGXINFO()
      if (is.null(df)) {
        return(NULL)
      }

      ## filter on actual PGX present in folder
      pgxdir <- getPGXDIR()
      pgxfiles <- dir(pgxdir, pattern = ".pgx$")
      sel <- sub("[.]pgx$", "", df$dataset) %in% sub("[.]pgx$", "", pgxfiles)
      df <- df[sel, , drop = FALSE]

      ## Apply table selection filters
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

    }

    selectedPGX <- shiny::reactive({
      req(pgxtable_data())
      sel <- rl$selected_row
      if (is.null(sel) || length(sel) == 0) {
        return(NULL)
      }
      df <- pgxtable_data()
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


    # DOWNLOAD PGX FILE #
    observeEvent(rl$download_pgx, {
      shinyjs::click(id = 'download_pgx_btn')
      rl$download_pgx_idx <- rl$download_pgx
      rl$download_pgx <- NULL
    }, ignoreNULL = TRUE)

    output$download_pgx_btn <- shiny::downloadHandler(
      ## filename = "userdata.pgx",
      filename = function() {
        sel <- row_idx <- as.numeric(stringr::str_split(rl$download_pgx_idx, '_row_')[[1]][2])
        df <- getFilteredPGXINFO()
        pgxfile <- as.character(df$dataset[sel])
        pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx")
        pgxfile
      },
      content = function(file) {
        sel <- row_idx <- as.numeric(stringr::str_split(rl$download_pgx_idx, '_row_')[[1]][2])
        df <- getFilteredPGXINFO()
        pgxfile <- as.character(df$dataset[sel])
        pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx")
        pgx <- loadPGX(pgxfile)
        temp <- tempfile()
        save(pgx, file = temp)
        file.copy(temp, file)
      }
    )

    # DOWNLOAD DATA AS ZIP FILE #
    observeEvent(rl$download_zip, {
      shinyjs::click(id = 'download_zip_btn')
      rl$download_zip_idx <- rl$download_zip
      rl$download_zip <- NULL
    }, ignoreNULL = TRUE)

    output$download_zip_btn <- shiny::downloadHandler(
      ## filename = "userdata.zip",
      filename = function() {
        sel <- row_idx <- as.numeric(stringr::str_split(rl$download_zip_idx, '_row_')[[1]][2])
        df <- getFilteredPGXINFO()
        pgxfile <- as.character(df$dataset[sel])
        pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx") ## add/replace .pgx
        newfile <- sub("pgx$", "zip", pgxfile)
        newfile
      },
      content = function(file) {
        sel <- row_idx <- as.numeric(stringr::str_split(rl$download_zip_idx, '_row_')[[1]][2])
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
          files = paste0(pgxname, "/", c("counts.csv", "samples.csv",
            "contrasts.csv", "normalized.csv")),
          root = tmp
        )
        file.copy(zipfile, file)
        remove(pgx)
        gc()
      }
    )
    
    shiny::observeEvent(rl$delete_pgx, {
      
      row_idx <- as.numeric(stringr::str_split(rl$delete_pgx, '_row_')[[1]][2])
      df <- getFilteredPGXINFO()
      pgxfile <- as.character(df$dataset[row_idx])
      pgxname <- sub("[.]pgx$", "", pgxfile)
      pgxfile <- paste0(pgxname, ".pgx") ## add/replace .pgx
      
      pgx.path <- getPGXDIR()
      pgxfile1 <- file.path(pgx.path, pgxfile)
      sel <- NULL

      deletePGX <- function(x) {
        if (input$confirmdelete) {
          pgxfile2 <- paste0(pgxfile1, "_") ## mark as deleted
          file.rename(pgxfile1, pgxfile2)
          ## !!!! we should also delete entry in PGXINFO and allFC !!!
          ## playbase::pgx.deleteInfoPGX(pgxinfo, pgxname)
          info <- read.csv(file.path(pgx.path,"datasets-info.csv"),row.names=1)
          idx <- match(pgxname,info$dataset)
          if(length(idx)) {
            info <- info[-idx,]
            write.csv(info, file.path(pgx.path,"datasets-info.csv"))
          }
          r_global$reload_pgxdir <- r_global$reload_pgxdir + 1
        }
      }

      if(enable_delete) {
        shinyalert::shinyalert(
          "Delete this dataset?",
          paste("Are you sure you want\nto delete '", pgxfile, "'?"),
          confirmButtonText = "Delete",
          showCancelButton = TRUE,
          callbackR = deletePGX,
          inputId = "confirmdelete"
        )
      } else {
          shinyalert::shinyalert(
            title = "Oops!",
            text = "Delete is disabled on this server"
          )
      }

      rl$delete_pgx <- NULL

    }, ignoreNULL = TRUE, ignoreInit = TRUE)


    ## ================================================================================
    ## ========================== LOAD DATA FROM LIST =================================
    ## ================================================================================

    load_react <- reactive({
      btn <- input$loadbutton
      btn2 <- rl$found_example_trigger
      logged <- isolate(auth$logged()) ## avoid reloading when logout/login
      !is.null(btn) && logged
    })

    loadRowManual <- reactiveVal(0)

    observeEvent(r_global$load_data_from_upload, {
        data_names <- as.character(pgxtable_data()$dataset)
        data_names <- sub("[.]pgx$","",data_names)
        upload_pgx <- sub("[.]pgx$","",r_global$load_data_from_upload)
        load_row <- which(data_names == upload_pgx)[1]
        rl$selected_row <- load_row
        loadRowManual(loadRowManual() + 1)
        r_global$load_data_from_upload <- NULL
    }, ignoreNULL = TRUE)

    observeEvent(loadRowManual(), {
        shinyjs::click('loadbutton')
    }, ignoreInit = TRUE)

    shiny::observeEvent(load_react(), {
      if (!load_react()) {
        return(NULL)
      }

      on.exit({
        bigdash.showTabsGoToDataView(session)  ## in ui-bigdashplus.R
      })

      pgxfile <- NULL

      ## Observe button press (over-rides URL query)
      btn <- input$loadbutton
      if (!is.null(btn) && btn != 0) {
        pgxfile <- selectedPGX()
      }
      ## Observe "try example dataset" press
      if (!is.null(rl$found_example_trigger)) {
        pgxfile <- selectedPGX()
      }

      ## check if file is there
      if (is.na(pgxfile) || is.null(pgxfile) || pgxfile == "" || length(pgxfile) == 0) {
        message("[LoadingBoard@load_react] ERROR file not found : ", pgxfile, "\n")
        return(NULL)
      }

      ## During loading show loading pop-up modal
      pgx.showCartoonModal()

      loadAndActivatePGX(pgxfile)

      ## notify new data uploaded
      r_global$loadedDataset <- r_global$loadedDataset + 1
      rl$found_example_trigger <- NULL
    })


    ## ---------------------------------------------------------------------
    ## ----------------- Load and activate PGX object ----------------------
    ## ---------------------------------------------------------------------
    loadAndActivatePGX <- function(pgxfile) {

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
      ## ----------------- remove modal on exit?? -------------------------
      remove(loaded_pgx)
      info("[loading_server.R] copying pgx done!")
      gc()
    }


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

    # re-write datasets-info.csv when pgxTable_edited
    # also edit the pgx files
    observeEvent(rl$pgxTable_edited, {
      pgxdir <- getPGXDIR()
      fname <- file.path(pgxdir, 'datasets-info.csv')
      pgxinfotable <- rl$pgxTable_data
      write.csv(pgxinfotable, fname)

      ## also rewrite description in actual pgx file
      pgx_name <- pgxinfotable[rl$pgxTable_edited_row, 'dataset']
      pgx_file <- file.path(pgxdir, paste0(pgx_name, '.pgx'))
      pgx <- local(get(load(pgx_file, verbose = 0))) ## override any name

      col_edited <- colnames(pgxinfotable)[rl$pgxTable_edited_col]
      new_val <- rl$pgxTable_data[rl$pgxTable_edited_row, rl$pgxTable_edited_col]

      pgx[[col_edited]] <- new_val
      playbase::pgx.save(pgx, file = pgx_file)
      remove(pgx)

    }, ignoreInit = TRUE)

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
