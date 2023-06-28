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
      delete_trigger = 0
    )

    renewReceivedTable <- reactiveVal(0)

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
    observeEvent(pgxtable$rows_selected(), {
      rl$selected_row <- pgxtable$rows_selected()
    }, ignoreNULL = FALSE)

    # disable buttons when no row is selected; enable when one is selected
    observeEvent(rl$selected_row, {
      if (is.null(rl$selected_row)) {
        shinyjs::disable(id = 'loadbutton')
      } else {
        shinyjs::enable(id = 'loadbutton')
      }
    }, ignoreNULL = FALSE)


    ##-------------------------------------------------------------------
    ## util functions
    ##-------------------------------------------------------------------

    is_valid_email <- function(email) {
        is_personal <- grepl("gmail|ymail|outlook|yahoo|hotmail|mail.com$|icloud",email)
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

    ##-------------------------------------------------------------------
    ## share dataset with specific user
    ##-------------------------------------------------------------------
    output$receive_pgx_alert <- renderUI({

        if(auth$email()=="") return(NULL)
        pgx_received <- getReceivedFiles()

        if(length(pgx_received)>0) {
            df <- receivedPGXtable()
            dt <- DT::datatable(
                df,
                rownames = FALSE,
                escape = FALSE,
                selection = "none",
                ## class = "compact cell-border",
                class = "compact row-border",
                options = list(
                    dom = "t",
                    pageLength = 999
                )
             ) %>%
                DT::formatStyle(0, target = "row", fontSize = "14px", lineHeight = "90%")

            bs_alert(
               style = "danger",
               tagList(
                   shiny::HTML("<strong>A user has shared data with you!</strong>"),
                   dt
               )
            )
        }
    })

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
      pgx_name <- rl$pgxTable_data[selected_row, "dataset"]
      pgx_name
    })

    # put user dataset into shared folder
    observeEvent( rl$share_pgx, {

        ## reset for next observe
        rl$share_pgx <- NULL
        
        ## sharing folder has to exists
        if(!dir.exists(pgx_shared_dir)) {
          shinyalert::shinyalert(
            title = "Oops! Cannot share...",
            text = paste('This server does not support sharing.',
                         'Please contact your administrator.')
          )
          return()
        }
        
        ## user has to be logged in and have email for them to share with other users
        if (auth$email() == "") {
          shinyalert::shinyalert(
             title = "Oops! You're not logged in...",
             text = paste("You need to be logged in with a valid email",
                          "address to share pgx files with other users.")
           )
          return()
        }

        ## check how many are already in queue
        pp <- paste0("__from__",auth$email(),"__$")
        num_shared_queue <- length(dir(pgx_shared_dir, pattern=pp))
        if(num_shared_queue > opt$MAX_SHARED_QUEUE) {   ## NB opt is global...
          shinyalert::shinyalert(
            title = "Oops! Too many shared...",
            text = paste("You have already too many shared datasets in the waiting queue.",
                         "Please contact your administrator.")
          )
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
            output$error_alert <- renderText({'Email is not valid. Please use
          only work or business emails.'})
            return()
        }
        if (share_user==auth$email()) {
            output$error_alert <- renderText({'Error. You cannot share with yourself...'})
            return()
        }

        output$error_alert <- renderText({''})

        shiny::removeModal()

        pgx_name <- selectedPGX2()

        alert_val <- shinyalert::shinyalert(
            inputId = "share_confirm",
            title = "Please confirm",
            tagList(
                paste("Your dataset", pgx_name, "will be shared with user",
                      share_user,". Are you sure?")
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
            path_to_creds <- 'gmail_creds'
            sender <- auth$email()
            if (file.exists(path_to_creds)) {
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
        }

        rl$share_pgx <- NULL
    })

    ##======================================================================
    ## LOAD EXAMPLE TRIGGER
    ##======================================================================
    observeEvent(r_global$load_example_trigger, {

      data_names <- as.character(pgxtable$data()$dataset)

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
    })

    ##======================================================================
    ## PUBLIC DATASET FOLDER
    ##======================================================================
    if(enable_public_tabpanel) {

      ##-------------------------------------------------------------------
      ## observers
      ##-------------------------------------------------------------------

      observeEvent(pgxtable_public$rows_selected(), {
        rl$selected_row_public <- pgxtable_public$rows_selected()
      }, ignoreNULL = FALSE)

      # disable buttons when no row is selected; enable when one is selected
      observeEvent(rl$selected_row_public, {
        if (is.null(rl$selected_row_public)) {
          shinyjs::disable(id = 'importbutton')
        } else {
          shinyjs::enable(id = 'importbutton')
        }
      }, ignoreNULL = FALSE)

      observeEvent(
        input$importbutton, {
          selected_row <- rl$selected_row_public
          pgx_name <- rl$pgxTablePublic_data[selected_row, 'dataset']
          pgx_file <- file.path(pgx_public_dir, paste0(pgx_name, '.pgx'))
          pgx_path <- getPGXDIR()
          new_pgx_file <- file.path(pgx_path, paste0(pgx_name, '.pgx'))

          if(file.exists(new_pgx_file)) {
            shinyalert::shinyalert(
              title = "Oops! File exists...",
              paste('There is already a dataset called', pgx_name,
                'in your dataset folder. Please delete your file first.')
            )
            return()
          }

          file.copy(from = pgx_file, to = new_pgx_file)
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
          pgx_name <- rl$pgxTable_data[selected_row, "dataset"]

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
          rl$share_public_pgx <- NULL
        },
        ignoreNULL = TRUE
      )

      observeEvent(input$share_public_confirm, {

        if (input$share_public_confirm) {

          selected_row <- as.numeric(stringr::str_split(rl$share_public_pgx, "_row_")[[1]][2])
          pgx_name <- rl$pgxTable_data[selected_row, "dataset"]
          pgx_name <- sub("[.]pgx$", "", pgx_name)
          pgx_path <- getPGXDIR()
          pgx_file <- file.path(pgx_path, paste0(pgx_name, '.pgx'))
          new_pgx_file <- file.path(pgx_public_dir, paste0(pgx_name, '.pgx'))
          pgx_file <- file.path(pgx_path, paste0(pgx_name, ".pgx"))

          new_pgx_file <- file.path(
            pgx_public__dir,
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

      })
    }


    ## ================================================================================
    ## Modules
    ## ================================================================================
    loading_tsne_server(
      id = "tsne",
      pgx.dir = getPGXDIR,
      info.table = getFilteredPGXINFO,
      r_selected = reactive(pgxtable$rows_all()),
      watermark = WATERMARK
    )

    pgxtable <- loading_table_datasets_server(
      id = "pgxtable",
      rl = rl,
      enable_pgxdownload = enable_pgxdownload,
      enable_delete = enable_delete,
      enable_public_share = enable_public_share,
      enable_user_share = enable_user_share
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

    observeEvent(getPGXINFO(), {
      df <- getPGXINFO()
      if(is.null(df)) return()
      datatypes <- sort(setdiff(df$datatype, c(NA, "")))
      organisms <- sort(setdiff(df$organism, c(NA, "")))
      shiny::updateCheckboxGroupInput(session, "flt_datatype", choices = datatypes)
      shiny::updateCheckboxGroupInput(session, "flt_organism", choices = organisms)
    })

    ## -----------------------------------------------------------------------------
    ## READ initial PGX file info
    ## -----------------------------------------------------------------------------

    ## Get the pgx folder. If user folders are enabled, the user email
    ## is appended to the pgx dirname.
    getPGXDIR <- shiny::reactive({
      r_global$reload_pgxdir ## force reload

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
        warning("[LoadingBoard:getPGXINFO] user not logged in!")
        return(NULL)
      }

      pgxdir <- getPGXDIR()
      info <- NULL

      shiny::withProgress(message = "Checking datasets library...", value = 0.33, {
        FOLDER_UPDATE_STATUS <- playbase::pgx.scanInfoFile(
                     pgxdir, file = "datasets-info.csv", verbose = TRUE)
      })

      if(FOLDER_UPDATE_STATUS$INITDATASETFOLDER == TRUE) {
        pgx.showSmallModal("Updating your library<br>Please wait...")
        shiny::withProgress(message = "Updating your library...", value = 0.33, {
          info <- playbase::pgx.initDatasetFolder(
            pgxdir,
            pgxinfo = FOLDER_UPDATE_STATUS$pgxinfo,
            pgx.files = FOLDER_UPDATE_STATUS$pgx.files,
            pgxinfo.changed = FOLDER_UPDATE_STATUS$pgxinfo.changed,
            pgxfc.changed = FOLDER_UPDATE_STATUS$pgxfc.changed,
            info.file1 = FOLDER_UPDATE_STATUS$info.file1,
            pgx.missing = FOLDER_UPDATE_STATUS$pgx.missing,
            pgx.missing0 = FOLDER_UPDATE_STATUS$pgx.missing0,
            pgx.missing1 = FOLDER_UPDATE_STATUS$pgx.missing1,
            update.sigdb = FALSE, ## we do it later
            verbose = TRUE
          )

          ## before reading the info file, we need to update for new files
          shiny::removeModal(session)
          return(info)
        })

      }
      info <- playbase::pgxinfo.read(pgxdir, file = "datasets-info.csv")
      shiny::removeModal(session)
      return(info)
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
          FOLDER_UPDATE_STATUS <- playbase::pgx.scanInfoFile(pgx_public_dir, file = "datasets-info.csv", verbose = TRUE)
        })

      if(FOLDER_UPDATE_STATUS$INITDATASETFOLDER == TRUE) {
        pgx.showSmallModal("Updating datasets library<br>Please wait...")
        shiny::withProgress(message = "Updating datasets library...", value = 0.33, {
          info <- playbase::pgx.initDatasetFolder(
            pgx_public_dir,
            pgxinfo = FOLDER_UPDATE_STATUS$pgxinfo,
            pgx.files = FOLDER_UPDATE_STATUS$pgx.files,
            info.file1 = FOLDER_UPDATE_STATUS$info.file1,
            pgxinfo.changed = FOLDER_UPDATE_STATUS$pgxinfo.changed,
            pgxfc.changed = FOLDER_UPDATE_STATUS$pgxfc.changed,
            pgx.missing = FOLDER_UPDATE_STATUS$pgx.missing,
            pgx.missing0 = FOLDER_UPDATE_STATUS$pgx.missing0,
            pgx.missing1 = FOLDER_UPDATE_STATUS$pgx.missing1,
            update.sigdb = FALSE,  ## we do it later
            verbose = TRUE)
          ## before reading the info file, we need to update for new files
          shiny::removeModal(session)
          return(info)
        })
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
      req(pgxtable)
      sel <- rl$selected_row
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

    getReceivedFiles <- shiny::reactive({
      req(auth)
      if (!auth$logged()) return(c())
      if(auth$email()=="") return(c())      
      ## allow trigger for when a shared pgx is accepted / decline
      renewReceivedTable()
      pgxfiles <- dir(pgx_shared_dir, pattern = paste0('__to__',auth$email(),'__from__'))
      return(pgxfiles)
    })

    makebuttonInputs2 <- function(FUN, len, id, ...) {
        inputs <- character(length(len))
        for (i in seq_along(len)) {
            inputs[i] <- as.character(FUN(paste0(id, len[i]), ...))
        }
        inputs
    }

    receivedPGXtable <- shiny::eventReactive(
        c(r_global$nav, renewReceivedTable()),
    {
        req(r_global$nav == 'load-tab')
        shared_pgx_names <- getReceivedFiles()

        accept_btns <- makebuttonInputs2(
            FUN = actionButton,
            len = shared_pgx_names,
            id = "accept_pgx__",
            label = "",
            width = "50px",
            inline = TRUE,
            icon = shiny::icon("check"),
            class = "btn-inline btn-success",
            style = "padding:0px; margin:0px; font-size:85%;",
            onclick=paste0('Shiny.onInputChange(\"',ns("accept_pgx"),'\", this.id, {priority: "event"})')
        )

        decline_btns <- makebuttonInputs2(
            FUN = actionButton,
            len = shared_pgx_names,
            id = "decline_pgx__",
            label = "",
            width = "50px",
            inline = TRUE,
            icon = shiny::icon("x"),
            class = "btn-inline btn-danger",
            style = "padding:0px; margin:0px; font-size:85%;",
            onclick=paste0('Shiny.onInputChange(\"',ns("decline_pgx"),'\", this.id, {priority: "event"})')
        )

        # split the file name into user who shared and file name
        #shared_pgx_df <- data.frame(stringr::str_split(shared_pgx_names, '__from__'))
        # remove the last '__' from the file name
        #shared_pgx_df[2,] <- stringr::str_sub(shared_pgx_df[2,], end=-3)
        shared_dataset <- sub("__to__.*","",shared_pgx_names)
        shared_from <- gsub(".*__from__|__$","",shared_pgx_names)

        df <- data.frame(
            Dataset = shared_dataset,
            From = shared_from,
            Actions = paste(accept_btns, decline_btns)
        )
        df
    })

    # event when a shared pgx is accepted by a user
    observeEvent(input$accept_pgx, {
        # get pgx name and remove the __from__* tag
        pgx_name <- stringr::str_split(input$accept_pgx, 'accept_pgx__')[[1]][2]
        new_pgx_name <- stringr::str_split(pgx_name, '__from__')[[1]][1]
        new_pgx_name <- sub("__to__.*","",pgx_name)

        # rename the file to be a valid pgx file
        pgdir <- getPGXDIR()
        file_from <- file.path(pgx_shared_dir, pgx_name)
        file_to   <- file.path(pgdir, new_pgx_name)

        if(file.exists(file_to)) {
            shinyalert::shinyalert(
                "File already exists!",
                paste("You have already a dataset called '", new_pgx_name,
                      "'. Please delete it before accepting the new file."),
                confirmButtonText = "Cancel",
                showCancelButton = FALSE
            )
            return(NULL)
        } else {
            dbg("[loading_server.R] accept_pgx : renaming file from = ",file_from,"to = ",file_to)
            if(!file.rename(file_from, file_to)) {
              info("[loading_server.R] accept_pgx : rename failed. trying file.copy ")
              ## file.rename does not allow "cross-device link"
              file.copy(file_from, file_to)
              file.remove(file_from)
            }
        }

        # reload pgx dir so the newly accepted pgx files are registered in user table
        r_global$reload_pgxdir <- r_global$reload_pgxdir + 1

        # remove the accepted pgx from the table
        renewReceivedTable(renewReceivedTable() + 1)
    }, ignoreInit = TRUE)

    # event when a shared pgx is declined by a user
    observeEvent(input$decline_pgx, {
        pgx_name <- stringr::str_split(input$decline_pgx, 'decline_pgx__')[[1]][2]
        ## pgdir <- getPGXDIR()
        shared_file <- file.path(pgx_shared_dir, pgx_name)
        dbg("[loading_server.R] decline_pgx : removing shared_file = ",shared_file)
        file.remove(shared_file)

        # remove the declined pgx from the table
        renewReceivedTable(renewReceivedTable() + 1)
    }, ignoreInit = TRUE)

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

      not.anonymous <- (!is.na(auth$name()) && auth$name() != "")
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
        data_names <- as.character(pgxtable$data()$dataset)
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

      ## Observe URL query
      ## query <- parseQueryString(session$clientData$url_search)
      ## if (!is.null(query[["pgx"]])) {
      ##   pgxfile <- query[["pgx"]]
      ##   pgxfile <- basename(pgxfile) ## for security
      ##   pgxfile <- paste0(sub("[.]pgx$", "", pgxfile), ".pgx") ## add/replace .pgx
      ## }

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

      activatePGX(pgxfile)

      ## notify new data uploaded
      r_global$loadedDataset <- r_global$loadedDataset + 1
      rl$found_example_trigger <- NULL
    })


    ## ---------------------------------------------------------------------
    ## ----------------- Load and activate PGX object ----------------------
    ## ---------------------------------------------------------------------
    activatePGX <- function(pgxfile) {

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
      pgx <- getFilteredPGXINFO()
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

    observeEvent(
      c(getFilteredPGXINFO(), r_global$reload_pgxdir), {
        df <- getFilteredPGXINFO()
        df$dataset <- sub("[.]pgx$", "", df$dataset)
        df$conditions <- gsub("[,]", " ", df$conditions)
        df$conditions <- sapply(as.character(df$conditions), andothers, split = " ", n = 5)
        df$description <- playbase::shortstring(as.character(df$description), 200)
        df$nsets <- NULL
        df$organism <- NULL

        rl$pgxTable_data <- df
      }
    )

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
      write.csv(rl$pgxTable_data, fname)

      ## also rewrite description in actual pgx file
      pgx_name <- rl$pgxTable_data[rl$pgxTable_edited_row, 'dataset']
      pgx_file <- file.path(pgxdir, paste0(pgx_name, '.pgx'))
      pgx <- local(get(load(pgx_file, verbose = 0))) ## override any name

      col_edited <- colnames(rl$pgxTable_data)[rl$pgxTable_edited_col]
      new_val <- rl$pgxTable_data[rl$pgxTable_edited_row, rl$pgxTable_edited_col]

      pgx[[col_edited]] <- new_val
      save(pgx, file = pgx_file)
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
