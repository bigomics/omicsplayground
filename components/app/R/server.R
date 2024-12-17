##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' The main application Server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
app_server <- function(input, output, session) {
  
  info("[ui.R] >>> creating SERVER")
  message("\n===========================================================")
  message("======================== SERVER ===========================")
  message("===========================================================\n")

  VERSION <- scan(file.path(OPG, "VERSION"), character())[1]

  info("[SERVER] getwd = ", normalizePath(getwd()))
  info("[SERVER] SESSION = ", session$token)

  ## default disconnected screen
  sever::sever(sever_disconnected(), bg_color = "#004c7d")

  setwd(WORKDIR) ## for some reason it can change!! (defined in global.R)
  server.start_time <- Sys.time()
  session.start_time <- -1

  ## show warning for mobile
  is_mobile <- reactive({
    bw <- shinybrowser::get_width()
    bh <- shinybrowser::get_height()
    dbg("[SERVER] shinybrowser: ", bw, "x", bh)
    is_portrait <- bh > 1.2 * bw
    is_small <- min(bw, bh) < 768
    (is_portrait && is_small)
  })

  observeEvent(is_mobile(), {
    if (is_mobile()) {
      shinyalert::shinyalert(
        title = "Sorry, not for mobile...",
        text = "Omics Playground is not yet optimized for mobile. For the best experience, please use it from your desktop PC",
        immediate = TRUE,
        showCancelButton = FALSE,
        showConfirmButton = TRUE
      )
    }
  })

  ## -------------------------------------------------------------
  ## Authentication
  ## -------------------------------------------------------------

  authentication <- opt$AUTHENTICATION
  auth <- NULL ## shared in module
  credentials_file <- file.path(ETC, "CREDENTIALS")
  has.credentials <- file.exists(credentials_file)
  no.credentials <- (!isTRUE(opt$USE_CREDENTIALS) || !has.credentials)
  if (no.credentials && authentication != "password") {
    credentials_file <- NULL
  }

  user_database <- file.path(ETC, "user_details.sqlite")
  has.user_database <- file.exists(user_database)
  if (!has.user_database && authentication != "login-code") {
    user_database <- NULL
  }

  if (authentication == "password") {
    auth <- PasswordAuthenticationModule(
      id = "auth",
      credentials_file = credentials_file,
      allow_personal = opt$ALLOW_PERSONAL_EMAIL,
      domain = opt$DOMAIN,
      blocked_domain = opt$BLOCKED_DOMAIN
    )
    ## } else if (authentication == "firebase") {
    ##   auth <- FirebaseAuthenticationModule(
    ##     id = "auth",
    ##     domain = opt$DOMAIN,
    ##     firebase.rds = "firebase.rds",
    ##     credentials_file = credentials_file,
    ##     allow_personal = opt$ALLOW_PERSONAL_EMAIL,
    ##     allow_new_users = opt$ALLOW_NEW_USERS
    ##   )
    ## } else if (authentication == "email-link") {
    ##   auth <- EmailLinkAuthenticationModule(
    ##     id = "auth",
    ##     pgx_dir = PGX.DIR,
    ##     domain = opt$DOMAIN,
    ##     firebase.rds = "firebase.rds",
    ##     credentials_file = credentials_file,
    ##     allow_personal = opt$ALLOW_PERSONAL_EMAIL,
    ##     allow_new_users = opt$ALLOW_NEW_USERS
    ##   )
  } else if (authentication == "login-code") {
    auth <- LoginCodeAuthenticationModule(
      id = "auth",
      mail_creds = file.path(ETC, "gmail_creds"),
      domain = opt$DOMAIN,
      user_database = user_database,
      blocked_domain = opt$BLOCKED_DOMAIN,
      allow_personal = opt$ALLOW_PERSONAL_EMAIL,
      allow_new_users = opt$ALLOW_NEW_USERS,
      redirect_login = FALSE
    )
  } else if (authentication == "login-code-redirect") {
    auth <- LoginCodeAuthenticationModule(
      id = "auth",
      mail_creds = file.path(ETC, "gmail_creds"),
      domain = opt$DOMAIN,
      blocked_domain = opt$BLOCKED_DOMAIN,
      user_database = user_database,
      allow_personal = opt$ALLOW_PERSONAL_EMAIL,
      allow_new_users = opt$ALLOW_NEW_USERS,
      redirect_login = TRUE
    )
  } else if (authentication == "shinyproxy") {
    username <- Sys.getenv("SHINYPROXY_USERNAME")
    auth <- NoAuthenticationModule(
      id = "auth",
      show_modal = TRUE,
      username = username,
      email = username
    )
  } else if (authentication == "none") {
    auth <- NoAuthenticationModule(
      id = "auth",
      show_modal = TRUE
    )
  } else if (authentication == "none2") {
    ## no authentication but also not showing main modal (enter)
    username <- Sys.getenv("PLAYGROUND_USERNAME")
    auth <- NoAuthenticationModule(
      id = "auth",
      show_modal = FALSE,
      username = username,
      email = username
    )
  } else if (authentication == "apache-cookie") {
    auth <- AuthenticationModuleApacheCookie(
      id = "auth",
      show_modal = FALSE
    )
  } else {
    ## stop everything
    stop("unsupported authorization method", authentication)
  }

  ## -------------------------------------------------------------
  ## Call modules
  ## -------------------------------------------------------------

  env <- list() ## communication "environment"

  ## Global reactive value for PGX object
  PGX <- reactiveValues()

  ## Global reactive values for app-wide triggering
  load_example <- reactiveVal(NULL)
  load_uploaded_data <- reactiveVal(NULL)
  labeltype <- reactiveVal("feature") # can be feature (rownames counts), symbol or name
  reload_pgxdir <- reactiveVal(0)
  inactivityCounter <- reactiveVal(0)
  new_upload <- reactiveVal(0)

  ## Default boards ------------------------------------------
  WelcomeBoard("welcome",
    auth = auth,
    load_example = load_example,
    new_upload = new_upload
  )

  env$user_profile <- UserProfileBoard(
    "user_profile",
    auth = auth,
    nav_count = reactive(nav$count)
  )

  AppSettingsBoard(
    "app_settings",
    auth = auth,
    pgx = PGX
  )

  env$user_settings <- list(
    enable_beta = shiny::reactive(input$enable_beta),
    enable_info = shiny::reactive(input$enable_info)
  )

  ## Do not display "Welcome" tab on the menu
  bigdash.hideMenuItem(session, "welcome-tab")
  shinyjs::runjs("sidebarClose()")

  ## Modules needed from the start
  recompute_pgx <- shiny::reactiveVal(NULL)

  env$load <- LoadingBoard(
    id = "load",
    pgx = PGX,
    auth = auth,
    pgx_topdir = PGX.DIR,
    load_example = load_example,
    reload_pgxdir = reload_pgxdir,
    current_page = reactive(input$nav),
    load_uploaded_data = load_uploaded_data,
    recompute_pgx = recompute_pgx,
    new_upload = new_upload
  )

  ## Modules needed from the start
  if (opt$ENABLE_UPLOAD) {
    upload_datatype <- UploadBoard(
      id = "upload",
      pgx_dir = PGX.DIR,
      pgx = PGX,
      auth = auth,
      reload_pgxdir = reload_pgxdir,
      load_uploaded_data = load_uploaded_data,
      recompute_pgx = recompute_pgx,
      inactivityCounter = inactivityCounter,
      new_upload = new_upload
    )


    shiny::observeEvent(upload_datatype(), {
      if (grepl("proteomics", upload_datatype(), ignore.case = TRUE)) {
        shiny.i18n::update_lang("proteomics", session)
      } else if (tolower(upload_datatype()) == "metabolomics") {
        shiny.i18n::update_lang("metabolomics", session)
      } else {
        shiny.i18n::update_lang("RNA-seq", session)
      }
    })
  }

  ## Modules needed after dataset is loaded (deferred) --------------
  observeEvent(env$load$is_data_loaded(), {
    if (env$load$is_data_loaded() == 1) {
      
      additional_ui_tabs <- list(
        dataview = bigdash::bigTabItem(
          "dataview-tab",
          DataViewInputs("dataview"),
          DataViewUI("dataview")
        ),
        diffexpr = bigdash::bigTabItem(
          "diffexpr-tab",
          ExpressionInputs("diffexpr"),
          ExpressionUI("diffexpr")
        ),
        enrich = bigdash::bigTabItem(
          "enrich-tab",
          EnrichmentInputs("enrich"),
          EnrichmentUI("enrich")
        ),
        clustersamples = bigdash::bigTabItem(
          "clustersamples-tab",
          ClusteringInputs("clustersamples"),
          ClusteringUI("clustersamples")
        ),
        clusterfeatures = bigdash::bigTabItem(
          "clusterfeatures-tab",
          FeatureMapInputs("clusterfeatures"),
          FeatureMapUI("clusterfeatures")
        ),
        wgcna = bigdash::bigTabItem(
          "wgcna-tab",
          WgcnaInputs("wgcna"),
          WgcnaUI("wgcna")
        ),
        pcsf = bigdash::bigTabItem(
          "pcsf-tab",
          PcsfInputs("pcsf"),
          PcsfUI("pcsf")
        ),
        corr = bigdash::bigTabItem(
          "corr-tab",
          CorrelationInputs("corr"),
          CorrelationUI("corr")
        ),
        pathway = bigdash::bigTabItem(
          "pathway-tab",
          PathwayInputs("pathway"),
          PathwayUI("pathway")
        ),
        wordcloud = bigdash::bigTabItem(
          "wordcloud-tab",
          WordCloudInputs("wordcloud"),
          WordCloudUI("wordcloud")
        ),
        drug = bigdash::bigTabItem(
          "drug-tab",
          DrugConnectivityInputs("drug"),
          DrugConnectivityUI("drug")
        ),
        isect = bigdash::bigTabItem(
          "isect-tab",
          IntersectionInputs("isect"),
          IntersectionUI("isect")
        ),
        sig = bigdash::bigTabItem(
          "sig-tab",
          SignatureInputs("sig"),
          SignatureUI("sig")
        ),
        bio = bigdash::bigTabItem(
          "bio-tab",
          BiomarkerInputs("bio"),
          BiomarkerUI("bio")
        ),
        cmap = bigdash::bigTabItem(
          "cmap-tab",
          ConnectivityInputs("cmap"),
          ConnectivityUI("cmap")
        ),
        comp = bigdash::bigTabItem(
          "comp-tab",
          CompareInputs("comp"),
          CompareUI("comp")
        ),
        cell = bigdash::bigTabItem(
          "cell-tab",
          SingleCellInputs("cell"),
          SingleCellUI("cell")
        ),
        tcga = bigdash::bigTabItem(
          "tcga-tab",
          TcgaInputs("tcga"),
          TcgaUI("tcga")
        )
        ## mofa = bigdash::bigTabItem(
        ##   "mofa-tab",
        ##   MofaInputs("mofa"),
        ##   MofaUI("mofa")
        ## ),        
        ## mgsea = bigdash::bigTabItem(
        ##   "mgsea-tab",
        ##   MGseaInputs("mgsea"),
        ##   MGseaUI("mgsea")
        ## ),        
        ## snf = bigdash::bigTabItem(
        ##   "snf-tab",
        ##   SNF_Inputs("snf"),
        ##   SNF_UI("snf")
        ## )        
      )

      insertBigTabUI <- function(ui) {
        ##if( inherits(class(ui),"list")) ui <- list(ui)
        for(i in 1:length(ui)) {
          shiny::insertUI(
            selector = "#big-tabs",
            where = "beforeEnd",
            ui = ui[[i]],
            immediate = TRUE
          )
        }
      }
      insertBigTabItem <- function(tab) {
        insertBigTabUI( additional_ui_tabs[tab] ) 
      }

      shiny::withProgress(
        message = "Preparing your dashboard (server)...",
        value = 0,
        {
          if (MODULES_ENABLED["DataView"]) {
            info("[SERVER] calling DataView module")
            insertBigTabItem("dataview")
            DataViewBoard("dataview", pgx = PGX, labeltype = labeltype)
          }
          shiny::incProgress(0.1)

          if(MODULES_ENABLED['Clustering']) {
            info("[SERVER] calling ClusteringBoard module")
            insertBigTabItem("clustersamples")
            ClusteringBoard("clustersamples", pgx = PGX, labeltype = labeltype)

            info("[SERVER] calling FeatureMapBoard module")
            insertBigTabItem("clusterfeatures")
            FeatureMapBoard("clusterfeatures", pgx = PGX, labeltype = labeltype)
          }
          shiny::incProgress(0.1)
          
          if(MODULES_ENABLED['Expression']) {
            
            info("[SERVER] calling ExpressionBoard module")
            insertBigTabItem("diffexpr")
            ExpressionBoard("diffexpr", pgx = PGX, labeltype = labeltype) -> env$diffexpr

            info("[SERVER] calling CorrelationBoard module")
            insertBigTabItem("corr")
            CorrelationBoard("corr", pgx = PGX, labeltype = labeltype)
            
            info("[SERVER] calling BiomarkerBoard module")
            insertBigTabItem("bio")
            BiomarkerBoard("bio", pgx = PGX)

          }
          shiny::incProgress(0.1)
          
          if(MODULES_ENABLED['GeneSets']) {
            
            info("[SERVER] calling EnrichmentBoard module")
            insertBigTabItem("enrich")
            EnrichmentBoard("enrich",
              pgx = PGX,
              selected_gxmethods = env$diffexpr$selected_gxmethods
            ) -> env$enrich
            
            info("[SERVER] calling SignatureBoard module")
            insertBigTabItem("sig")
            SignatureBoard("sig",
              pgx = PGX,
              selected_gxmethods = env$diffexpr$selected_gxmethods
            )

            info("[SERVER] calling PathwayBoard module")
            insertBigTabItem("pathway")
            PathwayBoard("pathway",
              pgx = PGX,
              selected_gsetmethods = env$enrich$selected_gsetmethods
            )
          
            info("[SERVER] calling WordCloudBoard module")
            insertBigTabItem("wordcloud")
            WordCloudBoard("wordcloud", pgx = PGX)
          }
          shiny::incProgress(0.1)

          if(MODULES_ENABLED['Compare']) {

            info("[SERVER] calling IntersectionBoard module")
            insertBigTabItem("isect")
            IntersectionBoard("isect",
              pgx = PGX,
              selected_gxmethods = env$diffexpr$selected_gxmethods,
              selected_gsetmethods = env$enrich$selected_gsetmethods
            )

            info("[SERVER] calling CompareBoard module")
            insertBigTabItem("comp")
            CompareBoard("comp", pgx = PGX, pgx_dir = reactive(auth$user_dir),
                         labeltype = labeltype)

            info("[SERVER] calling ConnectivityBoard module")
            insertBigTabItem("cmap")
            ConnectivityBoard("cmap",
              pgx = PGX,
              auth = auth,
              reload_pgxdir = reload_pgxdir
            )
          }
          shiny::incProgress(0.1)

          if (MODULES_ENABLED["SystemsBio"]) {
            
            info("[SERVER] calling DrugConnectivityBoard module")
            insertBigTabItem("drug")
            DrugConnectivityBoard("drug", pgx = PGX)
            
            info("[SERVER] calling SingleCellBoard module")
            insertBigTabItem("cell")
            SingleCellBoard("cell", pgx = PGX)

            info("[SERVER] calling PcsfBoard module")
            insertBigTabItem("pcsf")
            PcsfBoard("pcsf", pgx = PGX)

            info("[SERVER] calling WgcnaBoard module")
            insertBigTabItem("wgcna")
            WgcnaBoard("wgcna", pgx = PGX)

            info("[SERVER] calling TcgaBoard module")
            insertBigTabItem("tcga")
            TcgaBoard("tcga", pgx = PGX)

          }
          shiny::incProgress(0.1)
          
          if (MODULES_ENABLED["MultiOmics"] && exists("MODULE.multiomics")) {
            info("[SERVER] initializing MultiOmics module")
            mod <- MODULE.multiomics
            insertBigTabUI( mod$module_ui() ) 
            mod$module_server(PGX)
          }
          
          info("[SERVER] calling modules done!")
        }
      )
    }

    if (env$load$is_data_loaded() == 1) {
      # this is a function - like "handleSettings()" in bigdash- needed to
      # make the settings sidebar show up for the inserted tabs
      shinyjs::runjs(
        "  $('.big-tab')
    .each((index, el) => {
      let settings = $(el)
        .find('.tab-settings')
        .first();

      $(settings).data('target', $(el).data('name'));
      $(settings).appendTo('#settings-content');
    });"
      )
    }

    bigdash.selectTab(session, selected = "dataview-tab")
    bigdash.openSettings()

    ## remove loading modal from LoadingBoard
    shiny::removeModal()

    bigdash.showTabsGoToDataView(session)
  })


  ## --------------------------------------------------------------------------
  ## Current navigation
  ## --------------------------------------------------------------------------

  shinyjs::onclick("logo-bigomics", {
    shinyjs::runjs("console.info('logo-bigomics clicked')")
    bigdash.selectTab(session, selected = "welcome-tab")
    shinyjs::runjs("sidebarClose()")
    shinyjs::runjs("settingsClose()")
  })

  observeEvent(input$menu_createreport, {
    shinyjs::click("load-generate_report-show_report_modal")
  })

  output$current_user <- shiny::renderText({
    ## trigger on change of user
    shiny::req(auth$logged)
    if (is.null(auth$logged) || !auth$logged) {
      return("(not logged in)")
    }

    user <- auth$email
    if (is.null(user) || user %in% c("", NA)) user <- auth$username
    if (is.null(user) || user %in% c("", NA)) user <- "User"
    user
  })

  ## this avoids flickering
  welcome_detector <- reactiveVal(NULL)
  observeEvent(input$nav, {
    welcome_detector(input$nav == "welcome-tab")
  })

  output$current_dataset <- shiny::renderUI({
    shiny::req(auth$logged)
    has.pgx <- !is.null(PGX$name) && length(PGX$name) > 0
    nav.welcome <- welcome_detector()
    if (isTRUE(auth$logged) && has.pgx && !nav.welcome) {
      ## trigger on change of dataset
      pgx.name <- gsub(".*\\/|[.]pgx$", "", PGX$name)
      tag <- shiny::actionButton(
        "dataset_click", pgx.name,
        class = "quick-button",
        style = "border: none; color: black; font-size: 0.9em;"
      )
    } else {
      tag <- HTML(paste("Omics Playground", VERSION))
    }
    tag
  })

  observeEvent(input$dataset_click, {
    shiny::req(PGX$name)
    pgx.name <- gsub(".*\\/|[.]pgx$", "", PGX$name)
    ##    fields <- c("name", "datatype", "description", "date", "norm_method", "imputation_method", "bc_method", "remove_outliers")
    fields <- c("name", "datatype", "description", "date", "settings")
    fields <- intersect(fields, names(PGX))
    body <- ""
    listcollapse <- function(lst) paste0(names(lst), "=", lst, collapse = "; ")
    for (f in fields) {
      if (length(PGX[[f]]) > 1) {
        for (n in names(PGX[[f]])) {
          val <- PGX[[f]][[n]]
          if (length(val) > 1) val <- listcollapse(val)
          txt1 <- paste0("<b>", paste0(f, ".", n), ":</b>&nbsp; ", val, "<br>", collapse = "")
          body <- paste(body, txt1)
        }
      } else {
        txt1 <- paste0("<b>", f, ":</b>&nbsp; ", PGX[[f]], "<br>")
        body <- paste(body, txt1)
      }
    }

    shiny::showModal(shiny::modalDialog(
      title = pgx.name,
      div(HTML(body), style = "font-size: 1.1em;"),
      footer = NULL,
      size = "l",
      easyClose = TRUE,
      fade = FALSE
    ))
  })


  ## count the number of times a navtab is clicked during the session
  nav <- reactiveValues(count = c())
  observeEvent(input$nav, {
    i <- input$nav
    if (is.null(nav$count[i]) || is.na(nav$count[i])) {
      nav$count[i] <<- 1
    } else {
      nav$count[i] <<- nav$count[i] + 1
    }
  })

  ## --------------------------------------------------------------------------
  ## upon change of user OR beta toggle OR new pgx
  ## --------------------------------------------------------------------------

  shiny::observeEvent(
    {
      auth$logged
      env$user_settings$enable_beta()
      PGX$name
    },
    {
      ## trigger on change dataset
      dbg("[SERVER] trigger on new PGX")

      ## write GLOBAL variables
      LOADEDPGX <<- PGX$name
      DATATYPEPGX <<- tolower(PGX$datatype)

      ## change language
      if (grepl("proteomics", DATATYPEPGX, ignore.case = TRUE)) {
        lang <- "proteomics"
      } else if (DATATYPEPGX == "metabolomics") {
        lang <- "metabolomics"
      } else {
        lang <- "RNA-seq"
      }
      dbg("[SERVER] changing 'language' to", lang)
      shiny.i18n::update_lang(lang, session)

      # choose the default labeltype based on datatype

      if (PGX$datatype == "metabolomics") {
        labeltype("gene_title")
      } else {
        labeltype("feature") # probe is feature (rownames of counts)
      }

      ## show beta feauture
      show.beta <- env$user_settings$enable_beta()
      if (is.null(show.beta) || length(show.beta) == 0) show.beta <- FALSE
      is.logged <- auth$logged

      ## hide all main tabs until we have an object
      if (is.null(PGX) || is.null(PGX$name) || !is.logged) {
        warning("[SERVER] !!! no data. hiding menu.")
        shinyjs::runjs("sidebarClose()")
        shinyjs::runjs("settingsClose()")
        bigdash.selectTab(session, selected = "welcome-tab")
        return(NULL)
      }

      ## show all main tabs
      shinyjs::runjs("sidebarOpen()")
      shinyjs::runjs("settingsOpen()")

      ## do we have libx libraries?
      has.libx <- dir.exists(file.path(OPG, "libx"))

      ## Beta features
      info("[SERVER] disabling beta features")
      bigdash.toggleTab(session, "tcga-tab", show.beta && has.libx)
      toggleTab("drug-tabs", "Connectivity map (beta)", show.beta) ## too slow
      toggleTab("pathway-tabs", "Enrichment Map (beta)", show.beta) ## too slow
      
      
      ## Dynamically show upon availability in pgx object
      info("[SERVER] disabling extra features")
      tabRequire(PGX, session, "wgcna-tab", "wgcna", TRUE)
      tabRequire(PGX, session, "drug-tab", "drugs", TRUE)
      tabRequire(PGX, session, "wordcloud-tab", "wordcloud", TRUE)
      tabRequire(PGX, session, "cell-tab", "deconv", TRUE)
      gset_tabs <- c("enrich-tab", "pathway-tab", "isect-tab", "sig-tab")
      for (tab_i in gset_tabs) {
        tabRequire(PGX, session, tab_i, "gsetX", TRUE)
        tabRequire(PGX, session, tab_i, "gset.meta", TRUE)
      }

      ## Hide PCSF and WGCNA for metabolomics
      # WGCNA will be abailable upon gmt refactoring
      if (DATATYPEPGX == "metabolomics") {
        info("[SERVER] disabling WGCNA and PCSF for metabolomics data")
        bigdash.hideTab(session, "pcsf-tab")
        bigdash.hideTab(session, "wgcna-tab")
        bigdash.hideTab(session, "cmap-tab")
      }

      info("[SERVER] trigger on change dataset done!")
    }
  )

  # populate labeltype selector based on pgx$genes

  observeEvent(
    {
      PGX$genes
    },
    {
      req(PGX$genes)

      clean_genes_matrix <- PGX$genes

      # remove NA columns
      clean_genes_matrix <- clean_genes_matrix[, !apply(is.na(clean_genes_matrix), 2, all), drop = FALSE]

      # remove columns with only 1 unique value
      clean_genes_matrix <- clean_genes_matrix[, sapply(clean_genes_matrix, function(x) length(unique(x)) > 1), drop = FALSE]

      # remove duplicated columns
      clean_genes_matrix <- clean_genes_matrix[, !duplicated(t(clean_genes_matrix)), drop = FALSE]


      # improve naming of label types (gene_title -> name) and remove pos, map, tx_len
      label_types_available <- colnames(clean_genes_matrix)

      names(label_types_available) <- label_types_available

      # rename gene_title name to name
      names(label_types_available)[names(label_types_available) == "gene_title"] <- "name"

      # remove pos, map, tx_len, source (not interesting for label types)
      label_types_available <- label_types_available[!grepl("pos|map|tx_len|source", names(label_types_available))]

      # if available, rename chr0 to Chromossome and chr to locus
      if ("chr0" %in% names(label_types_available)) {
        names(label_types_available)[names(label_types_available) == "chr0"] <- "Chromossome"
      }

      if ("chr" %in% names(label_types_available)) {
        names(label_types_available)[names(label_types_available) == "chr"] <- "Locus"
      }

      shiny::updateSelectInput(
        session,
        "selected_labeltype",
        choices = label_types_available,
        selected = labeltype()
      )
    }
  )



  # change label type based on selected input
  shiny::observeEvent(
    {
      input$selected_labeltype
    },
    {
      labeltype(input$selected_labeltype)
    }
  )

  # if labeltype is updated (via different pgx data types), we need to update the selector choice
  shiny::observeEvent(
    {
      labeltype()
    },
    {
      # if input$selected_labeltype does not match labeltype, we need to update the selector
      if (input$selected_labeltype != labeltype()) {
        shiny::updateSelectInput(
          session,
          "selected_labeltype",
          selected = labeltype()
        )
      }
    }
  )

  ## -------------------------------------------------------------
  ## Session Timers
  ## -------------------------------------------------------------
  ## invite module (from menu)
  invite <- InviteFriendModule(
    id = "invite",
    auth = auth,
    callbackR = inviteCallback
  )
  inviteCallback <- function() {
    ## After succesful invite, we extend the session
    dbg("[MAIN] inviteCB called!")
    if (isTRUE(opt$TIMEOUT > 0)) {
      session_timer$reset()
      shinyalert::shinyalert(
        text = "Thanks! We have invited your friend and your session has been extended.",
        timer = 4000
      )
    } else {
      shinyalert::shinyalert(
        text = "Thanks! We have invited your friend.",
        timer = 4000
      )
    }
  }

  session_timer <- NULL
  if (isTRUE(opt$TIMEOUT > 0)) {
    #' Session timer. Closes session after TIMEOUT (seconds) This
    #' is acitve for free users. Set TIMEOUT=0 to disable session
    #' timer.
    session_timer <- TimerModule(
      "session_timer",
      condition = reactive(auth$logged),
      timeout = opt$TIMEOUT,
      warn_before = round(0.2 * opt$TIMEOUT),
      max_warn = 1,
      warn_callback = warn_timeout,
      timeout_callback = session_timeout
    )

    warn_timeout <- function() {
      if (!auth$logged) {
        return(NULL)
      }
      shinyalert::shinyalert(
        title = "Warning!",
        text = "Your FREE session is expiring soon",
        html = TRUE,
        immediate = TRUE,
        timer = 3000
      )
    }

    ## At the end of the timeout the user can choose type of referral
    ## modal and gain additional analysis time. We reset the timer.
    session_timeout <- function() {
      if (!auth$logged) {
        return(NULL)
      }
      shinyalert::shinyalert(
        title = "FREE session expired!",
        text = "Sorry. Your free session has expired. To extend your session you can refer Omics Playground to a friend. Do you want to logout or invite a friend?",
        html = TRUE,
        immediate = TRUE,
        timer = 60 * 1000,
        showCancelButton = TRUE,
        showConfirmButton = TRUE,
        cancelButtonText = "Logout",
        confirmButtonText = "Invite friend",
        confirmButtonCol = "#AEDEF4",
        callbackR = timeout_choice
      )
    }

    ## This handles the timeout choice
    timeout_choice <- function(x) {
      dbg("[MAIN:timeout_response] x = ", x)
      if (x == FALSE) {
        ## run logout sequence
        userLogoutSequence(auth, action = "user.timeout")
        sever::sever(sever_ciao(), bg_color = "#004c7d")
        session$close()
        ## session$reload()
      }
      if (x == TRUE) {
        invite$click()
      }
    }
  } ## end of if TIMEOUT>0


  #' Idle timer. Closes session if no one is logged in after a certain
  #' period. This frees up the R process from users that are uselessly
  #' waiting at the login prompt
  idle_timer <- TimerModule(
    "idle_timer",
    condition = reactive(!auth$logged),
    timeout = 600, ## max idle time in seconds
    timeout_callback = idle_timeout_callback
  )

  idle_timeout_callback <- function() {
    info("[SERVER] ********** closing idle login session **************")
    sever::sever(sever_ciao("Knock, knock — Anybody there?"), bg_color = "#004c7d")
    session$close()
  }

  ## -------------------------------------------------------------
  ## User locking + login/logout access logging
  ## -------------------------------------------------------------
  info("[SERVER] ENABLE_USER_LOCK = ", opt$ENABLE_USER_LOCK)
  lock <- NULL
  if (isTRUE(opt$ENABLE_USER_LOCK)) {
    # Initialize the reactiveTimer to update every 30 seconds. Set max
    # idle time to 2 minutes.
    lock <- FolderLock$new(
      poll_secs = 30,
      max_idle = 120,
      show_success = FALSE,
      show_details = FALSE
    )
    lock$start_shiny_observer(auth, session = session)
  }

  #' Track which users are online by repeatedly writing every delta
  # seconds a small ID file ' in the ONLINE_DIR folder.
  if (isTRUE(opt$ENABLE_HEARTBEAT)) {
    ONLINE_DIR <- file.path(ETC, "online")
    heartbeat <- pgx.start_heartbeat(auth, session, delta = 300, online_dir = ONLINE_DIR)
    observe({
      heartbeat()
    }) ## run indefinitely
  }

  ## -------------------------------------------------------------
  ## About
  ## -------------------------------------------------------------

  observeEvent(input$navbar_about, {
    authors <- c(
      "Ana Nufer, Antonino Zito, Axel Martinelli, Carson Sievert, Cédric Scherer, Gabriela Scorici, Griffin Seidel, Ivo Kwee, John Coene, Layal Abo Khayal, Marco Sciaini, Matt Leech, Mauro Miguel Masiero, Murat Akhmedov, Nick Cullen, Santiago Caño Muñiz, Shalini Pandurangan, Stefan Reifenberg, Xavier Escribà Montagut"
    )
    authors <- paste(sort(authors), collapse = ", ")

    shiny::showModal(
      shiny::modalDialog(
        div(
          h2("BigOmics Playground"),
          h5(VERSION),
          h5("Advanced omics analysis for everyone"), br(), br(),
          p("Created with love and proudly presented to you by BigOmics Analytics from Ticino, the sunny side of Switzerland."),
          p(tags$a(href = "https://www.bigomics.ch", "www.bigomics.ch")),
          style = "text-align:center; line-height: 1em;"
        ),
        footer = div(
          "© 2000-2024 BigOmics Analytics, Inc.",
          br(), br(),
          paste("Credits:", authors),
          style = "font-size: 0.8em; line-height: 0.9em; text-align:center;"
        ),
        size = "m",
        easyClose = TRUE,
        fade = FALSE
      )
    )
  })

  ## -------------------------------------------------------------
  ## Session login sequence
  ## -------------------------------------------------------------

  onSessionStart <- isolate({
    message("*********************************************************")
    message(paste("***** new session ", session$token, "*****"))
    message("*********************************************************")
    ## users$count = users$count + 1
    ACTIVE_SESSIONS <<- c(ACTIVE_SESSIONS, session$token)
    dbg("[SERVER] number of active sessions = ", length(ACTIVE_SESSIONS))

    if (length(ACTIVE_SESSIONS) > MAX_SESSIONS) {
      dbg("ERROR: Too many sessions. stopping session!!!\n")
      pgx.record_access(
        user = isolate(auth$email),
        action = "server.full",
        comment = "too many sessions. server at capacity",
        session = session
      )
      sever::sever(sever_serverfull(opt$HOSTNAME), bg_color = "#004c7d") ## lightblue=2780e3
      session$close()
    }
  })

  ## upon change of user
  observeEvent(auth$logged, {
    if (auth$logged) {
      message("--------- user login ----------")
      message("username       = ", auth$username)
      message("email          = ", auth$email)
      message("level          = ", auth$level)
      message("limit          = ", auth$limit)
      message("user_dir       = ", auth$user_dir)
      message("-------------------------------")

      pgx.record_access(auth$email, action = "login", session = session)
    } else {
      ## clear PGX data as soon as the user logs out
      clearPGX()
    }
  })


  ## -------------------------------------------------------------
  ## Session logout sequences
  ## -------------------------------------------------------------

  clearPGX <- function() {
    ## clear PGX data
    pgx.names <- isolate(names(PGX))
    length.pgx <- length(pgx.names)
    if (length.pgx > 0) {
      for (i in 1:length.pgx) {
        PGX[[pgx.names[i]]] <<- NULL
      }
    }
  }

  userLogoutSequence <- function(auth, action) {
    message("--------- user logout ---------")
    message("username       = ", isolate(auth$username))
    message("email          = ", isolate(auth$email))
    message("user_dir       = ", isolate(auth$user_dir))
    message("------------------------------")

    ## stop all timers
    if (!is.null(session_timer)) session_timer$run(FALSE)

    ## This removes user heartbeat and lock files
    if (isTRUE(opt$ENABLE_HEARTBEAT)) {
      hbfile <- isolate(heartbeat()) ## does not work because auth has been reset
      if (file.exists(hbfile)) file.remove(hbfile)
    }
    if (!is.null(lock)) {
      lock$remove_lock()
    }

    ## record tab navigation count and time
    nav.count <- isolate(nav$count)
    nav_count.str <- paste(paste0(names(nav.count), "=", nav.count), collapse = ";")
    nav_count.str <- gsub("-tab", "", nav_count.str)

    ## record num datasets
    pgxdir <- shiny::isolate(auth$user_dir)
    num_pgxfiles <- length(dir(pgxdir, pattern = ".pgx$"))

    pgx.record_access(
      user = isolate(auth$email),
      action = action,
      session = session,
      comment = nav_count.str,
      comment2 = isolate(PLOT_DOWNLOAD_LOGGER$str),
      num_datasets = num_pgxfiles,
      ip = session$request$HTTP_X_REAL_IP
    )

    ## reset (logout) user. This should already have been done with
    ## the JS call but this is a cleaner (preferred) shiny method.
    dbg("[SERVER:userLogout] >>> resetting USER")
    isolate(auth$resetUSER())

    ## clear PGX data as soon as the user logs out (if not done)
    clearPGX()
  }


  ## This will be called upon user logout *after* the logout() JS call
  observeEvent(input$userLogout, {
    dbg("[SERVER:userLogout] user logout triggered!")

    ## run logout sequence
    userLogoutSequence(auth, action = "user.logout")

    ## this triggers a fresh session. good for resetting all
    ## parameters.
    ## (IK 16-07-2023: some bug for firebase-based login, reload
    ## loop. To be fixed!!!
    ## dbg("[SERVER:userLogout] >>> reloading session")
    ## session$reload()
  })

  ## This code listens to the JS quit signal
  observeEvent(input$quit, {
    ## Choose between reloading or closing the session. Close is
    ## better because then the user does not allocate an idle session
    ## after exit.
    dbg("[SERVER:quit] exit session... ")
    sever::sever(sever_ciao(), bg_color = "#004c7d")
    ## session$close()
    session$reload()
  })

  ## This code will be run after the client has disconnected
  ## Note!!!: Strange behaviour, sudden session ending.
  session$onSessionEnded(
    function() {
      message("******** doing session cleanup ********")

      ## Remove from active sessions
      s <- session$token
      dbg("[SERVER] removing from active sessions :", s)
      ACTIVE_SESSIONS <<- setdiff(ACTIVE_SESSIONS, s)

      ## run log sequence
      userLogoutSequence(isolate(auth), action = "session.logout")

      ## we do extra logout actions for shinyproxy
      if (opt$AUTHENTICATION == "shinyproxy") {
        session$sendCustomMessage("shinyproxy-logout", list())
      }
    }
  )

  ## -------------------------------------------------------------
  ## UI observers
  ## -------------------------------------------------------------

  # This code will run when there is a shiny error. Then this
  # error will be shown on the app. Note that errors that are
  # not related to Shiny are not caught (e.g. an error on the
  # global.R file is not caught by this)
  options(shiny.error = function() {
    # The error message is on the parent environment, it is
    # not passed to the function called on error
    parent_env <- parent.frame()
    error <- parent_env$e
    err_traceback <- NULL

    if (!is.null(error)) {
      err_traceback <- capture.output(
        printStackTrace(
          error,
          full = get_devmode_option(
            "shiny.fullstacktrace",
            FALSE
          ),
          offset = getOption("shiny.stacktraceoffset", TRUE)
        ),
        type = "message"
      )
      err_traceback <- append(error$message, err_traceback)
    }

    # Errors to ignore
    if (error$message %in% c("figure margins too large", "invalid graphics state")) {
      return()
    }
    # Get inputs to reproduce state
    board_inputs <- names(input)[grep(substr(input$nav, 1, nchar(input$nav) - 4), names(input))]

    # Remove pdf + download + card_selector + copy_info + unnecessary table inputs
    board_inputs <- board_inputs[-grep("pdf_width|pdf_height|pdf_settings|downloadOption|card_selector|copy_info|_rows_current|_rows_all", board_inputs)]

    input_values <- lapply(board_inputs, function(x) {
      value <- input[[x]]
      return(paste0(x, ": ", value))
    }) |> unlist()

    err_traceback <- append(input_values, err_traceback)

    # clean up and concatenate err_traceback

    err_traceback <- paste(err_traceback, collapse = "\n ")

    # Skip if error is repeated
    # Initialize variable as well
    if (!exists("err_prev")) {
      err_prev <<- ""
    }
    if (error$message == err_prev) {
      return()
    }

    # Save globally last error
    err_prev <<- error$message

    pgx_name <- NULL
    user_email <- auth$email
    user_tab <- input$nav
    raw_dir <- raw_dir()

    if (!is.null(PGX) && !is.null(PGX$name)) {
      pgx_name <- PGX$name
    } else {
      pgx_name <- "No PGX loaded when error occurred"
    }

    if (is.null(raw_dir)) {
      raw_dir <- "Not a data upload error"
    }

    credential <- file.path(ETC, "hubspot_creds")

    # write dbg statement
    dbg("[SERVER] shiny.error triggered")

    sendErrorLogToCustomerSuport(user_email, pgx_name, raw_dir, error = err_traceback, path_to_creds = credential)
    sever::sever(sever_crash(error), bg_color = "#004c7d")
  })

  # this function sets 'enable_info' based on the user settings
  # and is used by all the bs_alert functions with conditional=T
  observeEvent(env$user_settings$enable_info(), {
    session$sendCustomMessage(
      "enableInfo",
      list(
        id = "enable_info",
        value = env$user_settings$enable_info()
      )
    )
  })

  # this function triggers JS upon 'enable_info' state
  observeEvent(env$user_settings$enable_info(), {
    if (env$user_settings$enable_info() == TRUE) {
      shinyjs::runjs(
        '$(".btn-close-bs-conditional").closest(".bslib-gap-spacing.html-fill-container").css("display", "");$(".btn-close-bs-conditional").parent().css("display", "");'
      )
    } else {
      shinyjs::runjs(
        '$(".btn-close-bs-conditional").closest(".bslib-gap-spacing.html-fill-container").css("display", "none");$(".btn-close-bs-conditional").parent().css("display", "none");'
      )
    }
  })


  ## clean up any remanining UI from previous aborted processx
  shiny::removeUI(selector = "#current_dataset > #spinner-container")

  ## Startup Message
  dbg("[MAIN] showing startup modal")
  observeEvent(auth$logged, {
    if (auth$logged) {
      shinyjs::delay(500, {
        ## read startup messages
        msg <- readLines(file.path(ETC, "MESSAGES"))
        msg <- msg[msg != "" & substr(msg, 1, 1) != "#"]
        msg <- c(msg[[1]], sample(msg, 4))
        STARTUP_MESSAGES <- msg
        shiny::showModal(
          ui.startupModal(
            id = "startup_modal",
            messages = STARTUP_MESSAGES,
            title = "BigOmics Highlights"
          )
        )
      })
    }
  })


  if (isTRUE(opt$ENABLE_INACTIVITY)) {
    # Resest inactivity counter when there is user activity (a click on the UI)
    observeEvent(input$userActivity, {
      inactivityCounter(0) # Reset counter on any user activity
    })

    inactivityControl <- start_inactivityControl(session, timeout = 1800, inactivityCounter)
    observe({
      inactivityControl()
    })
  }

  ## -------------------------------------------------------------
  ## report server times
  ## -------------------------------------------------------------
  server.init_time <- round(Sys.time() - server.start_time, digits = 4)
  info("[SERVER] server.init_time = ", server.init_time, " ", attr(server.init_time, "units"))
  total.lapse_time <- round(Sys.time() - main.start_time, digits = 4)
  info("[SERVER] total lapse time = ", total.lapse_time, " ", attr(total.lapse_time, "units"))
}
