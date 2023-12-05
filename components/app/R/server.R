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
  authentication <- opt$AUTHENTICATION

  ## -------------------------------------------------------------
  ## Authentication
  ## -------------------------------------------------------------

  auth <- NULL ## shared in module
  credentials_file <- file.path(ETC, "CREDENTIALS")
  has.credentials <- file.exists(credentials_file)
  no.credentials <- (!isTRUE(opt$USE_CREDENTIALS) || !has.credentials)
  if (no.credentials && authentication != "password") {
    credentials_file <- NULL
  }

  if (authentication == "password") {
    auth <- PasswordAuthenticationModule(
      id = "auth",
      credentials_file = credentials_file,
      allow_personal = opt$ALLOW_PERSONAL_EMAIL,
      domain = opt$DOMAIN
    )
  } else if (authentication == "firebase") {
    auth <- FirebaseAuthenticationModule(
      id = "auth",
      domain = opt$DOMAIN,
      firebase.rds = "firebase.rds",
      credentials_file = credentials_file,
      allow_personal = opt$ALLOW_PERSONAL_EMAIL,
      allow_new_users = opt$ALLOW_NEW_USERS
    )
  } else if (authentication == "email-link") {
    auth <- EmailLinkAuthenticationModule(
      id = "auth",
      pgx_dir = PGX.DIR,
      domain = opt$DOMAIN,
      firebase.rds = "firebase.rds",
      credentials_file = credentials_file,
      allow_personal = opt$ALLOW_PERSONAL_EMAIL,
      allow_new_users = opt$ALLOW_NEW_USERS
    )
  } else if (authentication == "login-code") {
    auth <- LoginCodeAuthenticationModule(
      id = "auth",
      mail_creds = file.path(ETC, "gmail_creds"),
      domain = opt$DOMAIN,
      credentials_file = credentials_file,
      allow_personal = opt$ALLOW_PERSONAL_EMAIL,
      allow_new_users = opt$ALLOW_NEW_USERS
    )
  } else if (authentication == "login-code-no-mail") {
    auth <- LoginCodeAuthenticationModule(
      id = "auth",
      mail_creds = file.path(ETC, "gmail_creds"),
      domain = opt$DOMAIN,
      credentials_file = credentials_file,
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
    auth <- AuthenticationModuleApacheCookie(id = "auth", show_modal = FALSE)
  } else {
    auth <- NoAuthenticationModule(id = "auth", show_modal = TRUE)
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
  reload_pgxdir <- reactiveVal(0)
  inactivityCounter <- reactiveVal(0)

  ## Default boards ------------------------------------------
  WelcomeBoard("welcome",
    auth = auth,
    load_example = load_example
  )

  env$user_profile <- UserProfileBoard(
    "user_profile",
    auth = auth,
    nav_count = reactive(nav$count)
  )

  UserSettingsBoard(
    "user_settings",
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
  recompute_info <- shiny::reactiveVal(NULL)
  env$load <- LoadingBoard(
    id = "load",
    pgx = PGX,
    auth = auth,
    pgx_topdir = PGX.DIR,
    load_example = load_example,
    reload_pgxdir = reload_pgxdir,
    current_page = reactive(input$nav),
    load_uploaded_data = load_uploaded_data,
    recompute_pgx = recompute_pgx
  )

  ## Modules needed from the start
  if (opt$ENABLE_UPLOAD) {
    UploadBoard(
      id = "upload",
      pgx_dir = PGX.DIR,
      pgx = PGX,
      auth = auth,
      reload_pgxdir = reload_pgxdir,
      load_uploaded_data = load_uploaded_data,
      recompute_pgx = recompute_pgx,
      recompute_info = recompute_info,
      inactivityCounter = inactivityCounter
    )
  }

  ## Chatbox
  if (opt$ENABLE_CHIRP) {
    shiny::observeEvent(input$chirp_button, {
      shinyjs::click(id = "actual-chirp-button")
    })
    r_chirp_name <- reactive({
      name <- auth$username
      if (is.null(name) || is.na(name) || name == "") name <- auth$email
      if (is.null(name) || is.na(name) || name == "") {
        name <- paste0("user", substring(session$token, 1, 3))
      }
      name <- getFirstName(name) ## in app/R/utils.R
    })
    shinyChatR::chat_server(
      "chatbox",
      csv_path = file.path(SHARE.DIR, "chirp_data.csv"),
      chat_user = r_chirp_name,
      nlast = 100
    )
  }

  ## Modules needed after dataset is loaded (deferred) --------------
  observeEvent(env$load$is_data_loaded(), {
    if (env$load$is_data_loaded() == 1) {
      additional_ui_tabs <- tagList(
        bigdash::bigTabItem(
          "dataview-tab",
          DataViewInputs("dataview"),
          DataViewUI("dataview")
        ),
        bigdash::bigTabItem(
          "clustersamples-tab",
          ClusteringInputs("clustersamples"),
          ClusteringUI("clustersamples")
        ),
        bigdash::bigTabItem(
          "clusterfeatures-tab",
          FeatureMapInputs("clusterfeatures"),
          FeatureMapUI("clusterfeatures")
        ),
        bigdash::bigTabItem(
          "wgcna-tab",
          WgcnaInputs("wgcna"),
          WgcnaUI("wgcna")
        ),
        bigdash::bigTabItem(
          "pcsf-tab",
          PcsfInputs("pcsf"),
          PcsfUI("pcsf")
        ),
        bigdash::bigTabItem(
          "diffexpr-tab",
          ExpressionInputs("diffexpr"),
          ExpressionUI("diffexpr")
        ),
        bigdash::bigTabItem(
          "corr-tab",
          CorrelationInputs("corr"),
          CorrelationUI("corr")
        ),
        bigdash::bigTabItem(
          "enrich-tab",
          EnrichmentInputs("enrich"),
          EnrichmentUI("enrich")
        ),
        bigdash::bigTabItem(
          "pathway-tab",
          PathwayInputs("pathway"),
          PathwayUI("pathway")
        ),
        bigdash::bigTabItem(
          "wordcloud-tab",
          WordCloudInputs("wordcloud"),
          WordCloudUI("wordcloud")
        ),
        bigdash::bigTabItem(
          "drug-tab",
          DrugConnectivityInputs("drug"),
          DrugConnectivityUI("drug")
        ),
        bigdash::bigTabItem(
          "isect-tab",
          IntersectionInputs("isect"),
          IntersectionUI("isect")
        ),
        bigdash::bigTabItem(
          "sig-tab",
          SignatureInputs("sig"),
          SignatureUI("sig")
        ),
        bigdash::bigTabItem(
          "bio-tab",
          BiomarkerInputs("bio"),
          BiomarkerUI("bio")
        ),
        bigdash::bigTabItem(
          "cmap-tab",
          ConnectivityInputs("cmap"),
          ConnectivityUI("cmap")
        ),
        bigdash::bigTabItem(
          "comp-tab",
          CompareInputs("comp"),
          CompareUI("comp")
        ),
        bigdash::bigTabItem(
          "tcga-tab",
          TcgaInputs("tcga"),
          TcgaUI("tcga")
        ),
        bigdash::bigTabItem(
          "cell-tab",
          SingleCellInputs("cell"),
          SingleCellUI("cell")
        ),
      )

      shiny::withProgress(message = "Preparing your dashboard (UI)...", value = 0, {
        shiny::insertUI(
          selector = "#big-tabs",
          where = "beforeEnd",
          ui = additional_ui_tabs,
          immediate = TRUE
        )
      })


      shiny::withProgress(message = "Preparing your dashboard (server)...", value = 0, {
        if (ENABLED["dataview"]) {
          info("[SERVER] calling module dataview")
          DataViewBoard("dataview", pgx = PGX)
        }

        if (ENABLED["clustersamples"]) {
          info("[SERVER] calling ClusteringBoard module")
          ClusteringBoard("clustersamples", pgx = PGX)
        }

        if (ENABLED["wordcloud"]) {
          info("[SERVER] calling WordCloudBoard module")
          WordCloudBoard("wordcloud", pgx = PGX)
        }
        shiny::incProgress(0.2)

        if (ENABLED["diffexpr"]) {
          info("[SERVER] calling ExpressionBoard module")
          ExpressionBoard("diffexpr", pgx = PGX) -> env$diffexpr
        }

        if (ENABLED["clusterfeatures"]) {
          info("[SERVER] calling FeatureMapBoard module")
          FeatureMapBoard("clusterfeatures", pgx = PGX)
        }

        if (ENABLED["enrich"]) {
          info("[SERVER] calling EnrichmentBoard module")
          EnrichmentBoard("enrich",
            pgx = PGX,
            selected_gxmethods = env$diffexpr$selected_gxmethods
          ) -> env$enrich
        }
        if (ENABLED["pathway"]) {
          info("[SERVER] calling PathwayBoard module")
          PathwayBoard("pathway",
            pgx = PGX,
            selected_gsetmethods = env$enrich$selected_gsetmethods
          )
        }

        shiny::incProgress(0.4)

        if (ENABLED["drug"]) {
          info("[SERVER] calling DrugConnectivityBoard module")
          DrugConnectivityBoard("drug", pgx = PGX)
        }

        if (ENABLED["isect"]) {
          info("[SERVER] calling IntersectionBoard module")
          IntersectionBoard("isect",
            pgx = PGX,
            selected_gxmethods = env$diffexpr$selected_gxmethods,
            selected_gsetmethods = env$enrich$selected_gsetmethods
          )
        }

        if (ENABLED["sig"]) {
          info("[SERVER] calling SignatureBoard module")
          SignatureBoard("sig",
            pgx = PGX,
            selected_gxmethods = env$diffexpr$selected_gxmethods
          )
        }

        if (ENABLED["corr"]) {
          info("[SERVER] calling CorrelationBoard module")
          CorrelationBoard("corr", pgx = PGX)
        }
        shiny::incProgress(0.6)

        if (ENABLED["bio"]) {
          info("[SERVER] calling BiomarkerBoard module")
          BiomarkerBoard("bio", pgx = PGX)
        }

        if (ENABLED["cmap"]) {
          info("[SERVER] calling ConnectivityBoard module")
          ConnectivityBoard("cmap",
            pgx = PGX,
            auth = auth,
            reload_pgxdir = reload_pgxdir
          )
        }

        if (ENABLED["cell"]) {
          info("[SERVER] calling SingleCellBoard module")
          SingleCellBoard("cell", pgx = PGX)
        }

        shiny::incProgress(0.8)
        if (ENABLED["tcga"]) {
          info("[SERVER] calling TcgaBoard module")
          TcgaBoard("tcga", pgx = PGX)
        }

        if (ENABLED["wgcna"]) {
          info("[SERVER] calling WgcnaBoard module")
          WgcnaBoard("wgcna", pgx = PGX)
        }

        if (ENABLED["pcsf"]) {
          info("[SERVER] calling PcsfBoard module")
          PcsfBoard("pcsf", pgx = PGX)
        }

        if (ENABLED["comp"]) {
          info("[SERVER] calling CompareBoard module")
          CompareBoard("comp", pgx = PGX, pgx_dir = reactive(auth$user_dir))
        }

        info("[SERVER] calling modules done!")
      })
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

  output$current_dataset <- shiny::renderText({
    shiny::req(auth$logged)
    has.pgx <- !is.null(PGX$name) && length(PGX$name) > 0
    nav.welcome <- welcome_detector()

    if (isTRUE(auth$logged) && has.pgx && !nav.welcome) {
      ## trigger on change of dataset
      name <- gsub(".*\\/|[.]pgx$", "", PGX$name)
    } else {
      name <- paste("Omics Playground", VERSION)
    }
    name
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

  ## --------------------------------------------------------------------------

  ## upon change of user OR beta toggle OR new pgx
  shiny::observeEvent(
    {
      auth$logged
      env$user_settings$enable_beta()
      PGX$name
    },
    {
      ## trigger on change dataset
      dbg("[SERVER] trigger on change dataset")

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
      bigdash.toggleTab(session, "pcsf-tab", show.beta) ## wgcna
      bigdash.toggleTab(session, "tcga-tab", show.beta && has.libx)
      toggleTab("drug-tabs", "Connectivity map (beta)", show.beta) ## too slow
      toggleTab("pathway-tabs", "Enrichment Map (beta)", show.beta) ## too slow

      ## Dynamically show upon availability in pgx object
      info("[SERVER] disabling extra features")
      tabRequire(PGX, session, "wgcna-tab", "wgcna", TRUE)
      ##      tabRequire(PGX, session, "cmap-tab", "connectivity", has.libx)
      tabRequire(PGX, session, "drug-tab", "drugs", TRUE)
      tabRequire(PGX, session, "wordcloud-tab", "wordcloud", TRUE)
      tabRequire(PGX, session, "cell-tab", "deconv", TRUE)

      ## DEVELOPER only tabs (still too alpha)
      info("[SERVER] disabling alpha features")
      #      toggleTab("cell-tabs", "iTALK", DEV) ## DEV only
      #      toggleTab("cell-tabs", "CNV", DEV) ## DEV only
      #      toggleTab("cell-tabs", "Monocle", DEV) ## DEV only

      info("[SERVER] trigger on change dataset done!")
    }
  )

  ## -------------------------------------------------------------
  ## Session Timers
  ## -------------------------------------------------------------

  session_timer <- NULL
  if (isTRUE(TIMEOUT > 0)) {
    #' Session timer. Closes session after TIMEOUT (seconds) This
    #' is acitve for free users. Set TIMEOUT=0 to disable session
    #' timer.
    session_timer <- TimerModule(
      "session_timer",
      condition = reactive(auth$logged),
      timeout = TIMEOUT,
      warn_before = round(0.15 * TIMEOUT),
      max_warn = 1
    )

    observeEvent(session_timer$warn_event(), {
      if (session_timer$warn_event() == 0) {
        return()
      } ## skip first atInit call
      shinyalert::shinyalert(
        title = "Warning!",
        text = "Your FREE session is expiring soon",
        html = TRUE,
        immediate = TRUE,
        timer = 3000
      )
    })

    ## At the end of the timeout the user can choose type of referral
    ## modal and gain additional analysis time. We reset the timer.
    r.timeout <- reactive({
      session_timer$timeout_event() && auth$logged
    })
    social <- SocialMediaModule("socialmodal", r.show = r.timeout)
    social$start_shiny_observer(session_timer$reset)
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
      poll_secs = 15,
      max_idle = 60,
      show_success = FALSE,
      show_details = FALSE
    )
    lock$start_shiny_observer(auth, session = session)
  }

  #' Track which users are online by repeatedly writing small ID file
  #' in the ONLINE_DIR folder.
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
      "Ana Nufer, Axel Martinelli, Carson Sievert, Cédric Scherer, Gabriela Scorici, Ivo Kwee, John Coene, Layal Abo Khayal, Marco Sciaini, Matt Leech, Mauro Miguel Masiero, Murat Akhmedov, Nick Cullen, Santiago Caño Muñiz, Shalini Pandurangan, Stefan Reifenberg, Xavier Escribà Montagut"
    )
    authors <- paste(sort(authors), collapse = ", ")

    shiny::showModal(
      shiny::modalDialog(
        div(
          h2("BigOmics Playground"),
          h5(VERSION),
          h5("Advanced omics analysis for everyone"), br(), br(),
          p("Created with love and proudly presented to you by BigOmics Analytics from Ticino,
               the sunny side of Switzerland."),
          p(tags$a(href = "https://www.bigomics.ch", "www.bigomics.ch")),
          style = "text-align:center; line-height: 1em;"
        ),
        footer = div(
          "© 2000-2023 BigOmics Analytics, Inc.",
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
      enable_upload <- auth$options$ENABLE_UPLOAD
      bigdash.toggleTab(session, "upload-tab", enable_upload)
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

    pgx.record_access(
      user = isolate(auth$email),
      action = action,
      session = session,
      comment = nav_count.str
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
  if (!is.null(opt$STARTUP_MESSAGE) && opt$STARTUP_MESSAGE[1] != "") {
    shinyalert::shinyalert(
      title = HTML(opt$STARTUP_TITLE),
      text = HTML(opt$STARTUP_MESSAGE),
      html = TRUE
    )
  }

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
