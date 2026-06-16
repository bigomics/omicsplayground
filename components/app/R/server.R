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
  
  ## Initialise the global colour theme (in-session only)
  init_color_theme()

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

  if (!copilot_packages_ok()) {
    shinyalert::shinyalert(
      title = "Missing packages",
      text = "Copilot will be disabled",
      immediate = TRUE,
      showCancelButton = FALSE,
      showConfirmButton = TRUE
    )      
  }

  ## -------------------------------------------------------------
  ## Authentication
  ## -------------------------------------------------------------

  authentication <- opt$AUTHENTICATION
  auth <- NULL ## shared in module
  credentials_file <- file.path(ETC, "CREDENTIALS")
  has.credentials <- file.exists(credentials_file)
  no.credentials <- (!isTRUE(opt$USE_CREDENTIALS) || !has.credentials)
  if (no.credentials && !(authentication %in% c("password", "shinyproxy-sso-admin"))) {
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
  } else if (authentication == "shinyproxy-sso") {
    auth <- AuthenticationModuleHeader(
      id = "auth"
    )
  } else if (authentication == "shinyproxy-sso-admin") {
    auth <- AuthenticationModuleHeader(
      id = "auth",
      credentials_file = credentials_file
    )
  } else {
    ## stop everything
    stop("unsupported authorization method", authentication)
  }

  ## -------------------------------------------------------------
  ## Call modules
  ## -------------------------------------------------------------

  ## Global reactive value for PGX object
  PGX <- reactiveValues()

  ## communication "environment"
  env <- list()
  env$trigger_on_change_dataset <- reactiveVal()

  ## Global reactive values for app-wide triggering
  load_example <- reactiveVal(NULL)
  load_uploaded_data <- reactiveVal(NULL)
  reload_pgxdir <- reactiveVal(0)
  inactivityCounter <- reactiveVal(0)
  new_upload <- reactiveVal(0)

  ## Default boards ------------------------------------------
  ## WelcomeBoard("welcome",
  ##   auth = auth,
  ##   load_example = load_example,
  ##   new_upload = new_upload
  ## )

  env$user_profile <- UserProfileBoard(
    "user_profile",
    auth = auth,
    nav_count = reactive(nav$count)
  )

  ## AppSettingsBoard(
  ##   "app_settings",
  ##   auth = auth,
  ##   pgx = PGX
  ## )

  if (isTRUE(opt$ENABLE_ADMIN)) {
    AdminPanelBoard(
      "admin_panel",
      auth = auth,
      credentials_file = credentials_file
    )
  }

  env$user_settings <- list(
    enable_beta = shiny::reactive(input$enable_beta),
    enable_info = shiny::reactive(input$enable_info)
  )

  ## Do not display "Welcome" tab on the menu
  #bigdash.hideMenuItem(session, "welcome-tab")
  ## Hide admin tab by default (will be shown for admin users after login)
  if (isTRUE(opt$ENABLE_ADMIN)) {
    bigdash.hideMenuItem(session, "admin-tab")
  }
  shinyjs::runjs("sidebarClose()")
  shinyjs::disable(selector = "a[data-value='Dashboard']")
  shinyjs::disable(selector = "a[data-value='Studio']")
  shinyjs::disable(selector = "a[data-value='Copilot']")
  shinyjs::disable(selector = "a[data-value='Copilot2']")  
  
  ## Modules needed from the start
  recompute_pgx <- shiny::reactiveVal(NULL)

  ## --------------------------------------------------------------------------
  ## Navigation
  ## --------------------------------------------------------------------------

  shiny::observeEvent( new_upload(), {
    dbg("[MAIN] reacted! : new_upload = ", new_upload())
    if(new_upload() > 0) {
      bslib::nav_select("app-sidebar", "Upload")
    }
  })

  shinyjs::onclick("logo-bigomics", {
    shinyjs::runjs("console.info('logo-bigomics clicked')")
    #bigdash.selectTab(session, selected = "welcome-tab")
    #shinyjs::runjs("sidebarClose()")
    #shinyjs::runjs("settingsClose()")
  })

  ## observeEvent(input$menu_createreport, {
  ##   shinyjs::click("load-generate_report-show_report_modal")
  ## })

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

  output$current_dataset <- shiny::renderUI({
    has.pgx <- !is.null(PGX$name) && length(PGX$name) > 0
    if (isTRUE(auth$logged) && has.pgx) {
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

  ## Show experiment info if dataset name is clicked.
  observeEvent(input$dataset_click, {
    shiny::req(PGX$name)
    has.infographic <- !is.null(PGX$wgcna$report$infographic)
    pgx.name <- gsub(".*\\/|[.]pgx$", "", PGX$name)    
    if(has.infographic) {
      img <- PGX$wgcna$report$infographic
      footer <- gsub("- |\n"," ",PGX$wgcna$report$bullets)
      footer <- paste("<b>WGCNA graphical abstract</b>. ",footer)
      ui.showImageModal(img, title=NULL, footer, width=1088)         
    } else {
      fields <- c("name", "description", "datatype", "date", "settings",
        "omicsplayground_version")
      body <- playbase::pgx.info(PGX, fields=fields, format="html")
      shiny::showModal(shiny::modalDialog(
        header = pgx.name,
        div(HTML(body), style = "font-size: 1.1em;"),
        footer = NULL,
        size = "l",
        easyClose = TRUE,
        fade = FALSE
      ))
    }
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

  observeEvent(input$navbar_about, ui.showAboutModal())
  
  ## --------------------------------------------------------------------------
  ## upon change of user OR beta toggle OR new pgx
  ## --------------------------------------------------------------------------

  shiny::observeEvent(
    {
      list(
        auth$logged,
        env$user_settings$enable_beta(),
        env$trigger_on_change_dataset()
      )
    },
    {
      ## trigger on change dataset
      shiny::req(PGX$X)
      info("[SERVER] trigger on change dataset")

      ## Set navbar color based on datatype
      ## if (PGX$datatype == "multi-omics") {
      ##   js_code <- sprintf("document.querySelector('.navbar').style.borderBottom = '2px solid #00923b'")
      ##   shinyjs::runjs(js_code)
      ## } else {
      ##   js_code <- sprintf("document.querySelector('.navbar').style.borderBottom = '2px solid #004ca7'")
      ##   shinyjs::runjs(js_code)
      ## }

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
      shiny.i18n::update_lang(lang, session)

      ## tab_control()
      ## ## hide all main tabs until we have an object
      ## if (is.null(PGX) || is.null(PGX$name) || !auth$logged) {
      ##   warning("[SERVER] !!! no data. hiding menu.")
      ##   shinyjs::runjs("sidebarClose()")
      ##   shinyjs::runjs("settingsClose()")
      ##   bigdash.selectTab(session, selected = "dataview-tab")        
      ##   return(NULL)
      ## }
      ## show all main tabs
      #shinyjs::runjs("sidebarOpen()")
      #shinyjs::runjs("settingsOpen()")

      # If is CRO dataset, no watermark\
      cro_emails <- get_cro_emails()
      if (!is.null(PGX$creator) && PGX$creator %in% cro_emails) {
        WATERMARK <<- FALSE
      } else {
        WATERMARK <<- auth$options$WATERMARK
      }

    }
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
  
  ## -------------------------------------------------------------
  ## Session Timers
  ## -------------------------------------------------------------
  ## invite module (from menu)
  invite <- InviteFriendModule(
    id = "invite",
    auth = auth,
    callbackR = inviteCallback
  )
  UpgradeModuleServer(
    id = "upgrade",
    auth = auth
  )
  inviteCallback <- function() {
    ## After succesful invite, we extend the session
    dbg("[MAIN] inviteCallback called!")
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
  #' seconds a small ID file ' in the ONLINE_DIR folder.
  if (isTRUE(opt$ENABLE_HEARTBEAT)) {
    ONLINE_DIR <- file.path(ETC, "online")
    heartbeat <- pgx.start_heartbeat(auth, session, delta = 300, online_dir = ONLINE_DIR)
    observe({
      heartbeat()
    }) ## run indefinitely
  }

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
      sever::sever(sever_max_sessions(opt$HOSTNAME), bg_color = "#004c7d") ## lightblue=2780e3
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

      ## Show/hide upload tab based on user-specific or global ENABLE_UPLOAD option
      ## User option takes precedence over global option if set
      enable_upload <- auth$options$ENABLE_UPLOAD
      if (is.null(enable_upload)) {
        enable_upload <- opt$ENABLE_UPLOAD
      }
      bigdash.toggleMenuItem(session, "upload-tab", isTRUE(enable_upload))
      dbg("[SERVER] ENABLE_UPLOAD for user = ", enable_upload)

      ## Show/hide admin tab based on user's ADMIN status AND global ENABLE_ADMIN option
      if (isTRUE(opt$ENABLE_ADMIN)) {
        is_admin <- isTRUE(auth$ADMIN)
        bigdash.toggleMenuItem(session, "admin-tab", is_admin)
        dbg("[SERVER] ADMIN status for user = ", is_admin)
      }
    } else {
      ## clear PGX data as soon as the user logs out
      clearPGX()
      ## Hide upload tab when logged out (will be re-evaluated on next login)
      bigdash.hideMenuItem(session, "upload-tab")
      ## Hide admin tab when logged out (only if admin is enabled globally)
      if (isTRUE(opt$ENABLE_ADMIN)) {
        bigdash.hideMenuItem(session, "admin-tab")
      }
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
      comment3 = isolate(REPORT_DOWNLOAD_LOGGER$str),
      comment4 = isolate(UPGRADE_LOGGER$str),
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
      if (opt$AUTHENTICATION %in% c("shinyproxy", "shinyproxy-sso", "shinyproxy-sso-admin")) {
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
  shiny::onUnhandledError(function(err) {
    error <- err
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
    board_inputs <- shiny::isolate(names(input)[grep(substr(input$nav, 1, nchar(input$nav) - 4),
      names(input))])

    # Remove pdf + download + card_selector + copy_info + unnecessary table inputs
    board_inputs <- shiny::isolate(board_inputs[-grep("pdf_width|pdf_height|pdf_settings|downloadOption|card_selector|copy_info|_rows_current|_rows_all", board_inputs)])

    input_values <- lapply(board_inputs, function(x) {
      value <- shiny::isolate(input[[x]])
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
    user_email <- shiny::isolate(auth$email)
    user_tab <- shiny::isolate(input$nav)
    raw_dir <- shiny::isolate(raw_dir())

    if (!is.null(PGX) && !is.null(shiny::isolate(PGX$name))) {
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

    if (inherits(error, "shiny.error.fatal")) {
      full_app_crash <- TRUE
    } else {
      full_app_crash <- FALSE
    }

    sendErrorLogToCustomerSuport(user_email, pgx_name, raw_dir, error = err_traceback, path_to_creds = credential, full_app_crash = full_app_crash)
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
  ## dbg("[MAIN] showing startup modal")
  ## observeEvent(auth$logged, {
  ##   if (auth$logged) {
  ##     shinyjs::delay(500, {
  ##       ## skip startup modal if the user has pending shared datasets:
  ##       ## the "New dataset received!" alert takes precedence
  ##       pgx_shared_dir <- stringr::str_replace_all(PGX.DIR, c("data" = "data_shared"))
  ##       has_received <- FALSE
  ##       if (!is.null(auth$email) && nzchar(auth$email) && dir.exists(pgx_shared_dir)) {
  ##         received <- dir(
  ##           path = pgx_shared_dir,
  ##           pattern = paste0("__to__", auth$email, "__from__.*__$"),
  ##           ignore.case = TRUE
  ##         )
  ##         has_received <- length(received) > 0
  ##       }
  ##       ## read startup messages
  ##       ## msg_file <- file.path(ETC, "MESSAGES")
  ##       ## if (!has_received && file.exists(msg_file)) {
  ##       ##   msg <- readLines(msg_file)
  ##       ##   msg <- msg[msg != "" & substr(msg, 1, 1) != "#"]
  ##       ##   if (FALSE  && length(msg) > 0) {
  ##       ##     msg <- c(msg[[1]], sample(msg, min(4, length(msg))))
  ##       ##     STARTUP_MESSAGES <- msg
  ##       ##     shiny::showModal(
  ##       ##       ui.startupModal(
  ##       ##         id = "startup_modal",
  ##       ##         messages = STARTUP_MESSAGES,
  ##       ##         title = "BigOmics Highlights"
  ##       ##       )
  ##       ##     )
  ##       ##   }
  ##       ## }        
  ##     })
  ##   }
  ## })

  if (isTRUE(opt$ENABLE_INACTIVITY)) {
    # Reset inactivity counter when there is user activity (a click on the UI)
    observeEvent(input$userActivity, {
      inactivityCounter(0) # Reset counter on any user activity
    })

    inactivityControl <- start_inactivityControl(session, timeout = opt$INACTIVITY_TIMEOUT, inactivityCounter)
    observe({
      inactivityControl()
    })
  }
  
  observeEvent(input$my_profile, {
    bslib::nav_select("app-sidebar", selected = "UserProfile")
  })
  
  ## -------------------------------------------------------------
  ## Standard modules
  ## -------------------------------------------------------------

  ## Single audited save path for voluntary report regeneration.
  ## Owner check = file must live in auth$user_dir. Public / shared
  ## datasets become a silent no-op.
  save_current_pgx <- function(pgx) {
    if (!isTRUE(auth$logged)) return(invisible(FALSE))
    if (is.null(pgx$name)) return(invisible(FALSE))
    pgxdir <- auth$user_dir
    if (!dir.exists(pgxdir)) return(invisible(FALSE))
    file <- paste0(sub("[.]pgx$", "", pgx$name), ".pgx")
    full <- file.path(pgxdir, file)
    if (!file.exists(full)) return(invisible(FALSE))
    pgx_obj <- if (methods::is(pgx, "reactivevalues")) {
      shiny::reactiveValuesToList(pgx)
    } else {
      pgx
    }
    try(playbase::pgx.save(pgx_obj, file = full), silent = TRUE)
    invisible(TRUE)
  }

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
    new_upload = new_upload,
    save_pgx = save_current_pgx,
    parent = session
  )

  ## Modules needed from the start
  ## NOTE: UploadBoard is always loaded to allow per-user ENABLE_UPLOAD options
  ## The upload tab visibility is controlled after login based on user/global options
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
   
  WelcomeBoard2("welcome2",
    auth = auth,
    load_example = load_example,
    new_upload = new_upload,
    parent = session
  )

  ## -------------------------------------------------------------
  ## Other servers and modules
  ## -------------------------------------------------------------

  opg_server(input, output, session, PGX, env, auth)
  
  CopilotServer("copilot", pgx = PGX, 
    layout = "fixed", maxturns = opt$LLM_MAXTURNS)

  if (copilot_packages_ok()) {
    CopilotBoardServer("copilot2", pgx = PGX, pgx_dir = PGX.DIR,
      auth = auth,
      maxturns = opt$LLM_MAXTURNS,
      tiers = opt$COPILOT_MODEL,
      is_data_loaded = NULL)
  }

  StudioServer("studio", pgx = PGX, save_pgx = save_current_pgx)
  
  AppSettingsBoard("app_settings", auth=auth, pgx=PGX) 
  
  if(opt$DEVMODE) {
    dbg("[SERVER] WARNING: DEVMODE modules enabled!")
    prism_server("prism")
    tools_server("tools", parent = session)
    RunMonitorServer("runmonitor")
  }
  
  ## -------------------------------------------------------------
  ## report server times
  ## -------------------------------------------------------------
  server.init_time <- round(Sys.time() - server.start_time, digits = 4)
  info("[SERVER] server.init_time = ", server.init_time, " ", attr(server.init_time, "units"))
  total.lapse_time <- round(Sys.time() - main.start_time, digits = 4)
  info("[SERVER] total lapse time = ", total.lapse_time, " ", attr(total.lapse_time, "units"))
}
