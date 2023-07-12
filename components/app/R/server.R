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
  message("\n=======================================================================")
  message("================================ SERVER ===============================")
  message("=======================================================================\n")

  VERSION <- scan(file.path(OPG, "VERSION"), character())[1]

  info("[server.R] getwd = ", normalizePath(getwd()))
  info("[server.R] HONCHO_URL = ", opt$HONCHO_URL)
  info("[server.R] SESSION = ", session$token)

  has.honcho <- FALSE
  if (has.honcho && !is.null(opt$HONCHO_URL) && opt$HONCHO_URL != "") {
    ## Determine is Honcho is alive
    curl.resp <- try(RCurl::getURL(paste0(opt$HONCHO_URL, "/__docs__/")))
    honcho.responding <- grepl("Swagger", curl.resp)
    honcho.token <- Sys.getenv("HONCHO_TOKEN", "")
    has.honcho <- (honcho.token != "" && honcho.responding)
  }
  if (1 && has.honcho) {
    info("[server.R] Honcho is alive! ")
    sever::sever(sever_screen2(session$token), bg_color = "#004c7d")
  } else {
    ## No honcho, no email....
    info("[server.R] No Honcho? No party..")
    sever::sever(sever_screen0(), bg_color = "#004c7d")
  }

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
  if ((is.null(opt$USE_CREDENTIALS) || !opt$USE_CREDENTIALS ||
    !has.credentials) && authentication != "password") {
    credentials_file <- NULL
  }

  if (authentication == "password") {
    auth <- PasswordAuthenticationModule(
      id = "auth",
      credentials_file = credentials_file
    )
  } else if (authentication == "firebase") {
    auth <- FirebaseAuthenticationModule(
      id = "auth",
      domain = opt$DOMAIN,
      firebase.rds = "firebase.rds"
    )
  } else if (authentication == "email-link") {
    auth <- EmailLinkAuthenticationModule(
      id = "auth",
      pgx_dir = PGX.DIR,
      domain = opt$DOMAIN,
      credentials_file = credentials_file,
      firebase.rds = "firebase.rds"
    )
  } else if (authentication == "login-code") {
    auth <- LoginCodeAuthenticationModule(
      id = "auth",
      mail_creds = file.path(ETC, "gmail_creds"),
      domain = opt$DOMAIN,
      credentials_file = credentials_file
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
    auth <- NoAuthenticationModule(id = "auth", show_modal = FALSE)
  } else {
    ## } else if(authentication == "none") {
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
  reload_pgxdir <- reactiveVal(NULL)

  ## Default boards ------------------------------------------
  WelcomeBoard("welcome",
    auth = auth,
    load_example = load_example
  )
  env$user_profile <- UserProfileBoard("user_profile", auth = auth)
  env$user_settings <- UserSettingsBoard("user_settings", auth = auth)

  ## Do not display "Welcome" tab on the menu
  bigdash.hideMenuItem(session, "welcome-tab")
  shinyjs::runjs("sidebarClose()")

  ## Modules needed from the start
  env$load <- LoadingBoard(
    id = "load",
    pgx = PGX,
    auth = auth,
    pgx_topdir = PGX.DIR,
    load_example = load_example,
    reload_pgxdir = reload_pgxdir,
    current_page = reactive(input$nav),
    load_uploaded_data = load_uploaded_data
  )

  ## Modules needed from the start
  if (opt$ENABLE_UPLOAD) {
    UploadBoard(
      id = "upload",
      pgx_dir = PGX.DIR,
      pgx = PGX,
      auth = auth,
      reload_pgxdir = reload_pgxdir,
      load_uploaded_data = load_uploaded_data
    )
  }

  ## Chatbox
  if (opt$ENABLE_CHIRP) {
    observeEvent(input$chirp_button, {
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
      ## db_file = file.path(ETC, "chirp_data.db"),
      csv_path = file.path(ETC, "chirp_data.csv"),
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
          FunctionalInputs("pathway"),
          FunctionalUI("pathway")
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
        )
      )

      shiny::withProgress(message = "Preparing your dashboard (UI)...", value = 0, {
        shiny::insertUI(
          selector = "#big-tabs",
          where = "beforeEnd",
          ui = additional_ui_tabs,
          immediate = TRUE
        )
      })
    }

    shiny::withProgress(message = "Preparing your dashboard (server)...", value = 0, {
      if (ENABLED["dataview"]) {
        info("[server.R] calling module dataview")
        DataViewBoard("dataview", pgx = PGX)
      }

      if (ENABLED["clustersamples"]) {
        info("[server.R] calling module clustersamples")
        ClusteringBoard("clustersamples", pgx = PGX)
      }

      if (ENABLED["wordcloud"]) {
        info("[server.R] calling WordCloudBoard module")
        WordCloudBoard("wordcloud", pgx = PGX)
      }
      shiny::incProgress(0.2)

      if (ENABLED["diffexpr"]) {
        info("[server.R] calling ExpressionBoard module")
        ExpressionBoard("diffexpr", pgx = PGX) -> env$diffexpr
      }

      if (ENABLED["clusterfeatures"]) {
        info("[server.R] calling FeatureMapBoard module")
        FeatureMapBoard("clusterfeatures", pgx = PGX)
      }

      if (ENABLED["enrich"]) {
        info("[server.R] calling EnrichmentBoard module")
        EnrichmentBoard("enrich",
          pgx = PGX,
          selected_gxmethods = env$diffexpr$selected_gxmethods
        ) -> env$enrich
      }
      if (ENABLED["pathway"]) {
        info("[server.R] calling FunctionalBoard module")
        FunctionalBoard("pathway",
          pgx = PGX,
          selected_gsetmethods = env$enrich$selected_gsetmethods
        )
      }

      shiny::incProgress(0.4)

      if (ENABLED["drug"]) {
        info("[server.R] calling DrugConnectivityBoard module")
        DrugConnectivityBoard("drug", pgx = PGX)
      }

      if (ENABLED["isect"]) {
        info("[server.R] calling IntersectionBoard module")
        IntersectionBoard("isect",
          pgx = PGX,
          selected_gxmethods = env$diffexpr$selected_gxmethods,
          selected_gsetmethods = env$enrich$selected_gsetmethods
        )
      }

      if (ENABLED["sig"]) {
        info("[server.R] calling SignatureBoard module")
        SignatureBoard("sig",
          pgx = PGX,
          selected_gxmethods = env$diffexpr$selected_gxmethods
        )
      }

      if (ENABLED["corr"]) {
        info("[server.R] calling CorrelationBoard module")
        CorrelationBoard("corr", pgx = PGX)
      }
      shiny::incProgress(0.6)

      if (ENABLED["bio"]) {
        info("[server.R] calling BiomarkerBoard module")
        BiomarkerBoard("bio", pgx = PGX)
      }

      if (ENABLED["cmap"]) {
        info("[server.R] calling ConnectivityBoard module")
        ConnectivityBoard("cmap", pgx = PGX, getPgxDir = auth$options$user_dir)
      }

      if (ENABLED["cell"]) {
        info("[server.R] calling SingleCellBoard module")
        SingleCellBoard("cell", pgx = PGX)
      }

      shiny::incProgress(0.8)
      if (ENABLED["tcga"]) {
        info("[server.R] calling TcgaBoard module")
        TcgaBoard("tcga", pgx = PGX)
      }

      if (ENABLED["wgcna"]) {
        info("[server.R] calling WgcnaBoard module")
        WgcnaBoard("wgcna", pgx = PGX)
      }

      if (ENABLED["pcsf"]) {
        info("[server.R] calling PcsfBoard module")
        PcsfBoard("pcsf", pgx = PGX)
      }

      if (ENABLED["comp"]) {
        info("[server.R] calling CompareBoard module")
        CompareBoard("comp", pgx = PGX)
      }

      info("[server.R] calling modules done!")
    })

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
    if (!auth$logged) {
      return("(not logged in)")
    }
    user <- auth$email
    if (is.null(user) || user %in% c("", NA)) user <- auth$username
    if (is.null(user) || user %in% c("", NA)) user <- "User"
    dbg("[server:output$current_user] user = ", user)
    user
  })

  output$current_dataset <- shiny::renderText({
    ## trigger on change of dataset
    name <- gsub(".*\\/|[.]pgx$", "", PGX$name)
    if (length(name) == 0) name <- paste("Omics Playground", VERSION)
    name
  })

  ## --------------------------------------------------------------------------
  ## Dynamically hide/show certain sections depending on USERMODE/object
  ## --------------------------------------------------------------------------

  ## upon change of user
  observeEvent(auth$logged, {
    if (auth$logged) {
      enable_upload <- auth$options$ENABLE_UPLOAD
      bigdash.toggleTab(session, "upload-tab", enable_upload)

      # check personal email for old users and ask them to change
      # their email
      if (auth$method %in% c("email-link", "firebase", "login-code")) {
        check_personal_email(auth, PGX.DIR)
      }
    } else {
      # clear PGX data as soon as the user logs out
      length.pgx <- length(names(PGX))
      if (length.pgx > 0) {
        for (i in 1:length.pgx) {
          PGX[[names(PGX)[i]]] <<- NULL
        }
      }
    }
  })

  ## upon change of user OR beta toggle OR new pgx
  shiny::observeEvent(
    {
      auth$logged
      env$user_settings$enable_beta()
      PGX$name
    },
    {
      ## trigger on change dataset
      dbg("[server.R] trigger on change dataset")

      ## show beta feauture
      show.beta <- env$user_settings$enable_beta()
      if (is.null(show.beta) || length(show.beta) == 0) show.beta <- FALSE
      is.logged <- auth$logged

      ## hide all main tabs until we have an object
      if (is.null(PGX) || is.null(PGX$name) || !is.logged) {
        warning("[server.R] !!! no data. hiding menu.")
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
      info("[server.R] disabling beta features")
      bigdash.toggleTab(session, "pcsf-tab", show.beta) ## wgcna
      bigdash.toggleTab(session, "tcga-tab", show.beta && has.libx)
      toggleTab("drug-tabs", "Connectivity map (beta)", show.beta) ## too slow
      toggleTab("pathway-tabs", "Enrichment Map (beta)", show.beta) ## too slow
      toggleTab("dataview-tabs", "Resource info", show.beta)

      ## Dynamically show upon availability in pgx object
      info("[server.R] disabling extra features")
      tabRequire(PGX, session, "wgcna-tab", "wgcna", TRUE)
      ##      tabRequire(PGX, session, "cmap-tab", "connectivity", has.libx)
      tabRequire(PGX, session, "drug-tab", "drugs", TRUE)
      tabRequire(PGX, session, "wordcloud-tab", "wordcloud", TRUE)
      tabRequire(PGX, session, "cell-tab", "deconv", TRUE)

      ## DEVELOPER only tabs (still too alpha)
      info("[server.R] disabling alpha features")
      #      toggleTab("cell-tabs", "iTALK", DEV) ## DEV only
      #      toggleTab("cell-tabs", "CNV", DEV) ## DEV only
      #      toggleTab("cell-tabs", "Monocle", DEV) ## DEV only

      info("[server.R] trigger on change dataset done!")
    }
  )

  ## -------------------------------------------------------------
  ## Session TimerModule
  ## -------------------------------------------------------------

  reset_timer <- function() {}
  run_timer <- function(run = TRUE) {}

  if (TIMEOUT > 0) {
    rv.timer <- reactiveValues(reset = 0, run = FALSE)
    reset_timer <- function() {
      dbg("[server.R] resetting timer")
      rv.timer$reset <- rv.timer$reset + 1
    }
    run_timer <- function(run = TRUE) {
      dbg("[server.R] run timer =", run)
      rv.timer$run <- run
    }
    WARN_BEFORE <- round(TIMEOUT / 6)

    info("[server.R] Creating TimerModule: TIMEOUT = ", TIMEOUT, "(s)")
    info("[server.R] WARN_BEFORE = ", WARN_BEFORE)

    timer <- TimerModule(
      "timer",
      timeout = TIMEOUT,
      warn_before = WARN_BEFORE,
      max_warn = 1,
      poll = Inf, ## not needed, just for timer output
      reset = reactive(rv.timer$reset),
      run = reactive(rv.timer$run)
    )

    observeEvent(timer$warn(), {
      if (timer$warn() == 0) {
        return()
      } ## skip first atInit call
      if (WARN_BEFORE < 60) {
        dt <- paste(round(WARN_BEFORE), "seconds!")
      } else {
        dt <- paste(round(WARN_BEFORE / 60), "minutes!")
      }
      showModal(modalDialog(
        HTML("<center><h4>Warning!</h4>Your FREE session is expiring in", dt, ".</center>"),
        size = "s",
        easyClose = TRUE
      ))
    })

    r.timeout <- reactive({
      timer$timeout() && auth$logged
    })

    ## Choose type of referral modal upon timeout:
    mod.timeout <- SocialMediaModule("socialmodal", r.show = r.timeout)

    observeEvent(mod.timeout$success(), {
      success <- mod.timeout$success()
      dbg("[server.R] success = ", success)
      if (success == 0) {
        info("[server.R] logout after no referral!!!")
        shinyjs::runjs("logoutInApp()")
      }
      if (success > 1) {
        info("[server.R] resetting timer after referral!!!")
        timeout.min <- round(TIMEOUT / 60)
        msg <- HTML("<center><h4>Thanks!</h4>Your FREE session has been extended.</center>")
        msg <- HTML(paste0("<center><h4>Ditch the ", timeout.min, "-minute limit</h4>
Upgrade today and experience advanced analysis features without the time limit.</center>"))

        showModal(modalDialog(
          msg,
          size = "m",
          easyClose = TRUE
        ))
        reset_timer()
      }
    })

    shiny::observeEvent(auth$logged, {
      ## trigger on change of USER
      logged <- auth$logged
      info("[server.R & TIMEOUT>0] change in user log status : logged = ", logged)

      ## --------- start timer --------------
      if (TIMEOUT > 0 && logged) {
        info("[server.R] starting session timer!!!")
        reset_timer()
        run_timer(TRUE)
      } else {
        info("[server.R] no timer!!!")
        run_timer(FALSE)
      }
    })
  } ## end of if TIMEOUT>0


  ## -------------------------------------------------------------
  ## About
  ## -------------------------------------------------------------

  observeEvent(input$navbar_about, {
    authors <- c(
      "Ivo Kwee", "Murat Akhmedov", "John Coene",
      "Stefan Reifenberg", "Marco Sciaini", "Cédric Scherer",
      "Mauro Miguel Masiero", "Nick Cullen", "Layal Abo Khayal",
      "Xavier Escribà Montagut", "Carson Sievert", "Matt Leech"
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
  ## Session logout functions
  ## -------------------------------------------------------------

  ## logout helper function
  logout.JScallback <- "logout()"
  if (opt$AUTHENTICATION == "shinyproxy") {
    logout.JScallback <- "function(x){logout();quit();window.location.assign('/logout');}"
  }

  ## This will be called upon user logout *after* the logout() JS call
  observeEvent(input$userLogout, {
    reset_timer()
    run_timer(FALSE)
  })

  ## This code listens to the JS quit signal
  observeEvent(input$quit, {
    dbg("[server.R:quit] closing session... ")
    session$close()
  })

  # This code will run when there is a shiny error. Then this
  # error will be shown on the app. Note that errors that are
  # not related to Shiny are not caught (e.g. an error on the
  # global.R file is not caught by this)
  options(shiny.error = function() {
    # The error message is on the parent environment, it is
    # not passed to the function called on error
    parent_env <- parent.frame()
    error <- parent_env$e
    sever::sever(sever_screen0(
      error
    ), bg_color = "#004c7d")
  })

  ## This code will be run after the client has disconnected
  ## Note!!!: Strange behaviour, sudden session ending.
  session$onSessionEnded(function() {
    message("******** doing session cleanup ********")
    ## fill me...
    if (opt$AUTHENTICATION == "shinyproxy") {
      session$sendCustomMessage("shinyproxy-logout", list())
    }
  })

  # this function sets 'input.enable_info' based on the user settings
  # and is used by all the bs_alert functions in a conditionalPanel
  observeEvent(env$user_settings$enable_info(), {
    session$sendCustomMessage(
      "enableInfo",
      list(
        id = "enable_info",
        value = env$user_settings$enable_info()
      )
    )
  })

  ## -------------------------------------------------------------
  ## report server times
  ## -------------------------------------------------------------
  server.init_time <- round(Sys.time() - server.start_time, digits = 4)
  message("[server.R] server.init_time = ", server.init_time, " ", attr(server.init_time, "units"))
  total.lapse_time <- round(Sys.time() - main.start_time, digits = 4)
  message("[server.R] total lapse time = ", total.lapse_time, " ", attr(total.lapse_time, "units"))


  ## clean up any remanining UI from previous aborted processx
  shiny::removeUI(selector = ".current-dataset > #spinner-container")

  ## Startup Message
  if (!is.null(opt$STARTUP_MESSAGE) && opt$STARTUP_MESSAGE != "") {
    shinyalert::shinyalert(
      title = HTML(opt$STARTUP_TITLE),
      text = HTML(opt$STARTUP_MESSAGE),
      html = TRUE
    )
  }
}
