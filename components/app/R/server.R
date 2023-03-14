##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


#' The main application Server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
app_server <- function(input, output, session) {

    message("\n=======================================================================")
    message("================================ SERVER =================================")
    message("=======================================================================\n")

    dbg("[server.R] 0: getwd = ",getwd())
    dbg("[server.R] 0: HONCHO_URL = ",opt$HONCHO_URL)
    dbg("[server.R] 0: SESSION = ",session$token)

    ## Determine is Honcho is alive
    ##honcho.responding <- grepl("Swagger",RCurl::getURL("http://localhost:8000/__docs__/"))
    ##curl.resp <- try(RCurl::getURL("http://localhost:8000/__docs__/"))
    curl.resp <- try(RCurl::getURL(paste0(opt$HONCHO_URL,"/__docs__/")))
    honcho.responding <- grepl("Swagger", curl.resp)
    honcho.responding
    honcho.token <- Sys.getenv("HONCHO_TOKEN", "")
    has.honcho <- (honcho.token!="" && honcho.responding)
    if(1 && has.honcho) {
        info("[server.R] Honcho is alive! ")
        sever::sever(sever_screen2(session$token), bg_color = "#004c7d")
    } else {
        ## No honcho, no email....
        info("[server.R] No Honcho? No party..")
        sever::sever(sever_screen0(), bg_color = "#004c7d") ## lightblue=2780e3
    }

    setwd(WORKDIR)  ## for some reason it can change!! (defined in global.R)
    server.start_time  <- Sys.time()
    session.start_time <- -1
    authentication <- opt$AUTHENTICATION

    limits <- c("samples" = opt$MAX_SAMPLES,
                "comparisons" = opt$MAX_COMPARISONS,
                "genes" = opt$MAX_GENES,
                "genesets" = opt$MAX_GENESETS,
                "datasets" = opt$MAX_DATASETS)
    pgx_dir <- PGX.DIR

    ## Parse and show URL query string
    if(0 && ALLOW_URL_QUERYSTRING) {
        observe({
            query <- parseQueryString(session$clientData$url_search)
            if(length(query)>0) {
                dbg("[server.R:parseQueryString] names.query =",names(query))
                for(i in 1:length(query)) {
                    dbg("[server.R:parseQueryString]",names(query)[i],"=>",query[[i]])
                }
            } else {
                dbg("[server.R:parseQueryString] no queryString!")
            }
            if(!is.null(query[['csv']])) {
                ## focus on this tab
                updateTabsetPanel(session, "load-tabs", selected="Upload data")
                updateTextAreaInput(session, "load-upload_panel-compute-upload_description",
                                    value = "CSV FILE DESCRIPTION")
            }

        })
        dbg("[server.R:parseQueryString] pgx_dir = ",pgx_dir)
    }

    ##-------------------------------------------------------------
    ## Authentication
    ##-------------------------------------------------------------

    auth <- NULL   ## shared in module
    if(authentication == "password") {
        auth <- shiny::callModule(
            PasswordAuthenticationModule, "auth",
            credentials.file = "CREDENTIALS")
    } else if(authentication == "firebase") {
        auth <- shiny::callModule(FirebaseAuthenticationModule, "auth")
    } else if(authentication == "shinyproxy") {
        username <- Sys.getenv("SHINYPROXY_USERNAME")
        ##email <- Sys.getenv("SHINYPROXY_EMAIL")
        auth <- shiny::callModule(NoAuthenticationModule, "auth",
                                  show_modal=TRUE,
                                  username=username, email=username)
    } else if(authentication == "none2") {
        auth <- shiny::callModule(NoAuthenticationModule, "auth",
                                  show_modal=FALSE)
    } else {
        ##} else if(authentication == "none") {
        auth <- shiny::callModule(NoAuthenticationModule, "auth",
                                  show_modal=TRUE)
    }
    dbg("[LoadingBoard] names.auth = ",names(auth))


    ##-------------------------------------------------------------
    ## Call modules
    ##-------------------------------------------------------------

    env <- list()  ## communication "environment"

    ## *** EXPERIMENTAL *** global reactive value replacing env list
    ## above create session global reactiveValue from list
    PGX <- reactiveValues()
    r_global <- reactiveValues(
        load_example_trigger = FALSE,
        reload_pgxdir = 0,
        loadedDataset = 0
    )

    ## Modules needed from the start
    env$load <- LoadingBoard(
        id = "load",
        pgx_dir = pgx_dir,
        pgx = PGX,
        limits = limits,
        auth = auth,
        enable_userdir = opt$ENABLE_USERDIR,
        enable_upload = opt$ENABLE_UPLOAD,
        enable_delete = opt$ENABLE_DELETE,
        enable_save = opt$ENABLE_SAVE,
        r_global = r_global
    )

    ## Modules needed from the start
    env$upload <- UploadBoard(
        id = "upload",
        pgx_dir = pgx_dir,
        pgx = PGX,
        auth = auth,
        limits = limits,
        enable_userdir = opt$ENABLE_USERDIR,
        enable_upload = opt$ENABLE_UPLOAD,
        enable_save = opt$ENABLE_SAVE,
        r_global = r_global
    )

    ## If user is logged off, we clear the data
    observeEvent( auth$logged(), {
        is.logged <- auth$logged()
        length.pgx <- length(names(PGX))
        if(!is.logged && length.pgx>0) {
            for(i in 1:length.pgx) {
                PGX[[names(PGX)[i]]] <<- NULL
            }
        }
    })

    is_data_loaded <- reactive({
        (env$load$loaded() || env$upload$loaded())
    })

    ## Default boards
    WelcomeBoard("welcome", auth=auth, r_global=r_global)
    env$user <- UserBoard("user", user=auth)

    ## Modules needed after dataset is loaded (deferred) --------------
    modules_loaded <- FALSE
    observeEvent( is_data_loaded(), {

        if( is_data_loaded()==0){
            return(NULL)
        }

        if(modules_loaded) {
            Sys.sleep(4)
            shiny::removeModal()  ## remove modal from LoadingBoard
            return(NULL)
        }
        modules_loaded <<- TRUE

        ## TEMPORARY SOLUTION. All modules should use PGX eventually.
        inputData <- reactive({
            if(all(sapply(PGX,is.null))) return(NULL)
            PGX
        })

        shiny::withProgress(message="Preparing your dashboards...", value=0, {

          if(ENABLED['dataview'])  {
            info("[server.R] calling module dataview")
            DataViewBoard("dataview", pgx=PGX)
          }

          if(ENABLED['clustersamples']) {
            info("[server.R] calling module clustersamples")
            ClusteringBoard("clustersamples", pgx=PGX)
          }

          if(ENABLED['wordcloud'])  {
            info("[server.R] calling WordCloudBoard module")
            WordCloudBoard("wordcloud", pgx=PGX)
          }
          shiny::incProgress(0.2)

          if(ENABLED['diffexpr'])   {
            info("[server.R] calling ExpressionBoard module")
            ExpressionBoard("diffexpr", pgx=PGX) -> env$diffexpr
          }

          if(ENABLED['clusterfeatures']) {
            info("[server.R] calling FeatureMapBoard module")
            FeatureMapBoard("clusterfeatures", pgx=PGX)
          }

          if(ENABLED['enrich']) {
            info("[server.R] calling EnrichmentBoard module")
            EnrichmentBoard("enrich", pgx = PGX,
              selected_gxmethods = env$diffexpr$selected_gxmethods ) -> env$enrich
          }
          if(ENABLED['pathway']) {
            info("[server.R] calling FunctionalBoard module")
            FunctionalBoard("pathway", pgx = PGX,
              selected_gsetmethods = env$enrich$selected_gsetmethods)
          }

          shiny::incProgress(0.4)

          if(ENABLED['drug']) {
            info("[server.R] calling DrugConnectivityBoard module")
            DrugConnectivityBoard("drug", pgx = PGX)
          }

          if(ENABLED['isect']) {
            info("[server.R] calling IntersectionBoard module")
            IntersectionBoard("isect", pgx = PGX,
              selected_gxmethods = env$diffexpr$selected_gxmethods,
              selected_gsetmethods = env$enrich$selected_gsetmethods)
          }

          if(ENABLED['sig']) {
            info("[server.R] calling SignatureBoard module")
            SignatureBoard("sig", pgx = PGX,
              selected_gxmethods = env$diffexpr$selected_gxmethods)
          }

          if(ENABLED['corr']) {
            info("[server.R] calling CorrelationBoard module")
            CorrelationBoard("corr", pgx = PGX)
          }
          shiny::incProgress(0.6)

          if(ENABLED['bio']) {
            info("[server.R] calling BiomarkerBoard module")
            BiomarkerBoard("bio", inputData = inputData)
          }

          if(ENABLED['cmap'])  {
            info("[server.R] calling ConnectivityBoard module")
            ConnectivityBoard("cmap", inputData = inputData)
          }

          if(ENABLED['cell']) {
            info("[server.R] calling SingleCellBoard module")
            SingleCellBoard("cell", inputData = inputData)
          }

          shiny::incProgress(0.8)
          if(ENABLED['tcga']) {
            info("[server.R] calling TcgaBoard module")
            TcgaBoard("tcga", inputData = inputData)
          }

          if(ENABLED['wgcna']) {
            info("[server.R] calling WgcnaBoard module")
            WgcnaBoard("wgcna", inputData = inputData)
          }

          if(ENABLED['comp']) {
            info("[server.R] calling CompareBoard module")
            CompareBoard("comp", inputData = inputData)
          }

          info("[server.R] calling modules done!")
        })

        ## remove modal from LoadingBoard
        shiny::removeModal()

        #show hidden tabs
        bigdash.showTabsGoToDataView(session)  # see ui-bigdashplus.R
        ##bigdash.hideTab(session, "upload-tab")

    })


    ##--------------------------------------------------------------------------
    ## Current navigation
    ##--------------------------------------------------------------------------

    output$current_user <- shiny::renderText({
        ## trigger on change of user
        user <- auth$email()
        if(user %in% c("",NA,NULL)) user <- "User"
        user
    })

    output$current_dataset <- shiny::renderText({
        ## trigger on change of dataset
        name <- gsub(".*\\/|[.]pgx$","",PGX$name)
        if(length(name)==0) name = "BigOmics Playground"
        name
    })

    output$current_section <- shiny::renderText({
        cdata <- session$clientData
        section <- sub("section-","",cdata[["url_hash"]])
        section
    })

    ##--------------------------------------------------------------------------
    ## Dynamically hide/show certain sections depending on USERMODE/object
    ##--------------------------------------------------------------------------

    ## toggleTab("load-tabs","Upload data", opt$ENABLE_UPLOAD)
    ##bigdash.toggleTab(session, "upload-tab", opt$ENABLE_UPLOAD)

    shiny::observeEvent({
        auth$logged()
        env$user$enable_beta()
        PGX$name
    }, {

        ## trigger on change dataset
        dbg("[server.R] trigger on change dataset")

        ## show beta feauture
        show.beta <- env$user$enable_beta()
        if(is.null(show.beta) || length(show.beta)==0) show.beta=FALSE
        is.logged <- auth$logged()

        ## hide all main tabs until we have an object
        if(is.null(PGX) || is.null(PGX$name) || !is.logged) {
            warning("[server.R] !!! no data. hiding menu.")
            lapply(MAINTABS, function(m) shiny::hideTab("maintabs",m))
            updateTabsetPanel(session, "maintabs", selected = "Home")
            toggleTab("load-tabs","Upload data",opt$ENABLE_UPLOAD)
            return(NULL)
        }

        ## show all main tabs
        lapply(MAINTABS, function(m) shiny::showTab("maintabs",m))

        ## Beta features
        info("[server.R] disabling beta features")
        toggleTab("drug-tabs","Connectivity map (beta)",show.beta)
        toggleTab("maintabs","TCGA survival (beta)",show.beta,req.file="tcga_matrix.h5")
        ##toggleTab("maintabs","Cluster features",show.beta)
        toggleTab("maintabs","WGCNA (beta)",show.beta)
        toggleTab("maintabs","Compare datasets (beta)",show.beta)

        ## DEVELOPER only tabs (still too alpha)
        info("[server.R] disabling alpha features")
        if(DEV) toggleTab("maintabs","DEV",DEV)
        toggleTab("corr-tabs","Functional",DEV)   ## too slow
        toggleTab("corr-tabs","Differential",DEV)
        toggleTab("dataview-tabs","Resource info",DEV)
        toggleTab("cell-tabs","iTALK",DEV)  ## DEV only
        toggleTab("cell-tabs","CNV",DEV)  ## DEV only
        toggleTab("cell-tabs","Monocle",DEV) ## DEV only
        toggleTab("corr-tabs","Functional",DEV)

        ## Dynamically show upon availability in pgx object
        info("[server.R] disabling extra features")
        tabRequire(PGX, "connectivity", "maintabs", "Similar experiments")
        tabRequire(PGX, "drugs", "maintabs", "Drug connectivity")
        tabRequire(PGX, "wordcloud", "maintabs", "Word cloud")
        tabRequire(PGX, "deconv", "maintabs", "CellProfiling")
        toggleTab("user-tabs","Visitors map",!is.null(ACCESS.LOG))

        info("[server.R] trigger on change dataset done!")
    })

    ##-------------------------------------------------------------
    ## Session TimerModule
    ##-------------------------------------------------------------

    reset_timer <- function() {}
    run_timer <- function(run=TRUE) {}

    if( TIMEOUT > 0 ) {

        rv.timer <- reactiveValues(reset=0, run=FALSE)
        reset_timer <- function() {
            dbg("[server.R] resetting timer")
            rv.timer$reset <- rv.timer$reset + 1
        }
        run_timer <- function(run=TRUE) {
            dbg("[server.R] run timer =",run)
            rv.timer$run <- run
        }
        WARN_BEFORE = round(TIMEOUT/6)

        info("[server.R] Creating TimerModule...")
        info("[server.R] TIMEOUT = ", TIMEOUT)
        info("[server.R] WARN_BEFORE = ", WARN_BEFORE)

        timer <- TimerModule(
            "timer",
            timeout = TIMEOUT,
            warn_before = WARN_BEFORE,
            max_warn = 1,
            poll = Inf,  ## not needed, just for timer output
            reset = reactive(rv.timer$reset),
            run = reactive(rv.timer$run)
        )

        observe({
            info("[server.R] timer = ",timer$timer())
            info("[server.R] lapse_time = ",timer$lapse_time())
        })

        observeEvent( timer$warn(), {
          info("[server.R] timer$warn = ",timer$warn())
          if(timer$warn()==0) return()  ## skip first atInit call
          if(WARN_BEFORE < 60) {
            dt <- paste(round(WARN_BEFORE),"seconds!")
          } else {
            dt <- paste(round(WARN_BEFORE/60),"minutes!")
          }
          showModal(modalDialog(
            HTML("<center><h4>Warning!</h4>Your FREE session is expiring in",dt,".</center>"),
            size = "s",
            easyClose = TRUE
          ))
        })

        r.timeout <- reactive({
          timer$timeout() && auth$logged()
        })

        ## Choose type of referral modal upon timeout:
        mod.timeout <- SocialMediaModule("socialmodal", r.show = r.timeout)
        ##mod.timeout <- SendReferralModule("sendreferral", r.user=auth$name, r.show=r.timeout)

        observeEvent( mod.timeout$success(), {
          success <- mod.timeout$success()
          dbg("[server.R] success = ",success)
          if(success==0) {
            info("[server.R] logout after no referral!!!")
            shinyjs::runjs("logout()")
          }
          if(success > 1) {
            info("[server.R] resetting timer after referral!!!")
            timeout.min <- round(TIMEOUT/60)
            msg = HTML("<center><h4>Thanks!</h4>Your FREE session has been extended.</center>")
            msg = HTML(paste0("<center><h4>Ditch the ",timeout.min,"-minute limit</h4>
Upgrade today and experience advanced analysis features without the time limit.</center>"))


            showModal(modalDialog(
              msg,

              size = "m",
              easyClose = TRUE
            ))
            reset_timer()
          }

        })

      shiny::observeEvent( auth$logged(), {
        ## trigger on change of USER
        logged <- auth$logged()
        info("[server.R & TIMEOUT>0] change in user log status : logged = ",logged)

        ##--------- start timer --------------
        if(TIMEOUT>0 && logged) {
          info("[server.R] starting session timer!!!")
          reset_timer()
          run_timer(TRUE)

        } else {
          info("[server.R] no timer!!!")
          run_timer(FALSE)
        }
      })

    } ## end of if TIMEOUT>0


    ##-------------------------------------------------------------
    ## Session logout functions
    ##-------------------------------------------------------------

    shiny::observe({

        ## trigger on change of USER
        logged <- auth$logged()
        info("[server.R] change in user log status : logged = ",logged)

        ##--------- force logout callback??? --------------
        if(opt$AUTHENTICATION!='firebase' && !logged) {
            ## Forcing logout ensures "clean" sessions. For firebase
            ## we allow sticky sessions.
            message("[server.R] user not logged in? forcing logout() JS callback...")
            shinyjs::runjs("logout()")
        }

    })

    ## logout helper function
    logout.JScallback = "logout()"
    if(opt$AUTHENTICATION=="shinyproxy") {
        logout.JScallback = "function(x){logout();quit();window.location.assign('/logout');}"
    }

    ## This will be called upon user logout *after* the logout() JS call
    observeEvent( input$userLogout, {
      reset_timer()
      run_timer(FALSE)
    })

    ## This code listens to the JS quit signal
    observeEvent( input$quit, {
        dbg("[server.R:quit] closing session... ")
        session$close()
    })

    ## This code will be run after the client has disconnected
    ## Note!!!: Strange behaviour, sudden session ending.
    session$onSessionEnded(function() {
        message("******** doing session cleanup ********")
        ## fill me...
        if(opt$AUTHENTICATION == "shinyproxy") {
            session$sendCustomMessage("shinyproxy-logout", list())
        }

    })

    ##-------------------------------------------------------------
    ## report server times
    ##-------------------------------------------------------------
    server.init_time <- round(Sys.time() - server.start_time, digits=4)
    message("[server.R] server.init_time = ",server.init_time," ",attr(server.init_time,"units"))
    total.lapse_time <- round(Sys.time() - main.start_time,digits=4)
    message("[server.R] total lapse time = ",total.lapse_time," ",attr(total.lapse_time,"units"))

}
