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

    dbg("[SERVER] 0: getwd = ",getwd())
    dbg("[SERVER] 0: HONCHO_URL = ",opt$HONCHO_URL)
    dbg("[SERVER] 0: SESSION = ",session$token)
    
    ## Logging of input/output events -------------------------------------
    log.path <- "../logs/"
    log.path <- file.path(OPG,"logs")
    dbg("[SERVER] shinylog log path = ",log.path)
    ## shinylogs::track_usage(storage_mode = shinylogs::store_rds(path = log.path))

    ## Determine is Honcho is alive
    ##honcho.responding <- grepl("Swagger",RCurl::getURL("http://localhost:8000/__docs__/"))
    ##curl.resp <- try(RCurl::getURL("http://localhost:8000/__docs__/"))
    curl.resp <- try(RCurl::getURL(paste0(opt$HONCHO_URL,"/__docs__/")))
    honcho.responding <- grepl("Swagger", curl.resp)
    honcho.responding      
    honcho.token <- Sys.getenv("HONCHO_TOKEN", "")
    has.honcho <- (honcho.token!="" && honcho.responding)
    if(1 && has.honcho) {
        dbg("[SERVER] Honcho is alive! ")    
        sever::sever(sever_screen2(session$token), bg_color = "#000000") 
    } else {
        ## No honcho, no email....
        dbg("[SERVER] No Honcho? No party..")          
        sever::sever(sever_screen0(), bg_color = "#000000") ## lightblue=2780e3
    }

    setwd(WORKDIR)  ## for some reason it can change!! (defined in global.R)
    dbg("[SERVER] 1: getwd = ",getwd())
    
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
                dbg("[SERVER:parseQueryString] names.query =",names(query))
                for(i in 1:length(query)) {
                    dbg("[SERVER:parseQueryString]",names(query)[i],"=>",query[[i]])
                }
            } else {
                dbg("[SERVER:parseQueryString] no queryString!")
            }            
            if(!is.null(query[['csv']])) {
                ## focus on this tab
                updateTabsetPanel(session, "load-tabs", selected="Upload data")
                updateTextAreaInput(session, "load-upload_panel-compute-upload_description",
                                    value = "CSV FILE DESCRIPTION")                
            }
            
        })
        dbg("[SERVER:parseQueryString] pgx_dir = ",pgx_dir)
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
        enable_save = opt$ENABLE_SAVE
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
        enable_save = opt$ENABLE_SAVE
    )   
  
    ## If user is logged off, we clear the data
    observeEvent( auth$logged(), {
        is.logged <- auth$logged()
        length.pgx <- length(names(PGX))
        dbg("[SERVER] *** clearing PGX ***")
        if(!is.logged && length.pgx>0) {
            for(i in 1:length.pgx) {
                PGX[[names(PGX)[i]]] <<- NULL
            }
        }
    })

    data_loaded <- reactive({
        (env$load$loaded() || env$upload$loaded())
    })
  
    ## Default boards
    WelcomeBoard("welcome", auth=auth) 
    UserBoard("user", user=auth) -> env$user  
    
    ## Modules needed after dataset is loaded (deferred) --------------
    modules_loaded <- FALSE
    observeEvent( data_loaded(), {        

        message("[SERVER:data.loaded] data_loaded = ",data_loaded())    
        if(data_loaded()==0){
            return(NULL)
        }
        
        if(modules_loaded) {
            Sys.sleep(4)
            shiny::removeModal()  ## remove modal from LoadingBoard
            return(NULL)
        }
        modules_loaded <<- TRUE

        ## load other modules if not yet loaded
        message("[SERVER] --------- calling shiny modules ----------")
        dbg("[SERVER] names(pgx) = ",names(PGX))        

        loadModule <- function(...) {
            id <- list(...)[[2]]
            if(ENABLED[id])  env[[id]] <<- shiny::callModule(...)
        }
        
        ## TEMPORARY SOLUTION. All modules should use PGX eventually.
        inputData <- reactive({
            if(all(sapply(PGX,is.null))) return(NULL)
            PGX
        })
        
        shiny::withProgress(message="initializing modules ...", value=0, {

            DataViewBoard("view", pgx=PGX)            
            ClusteringBoard("clust", pgx=PGX)
            WordCloudBoard("word", pgx=PGX)
            shiny::incProgress(0.2)

            ## *** DEVNOTE *** board below still need refactoring
            ExpressionBoard("expr", inputData=inputData) -> env$expr
            FeatureMapBoard("ftmap", inputData=inputData)
            EnrichmentBoard("enrich", inputData = inputData,
                            selected_gxmethods = env$expr$selected_gxmethods
                            ) -> env$enrich
            FunctionalBoard("func", inputData = inputData,
                            selected_gsetmethods = env$enrich$selected_gsetmethods)
            shiny::incProgress(0.4)
            DrugConnectivityBoard("drug", inputData = inputData)
            IntersectionBoard("isect", inputData = inputData,
                              selected_gxmethods = env$expr$selected_gxmethods,
                              selected_gsetmethods = env$enrich$selected_gsetmethods)
            SignatureBoard("sig", inputData = inputData,
                           selected_gxmethods = env$expr$selected_gxmethods)
            CorrelationBoard("cor", inputData = inputData)
            shiny::incProgress(0.6)            
            BiomarkerBoard("bio", inputData = inputData)
            ConnectivityBoard("cmap", inputData = inputData)
            SingleCellBoard("scell", inputData = inputData)
            shiny::incProgress(0.8)
            TcgaBoard("tcga", inputData = inputData)
            WgcnaBoard("wgcna", inputData = inputData)
            CompareBoard("comp", inputData = inputData)
            
        })

        message("[SERVER:data_loaded] --------- done! ----------")
        ## remove modal from LoadingBoard
        shiny::removeModal()

        session$sendCustomMessage(
            "show-tabs",
            list()
        )
    })
    

    ##--------------------------------------------------------------------------
    ## Current navigation
    ##--------------------------------------------------------------------------
    
    output$current_user <- shiny::renderText({
        ## trigger on change of user
        user <- auth$email()
        dbg("[SERVER:output$current_user] user = ",user)
        if(user %in% c("",NA,NULL)) user <- "User"
        user
    })
    
    output$current_dataset <- shiny::renderText({
        ## trigger on change of dataset
        name <- gsub(".*\\/|[.]pgx$","",PGX$name)
        dbg("[SERVER:output$current_dataset] dataset = ",name)
        if(length(name)==0) name = "BigOmics Playground"
        name
    })

    output$current_section <- shiny::renderText({
        cdata <- session$clientData
        section <- sub("section-","",cdata[["url_hash"]])
        dbg("[SERVER:output$current_section] section = ",section)
        section
    })
    
    ##--------------------------------------------------------------------------
    ## Dynamically hide/show certain sections depending on USERMODE/object
    ##--------------------------------------------------------------------------

    shiny::observeEvent({
        auth$logged()        
        env$user$enable_beta()
        PGX$name
    }, {

        message("[SERVER] !!! dataset changed. reconfiguring triggered!")
        ## trigger on change dataset

        ## show beta feauture
        show.beta <- env$user$enable_beta()
        dbg("[SERVER] show.beta = ",show.beta)
        if(is.null(show.beta) || length(show.beta)==0) show.beta=FALSE
        is.logged <- auth$logged()
        
        ## hide all main tabs until we have an object
        if(is.null(PGX) || is.null(PGX$name) || !is.logged) {
            message("[SERVER] !!! no data. hiding menu.")          
            lapply(MAINTABS, function(m) shiny::hideTab("maintabs",m))
            updateTabsetPanel(session, "maintabs", selected = "Home")                        
            toggleTab("load-tabs","Upload data",opt$ENABLE_UPLOAD)
            return(NULL)
        }
        
        message("[SERVER] dataset changed. reconfiguring menu...")
        ## show all main tabs
        lapply(MAINTABS, function(m) shiny::showTab("maintabs",m))
        
        ## Beta features
        toggleTab("drug-tabs","Connectivity map (beta)",show.beta)                
        toggleTab("maintabs","TCGA survival (beta)",show.beta,req.file="tcga_matrix.h5")
        ##toggleTab("maintabs","Cluster features",show.beta)
        toggleTab("maintabs","WGCNA (beta)",show.beta)
        toggleTab("maintabs","Compare datasets (beta)",show.beta)        
        
        ## DEVELOPER only tabs (still too alpha)
        if(DEV) toggleTab("maintabs","DEV",DEV)
        toggleTab("cor-tabs","Functional",DEV)   ## too slow
        toggleTab("cor-tabs","Differential",DEV)  
        toggleTab("view-tabs","Resource info",DEV)
        toggleTab("scell-tabs","iTALK",DEV)  ## DEV only
        toggleTab("scell-tabs","CNV",DEV)  ## DEV only        
        toggleTab("scell-tabs","Monocle",DEV) ## DEV only
        toggleTab("cor-tabs","Functional",DEV)
        
        ## Dynamically show upon availability in pgx object
        toggleTab("load-tabs","Upload data", opt$ENABLE_UPLOAD)            
        tabRequire(PGX, "connectivity", "maintabs", "Find similar experiments")
        tabRequire(PGX, "drugs", "maintabs", "Drug connectivity")
        tabRequire(PGX, "wordcloud", "maintabs", "Word cloud")
        tabRequire(PGX, "deconv", "maintabs", "CellProfiling")
        toggleTab("user-tabs","Visitors map",!is.null(ACCESS.LOG))                    

        message("[SERVER] reconfiguring menu done.")        

    })
    
    ##-------------------------------------------------------------
    ## Session TimerModule
    ##-------------------------------------------------------------

    reset_timer <- function() {}
    run_timer <- function(run=TRUE) {}
  
    if( TIMEOUT > 0 ) {

        rv.timer <- reactiveValues(reset=0, run=FALSE)
        reset_timer <- function() {
            dbg("[SERVER] resetting timer")
            rv.timer$reset <- rv.timer$reset + 1
        }
        run_timer <- function(run=TRUE) {
            dbg("[SERVER] run timer =",run)
            rv.timer$run <- run
        }
        WARN_BEFORE = round(TIMEOUT/6)

        message("[SERVER] Creating TimerModule...")
        message("[SERVER] TIMEOUT = ", TIMEOUT)
        message("[SERVER] WARN_BEFORE = ", WARN_BEFORE)                  
      
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
            message("[SERVER] timer = ",timer$timer())
            message("[SERVER] lapse_time = ",timer$lapse_time())        
        })

        observeEvent( timer$warn(), {
          message("[SERVER] timer$warn = ",timer$warn())
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
          message("[SERVER] success = ",success)          
          if(success==0) {
            message("[SERVER] logout after no referral!!!")
            shinyjs::runjs("logout()")    
          }
          if(success > 1) {
            message("[SERVER] resetting timer after referral!!!")
            showModal(modalDialog(
              HTML("<center><h4>Thanks!</h4>Your FREE session has been extended.</center>"),
              size = "s",
              easyClose = TRUE
            ))
            reset_timer()
          }

        })

      shiny::observeEvent( auth$logged(), {
        ## trigger on change of USER
        logged <- auth$logged()
        message("[SERVER] logged = ",logged)
        
        ##--------- start timer --------------
        if(TIMEOUT>0 && logged) {        
          message("[SERVER] starting session timer!!!")
          reset_timer()
          run_timer(TRUE)            
          
        } else {
          message("[SERVER] no timer!!!")            
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
        message("[SERVER] logged = ",logged)
        
        ##--------- force logout callback??? --------------
        if(opt$AUTHENTICATION!='firebase' && !logged) {
            ## Forcing logout ensures "clean" sessions. For firebase
            ## we allow sticky sessions.
            message("[SERVER] user not logged in? forcing logout() JS callback...")            
            shinyjs::runjs("logout()")    
        }
        
    })

    ## logout helper function
    logout.JScallback = "logout()"    
    if(opt$AUTHENTICATION=="shinyproxy") {
        logout.JScallback = "function(x){logout();quit();window.location.assign('/logout');}"
    }

    ## This will be called upon user logout after the logout() JS call
    observeEvent( input$userLogout, {
      message("[SERVER] >>>>>>>>> observe::input$userLogout reacted")
      reset_timer()
      run_timer(FALSE)
    })

    ## This code listens to the JS quit signal
    observeEvent( input$quit, {
        dbg("[SERVER:quit] !!!reacted!!!")
        dbg("[SERVER:quit] closing session... ")        
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
    message("[SERVER] server.init_time = ",server.init_time," ",attr(server.init_time,"units"))
    total.lapse_time <- round(Sys.time() - main.start_time,digits=4)
    message("[SERVER] total lapse time = ",total.lapse_time," ",attr(total.lapse_time,"units"))

}
