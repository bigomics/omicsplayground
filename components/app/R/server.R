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
    
    has.honcho <- Sys.getenv("HONCHO_TOKEN","")!="" &&
        !is.null(opt$HONCHO_URL) && opt$HONCHO_URL!=""
    if(1 && has.honcho) {
        sever::sever(sever_screen2(session$token), bg_color = "#000000") 
    } else {
        ## No honcho, no email....
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
  
    ## User board
    WelcomeBoard("welcome", auth=auth) 
    UserBoard("user", user=auth) -> env$user  
    
    ## Modules needed after dataset is loaded (deferred)
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
        ##pgx  <- env$load$inputData() 

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
    ## Session Timer (can we put it elsewhere?)
    ##-------------------------------------------------------------
    tm <- reactiveValues(timer=reactiveTimer(Inf), start=NULL)
    tm.warned <- FALSE
    shiny::observe({
        ## trigger on change of USER
        level <- auth$level()
        message("[SERVER] user LEVEL = ",level)
        logged <- auth$logged()
        message("[SERVER] logged = ",logged)
        
        ##--------- force logout callback??? --------------
        if(opt$AUTHENTICATION!='firebase' && !logged) {
            ## Forcing logout ensures "clean" sessions. For firebase
            ## we allow sticky sessions.
            message("[SERVER] user not logged in? forcing logout() JS callback...")            
            shinyjs::runjs("logout()")    
        }
        
        ##--------- start timer --------------
        ##if(tolower(level)=="free" && TIMEOUT>0 && logged) {
        if(TIMEOUT>0 && logged) {        
            message("[SERVER] starting session timer!!!")
            message("[SERVER] TIMEOUT = ", TIMEOUT)            
            tm$timer <- reactiveTimer(0.1*TIMEOUT*1000)  ## polling time
            tm$start <- Sys.time()
            tm.warned <<- FALSE
        } else {
            message("[SERVER] no timer!!!")            
            tm$timer <- reactiveTimer(Inf)
            tm$start <- NULL
            tm.warned <<- FALSE
        }       
    })
    
    ## Logout user after TIMEOUT
    shiny::observeEvent(tm$timer(), {
        message("[SERVER] timer reacted")
        if(is.null(tm$start) || TIMEOUT <=0 ) {
            message("[SERVER] timer is off")
            return()
        }
        if(!is.null(tm$start)) {
            secs.lapsed <- as.numeric(Sys.time() - tm$start, units='secs')            
            message("[SERVER] timer seconds lapsed = ",round(secs.lapsed,digits=2))
            if(secs.lapsed >= 0.80*TIMEOUT && !tm.warned) {
                message("[SERVER] timed out warning!!!")
                shinyalert::closeAlert()
                shinyalert::shinyalert(
                                title = "Warning!",
                                text = "Your FREE session is expiring soon."
                            )
                tm.warned <<- TRUE
            } else if(secs.lapsed >= TIMEOUT) {
                message("[SERVER] timed out!!!")
                shinyalert::closeAlert()
                js.cb = "function(x){logout();}"
                if(opt$AUTHENTICATION=="shinyproxy") {
                    js.cb = "function(x){logout();quit();window.location.assign('/logout');}"
                }
                showModal(
                    modalDialog(
                        title = "Session Expired",
                        p(
                            "Please enter three email addresses."
                        ),
                        fluidRow(
                            column(
                                4,
                                textInput(
                                    "name1",
                                    "Name"
                                ),
                                textInput(
                                    "email1",
                                    "Email"
                                )
                            ),
                            column(
                                4,
                                textInput(
                                    "name2",
                                    "Name"
                                ),
                                textInput(
                                    "email2",
                                    "Email"
                                )
                            ),
                            column(
                                4,
                                textInput(
                                    "name3",
                                    "Name"
                                ),
                                textInput(
                                    "email3",
                                    "Email"
                                )
                            )
                        ),
                        p(
                            class = "text-danger text-center",
                            id = "referral-global-error"
                        ),
                        footer = tagList(
                            actionButton(
                                "sendRefs",
                                "Send emails",
                                icon = icon("paper-plane")
                            ),
                            tags$a(
                                     class = "btn btn-danger",
                                     icon("times"),
                                     "Close",
                                     onClick = HTML(js.cb)
                                 )
                        )
                    )
                )
                tm$timer <- reactiveTimer(Inf)                
                tm$start <- NULL
                tm.warned <<- FALSE                
            }
        }
    })

    observeEvent(input$sendRefs, {
                                        # check inputs
        input_errors <- FALSE

                                        # check emails
        if(input$email1 == "") {
            session$sendCustomMessage(
                        "referral-input-error", 
                        list(
                            target = "email1",
                            message = "Missing email"
                        )
                    )
            input_errors <- TRUE
        }

        if(input$email2 == "") {
            session$sendCustomMessage(
                        "referral-input-error", 
                        list(
                            target = "email2",
                            message = "Missing email"
                        )
                    )
            input_errors <- TRUE
        }

        if(input$email3 == "") {
            session$sendCustomMessage(
                        "referral-input-error", 
                        list(
                            target = "email3",
                            message = "Missing email"
                        )
                    )
            input_errors <- TRUE
        }

                                        # check names
        if(input$name1 == "") {
            session$sendCustomMessage(
                        "referral-input-error", 
                        list(
                            target = "name1",
                            message = "Missing name"
                        )
                    )
            input_errors <- TRUE
        }
        
        if(input$name2 == "") {
            session$sendCustomMessage(
                        "referral-input-error", 
                        list(
                            target = "name2",
                            message = "Missing name"
                        )
                    )
            input_errors <- TRUE
        }

        if(input$name3 == "") {
            session$sendCustomMessage(
                        "referral-input-error", 
                        list(
                            target = "name3",
                            message = "Missing name"
                        )
                    )
            input_errors <- TRUE
        }

        emails <- trimws(
            c(
                input$email1,
                input$email2,
                input$email3
            )
        )

        if(!all(grepl("\\@", emails))) {
            session$sendCustomMessage(
                        "referral-global-error", 
                        list(
                            message = "Must enter valid email addresses"
                        )
                    )
            input_errors <- TRUE
        }

        if(length(unique(emails)) < 3) {
            session$sendCustomMessage(
                        "referral-global-error", 
                        list(
                            message = "Must enter different email addresses"
                        )
                    )
            input_errors <- TRUE
        }

        if(input_errors)
            return()

                                        # send emails
            body <- list(
                referrer = "The user",
                referrals = list(
                    list(
                        name = input$name1,
                        email = input$email1
                    ),
                    list(
                        name = input$name2,
                        email = input$email2
                    ),
                    list(
                        name = input$name3,
                        email = input$email3
                    )
                )
            )
            token <- Sys.getenv("HONCHO_TOKEN", "")
            uri <- sprintf("%s/referral?token=%s", opt$HONCHO_URL, token)
            response <- httr::POST(
                                  uri,
                                  body = body,
                                  encode = "json"
                              )

                                        # check response
            content <- httr::content(response)
            all_good <- lapply(content, function(ref) {
                return(ref$success)
            }) %>% 
                unlist() %>% 
                all()
            
            if(!all_good) {
                session$sendCustomMessage(
                            "referral-global-error", 
                            list(
                                message = "One or more of these email address was erroneous"
                            )
                        )
                return()
            }
                                        # remove modal
            removeModal()
    })
    
    ##-------------------------------------------------------------
    ## Session logout functions
    ##-------------------------------------------------------------
    
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
