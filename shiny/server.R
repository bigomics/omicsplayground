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
    
    ## Logging of input/output events
    shinylogs::track_usage(storage_mode = shinylogs::store_rds(path = "../logs/"))
    
    has.honcho <- Sys.getenv("HONCHO_TOKEN","")!="" &&
        !is.null(opt$HONCHO_URL) && opt$HONCHO_URL!=""
    if(1 && has.honcho) {
        ##sever::sever(sever_screen(), bg_color = "#000000") ## lightblue=2780e3
        sever::sever(sever_screen2(session$token), bg_color = "#000000") 
    } else {
        ## No honcho, no email....
        sever::sever(sever_screen0(), bg_color = "#000000") ## lightblue=2780e3
    }
    
    setwd(WORKDIR)  ## for some reason it can change!!
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
    if(ALLOW_URL_QUERYSTRING) {
        
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
                updateTabsetPanel(session, "load-tabs", selected = "Upload data")
                updateTextAreaInput(session, "load-upload_panel-compute-upload_description",
                                    value = "CSV FILE DESCRIPTION")                
            }
            
        })

    }
    dbg("[SERVER:parseQueryString] pgx_dir = ",pgx_dir)
    
    ##-------------------------------------------------------------
    ## Call modules
    ##-------------------------------------------------------------
    env <- list()  ## communication "environment"
    
    ## Modules needed from the start
        
    env[["load"]] <- shiny::callModule(
                                LoadingBoard, "load",
                                pgx_dir = pgx_dir,
                                limits = limits,
                                enable_userdir = opt$ENABLE_USERDIR,                                
                                authentication = authentication,
                                enable_upload = opt$ENABLE_UPLOAD,
                                enable_delete = opt$ENABLE_DELETE,                                 
                                enable_save = opt$ENABLE_SAVE
                                )   
    env[["user"]] <- shiny::callModule(UserBoard, "user", user = env[["load"]][["auth"]])
    ##shinyjs::runjs("logout()")    

    ## Modules needed after dataset is loaded (deferred)
    modules_loaded <- FALSE
    observeEvent( env[["load"]]$loaded(), {

        env.loaded <- env[["load"]]$loaded()
        message("[SERVER:env.loaded] env.loaded = ",env.loaded)    

        if(env[["load"]]$loaded()==0){
            message("[SERVER:env.loaded] env.loaded = FALSE")                                    
            return(NULL)
        }
        
        if(modules_loaded) {
            message("[SERVER:env.loaded] modules already loaded!")
            Sys.sleep(4)
            shiny::removeModal()  ## remove modal from LoadingBoard
            return(NULL)
        }
        modules_loaded <<- TRUE

        ## load other modules if
        message("[SERVER:env.loaded] --------- calling shiny modules ----------")
        shiny::withProgress(message="initializing modules ...", value=0, {
            if(ENABLED["view"])   env[["view"]]   <- shiny::callModule( DataViewBoard, "view", inputData = env[["load"]][["inputData"]])
            if(ENABLED["clust"])  env[["clust"]]  <- shiny::callModule( ClusteringBoard, "clust", inputData = env[["load"]][["inputData"]])
            if(ENABLED["ftmap"])  env[["ftmap"]]  <- shiny::callModule( FeatureMapBoard, "ftmap", inputData <- env[["load"]][["inputData"]])    
            shiny::incProgress(0.2)
            if(ENABLED["expr"])   env[["expr"]]   <- shiny::callModule( ExpressionBoard, "expr", inputData <- env[["load"]][["inputData"]])
            if(ENABLED["enrich"]) env[["enrich"]] <- shiny::callModule( EnrichmentBoard,
                                                                        "enrich",
                                                                        inputData = env[["load"]][["inputData"]],
                                                                        selected_gxmethods = env[["expr"]][["selected_gxmethods"]])
            if(ENABLED["func"])   env[["func"]]   <- shiny::callModule( FunctionalBoard, "func",
                                                                        inputData = env[["load"]][["inputData"]],
                                                                        selected_gxmethods = env[["expr"]][["selected_gxmethods"]],
                                                                        selected_gsetmethods = env[["enrich"]][["selected_gsetmethods"]])
            if(ENABLED["word"])   env[["word"]]   <- shiny::callModule( WordCloudBoard, "word", inputData = env[["load"]][["inputData"]],
                                                                        selected_gxmethods = env[["expr"]][["selected_gxmethods"]],
                                                                        selected_gsetmethods = env[["enrich"]][["selected_gsetmethods"]])
            shiny::incProgress(0.4)
            if(ENABLED["drug"])   env[["drug"]]   <- shiny::callModule( DrugConnectivityBoard, "drug", inputData = env[["load"]][["inputData"]])
            if(ENABLED["isect"])  env[["isect"]]  <- shiny::callModule( IntersectionBoard, "isect", inputData = env[["load"]][["inputData"]],
                                                                        selected_gxmethods = env[["expr"]][["selected_gxmethods"]],
                                                                        selected_gsetmethods = env[["enrich"]][["selected_gsetmethods"]])
            if(ENABLED["sig"])    env[["sig"]]    <- shiny::callModule( SignatureBoard, "sig", 
                                                                        inputData = env[["load"]][["inputData"]],
                                                                        selected_gxmethods = env[["expr"]][["selected_gxmethods"]])
            if(ENABLED["cor"])    env[["cor"]]    <- shiny::callModule( CorrelationBoard, "cor", inputData = env[["load"]][["inputData"]])
            shiny::incProgress(0.6)            
            if(ENABLED["bio"])    env[["bio"]]    <- shiny::callModule( BiomarkerBoard, "bio", inputData = env[["load"]][["inputData"]])
            if(ENABLED["cmap"])   env[["cmap"]]   <- shiny::callModule( ConnectivityBoard,
                                                                        "cmap",
                                                                        inputData = env[["load"]][["inputData"]])
            if(ENABLED["scell"])  env[["scell"]]  <- shiny::callModule( SingleCellBoard, "scell", inputData <- env[["load"]][["inputData"]])
            shiny::incProgress(0.8)
            if(ENABLED["tcga"])   env[["tcga"]]   <- shiny::callModule( TcgaBoard,
                                                                        "tcga",
                                                                        inputData = env[["load"]][["inputData"]])
            if(ENABLED["wgcna"])  env[["wgcna"]]  <- shiny::callModule( WgcnaBoard, "wgcna", inputData = env[["load"]][["inputData"]])
            if(ENABLED["comp"])   env[["comp"]]   <- shiny::callModule( CompareBoard, "comp", inputData = env[["load"]][["inputData"]])
            if(DEV) {            
                if(ENABLED["corsa"])  env[["corsa"]]  <- shiny::callModule( CorsaBoard, "corsa", env)
                if(ENABLED["system"]) env[["system"]] <- shiny::callModule( SystemBoard, "system", env)
                if(ENABLED["multi"])  env[["multi"]]  <- shiny::callModule( MultiLevelBoard, "multi", env)
                env[["qa"]] <- shiny::callModule( QuestionBoard, "qa", lapse = -1)
            }
        })
        message("[SERVER:env.loaded] --------- done! ----------")
        ## remove modal from LoadingBoard
        shiny::removeModal()
    })
    
    ## message("[SERVER] all boards called:",paste(names(env),collapse=" "))
    message("[SERVER] boards enabled:",paste(names(which(ENABLED)),collapse=" "))
    
    output$current_user <- shiny::renderText({
        ## trigger on change of user
        user <- env[["load"]][["auth"]]$email()
        user
    })
   
    output$current_dataset <- shiny::renderText({
        ## trigger on change of dataset
        pgx <- env[["load"]][["inputData"]]()
        name <- gsub(".*\\/|[.]pgx$","",pgx$name)
        if(length(name)==0) name = "(no data)"
        name
    })
    
    ##--------------------------------------------------------------------------
    ## Dynamically hide/show certain sections depending on USERMODE/object
    ##--------------------------------------------------------------------------
    shiny::observe({
        
        ## trigger on change dataset
        pgx  <- env[["load"]]$inputData() 
        show.beta <- env[["user"]]$enable_beta()
        dbg("[SERVER] show.beta = ",show.beta)
        if(is.null(show.beta) || length(show.beta)==0) show.beta=FALSE
        
        ## hide all main tabs until we have an object
        if(is.null(pgx)) {
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
        tabRequire(pgx, "connectivity", "maintabs", "Find similar experiments")
        tabRequire(pgx, "drugs", "maintabs", "Drug connectivity")
        tabRequire(pgx, "wordcloud", "maintabs", "Word cloud")
        tabRequire(pgx, "deconv", "maintabs", "CellProfiling")
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
        auth <- env[["load"]][["auth"]]
        level <- auth$level()
        message("[SERVER] user LEVEL = ",level)
        logged <- auth$logged()
        message("[SERVER] logged = ",logged)
        
        ##--------- force logout callback??? --------------
        if(opt$AUTHENTICATION!='firebase' && !logged) {
            ## Forcing logout ensures "clean" sessions. For firebase
            ## we allow sticky sessions.
            message("[SERVER] not logged? forcing logout() JS callback...")            
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
