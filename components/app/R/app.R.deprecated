##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2021 BigOmics Analytics Sagl. All rights reserved.
##

#########################################################################
##                                                                     ##
##              Main application for Omics Playground                  ##
##                                                                     ##
#########################################################################

message("\n\n")
message("###############################################################")
message("##################### OMICS PLAYGROUND ########################")
message("###############################################################")
message("\n")

main.start_time <- Sys.time()

WORKDIR = getwd()
message(">>>>> working directory = ",WORKDIR)
message(">>>>> LOADING INITIAL LIBS")

## some libraries that we often need and load fast
library(shiny)
library(shinyBS)
library(pryr)
library(grid)
library(ggplot2)

message("***********************************************")
message("***** RUNTIME ENVIRONMENT VARIABLES ***********")
message("***********************************************")

if(file.exists("Renviron.site")) {
    message("Loading local Renviron.site")    
    readRenviron("Renviron.site")
}

envcat <- function(var) message(var," = ",Sys.getenv(var))
envcat("SHINYPROXY_USERNAME")
envcat("SHINYPROXY_USERGROUPS")
envcat("PLAYGROUND_AUTHENTICATION")
envcat("PLAYGROUND_USERID")
envcat("PLAYGROUND_EXPIRY")
envcat("PLAYGROUND_QUOTA")
envcat("PLAYGROUND_LEVEL")
envcat("PLAYGROUND_HELLO")
envcat("OMICS_GOOGLE_TAG")

## --------------------------------------------------------------------
## -------------------------- INIT ------------------------------------
## --------------------------------------------------------------------

message("\n***********************************************")
message("*********** SETTING GLOBAL VARIABLES **********")
message("***********************************************")

source("global.R")  ## global variable
message("OPG =",OPG)
message("RDIR =",RDIR)
message("FILES =",FILES)
message("FILESX =",FILESX)
message("PGX.DIR =",PGX.DIR)
message("SHINYPROXY = ",SHINYPROXY)

DEV = (DEV && dir.exists("modulesx"))
##DEV = FALSE
if(DEV) {
    message('!!!!!!!!!! DEVELOPER MODE !!!!!!!!!!!!!!')
}

message("\n************************************************")
message("**************** READ FUNCTIONS ****************")
message("************************************************")

source(file.path(RDIR,"pgx-include.R"))    ## lots of libraries and source()
source(file.path(RDIR,"pgx-functions.R")) ## functions...
source(file.path(RDIR,"pgx-files.R"))     ## file functions
source(file.path(RDIR,"pgx-init.R"))
source(file.path(RDIR,"auth.R"))
source("app-init.R")

if(0) {
    ## pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)
    load("../data/geiger2016-arginine.pgx")
    load("../data/GSE10846-dlbcl-nc.pgx")
    load("../data/bojkova2020-sarscov2-RC2.pgx")
    load("../data/gtex-aging-n40svaNnm.pgx")
    load("../data/axel-test3.pgx")
    ngs = pgx.initialize(ngs)
}

message("\n************************************************")
message("************* parsing OPTIONS file *************")
message("************************************************")

if(!file.exists("OPTIONS")) stop("FATAL ERROR: cannot find OPTIONS file")
opt <- pgx.readOptions(file="OPTIONS")

## Check and set authentication method
if(Sys.getenv("PLAYGROUND_AUTHENTICATION")!="") {
    auth <- Sys.getenv("PLAYGROUND_AUTHENTICATION")
    message("[ENV] overriding PLAYGROUND_AUTHENTICATION = ",auth)
    opt$AUTHENTICATION = auth
}
if(1 && opt$AUTHENTICATION=="shinyproxy" && !in.shinyproxy()) {
    Sys.setenv("SHINYPROXY_USERNAME"="Test Person")  ## only for testing!!
}
if(1 && opt$AUTHENTICATION=="firebase" && !file.exists("firebase.rds")) {
    message("[ENV] WARNING: Missing firebase.rds file!!! reverting authentication to 'none'")
    opt$AUTHENTICATION = "none"
    ## opt$ENABLE_USERDIR = FALSE
    ## stop("[MAIN] FATAL Missing firebase.rds file")
}

## copy to global.R environment
WATERMARK <<- opt$WATERMARK
TIMEOUT   <<- as.integer(opt$TIMEOUT)  ## in seconds

## show options
message("\n",paste(paste(names(opt),"\t= ",sapply(opt,paste,collapse=" ")),collapse="\n"),"\n")

http.resp <- getFromNamespace("httpResponse", "shiny")

logHandler <- function(http.req){

    dbg("[MAIN.logHandler] >>>>> called! <<<<<")
    ##dbg("[MAIN.logHandler] names(http.req) = ",sort(names(http.req)))
    dbg("[MAIN.logHandler] http.req$PATH_INFO = ",http.req$PATH_INFO)

    if(!http.req$PATH_INFO == "/log") {
        return()
    }

    query <- shiny::parseQueryString(http.req$QUERY_STRING)
    dbg("[MAIN.logHandler] names(query) = ",names(query))
    dbg("[MAIN.logHandler] query$msg = ",query$msg)

    if(is.null(query$msg)) {
        dbg("[MAIN.logHandler] msg is NULL!")
        return(http.resp(400L, "application/json", jsonlite::toJSON(FALSE)))
    }

    if(query$msg == "") {
        dbg("[MAIN.logHandler] msg is empty!")
        return(http.resp(400L, "application/json", jsonlite::toJSON(FALSE)))
    }

    token <- Sys.getenv("HONCHO_TOKEN", "")
    if(token == "") {
        dbg("[MAIN.logHandler] missing HONCHO_TOKEN!")
        return(http.resp(403L, "application/json", jsonlite::toJSON(FALSE)))
    }

    uri <- sprintf("%s/log?token=%s", opt$HONCHO_URL, token)

    ## get the correct log file
    log.file = NULL
    the.log <- "Could not find log file!"
    id <- query$msg  ## use session id as message

    log.dirs <- "~/ShinyApps/log/*log /var/log/shiny-server/*log"
    suppressWarnings( log.file <- system(paste("grep -l -s",id,log.dirs),intern=TRUE) )
    log.file
    log.file <- tail(log.file,1)  ## take newest???
    log.file

    if(length(log.file)==0) {
        dbg("[MAIN.logHandler] could not resolve log file for session ID = ",id)
        return(http.resp(403L, "application/json", jsonlite::toJSON(FALSE)))
    }

    dbg("[logHandler] reading log.file = ",log.file)
    if(!is.null(log.file)) {
        ##the.log <- readr::read_file(log.file)
        ## truncate the log file
        the.log <- paste(system(paste("grep -B100 -A99999",id,log.file),intern=TRUE),collapse='\n')
    }

    dbg("[logHandler] sending log file... ")
    httr::POST(
        uri,
        body = list(
            msg = query$msg,
            log = "The log!",
            filename = "the_log.log"
        ),
        encode = "json"
    )

    http.resp(400L, "application/json", jsonlite::toJSON(TRUE))
}

run_application <- function(ui, server, ...){
    ## get handler
    handlerManager <- getFromNamespace("handlerManager", "shiny")

    ## add handler
    handlerManager$removeHandler("/log")
    handlerManager$addHandler(logHandler, "/log")

    shiny::shinyApp(ui, server, ...)
}

## --------------------------------------------------------------------
## ----------------- READ MODULES/BOARDS ------------------------------
## --------------------------------------------------------------------

modules <- dir("modules", pattern="Module.R$")
modules
for(m in modules) {
    message("[MAIN] loading module ",m)
    source(paste0("modules/",m))
}

BOARDS <- c("load","view","clust","expr","enrich","isect","func",
            "word","drug","sig","scell","cor","bio","cmap","ftmap",
            "wgcna", "tcga","multi","system","qa","corsa","comp","user")
if(is.null(opt$BOARDS_ENABLED))  opt$BOARDS_ENABLED = BOARDS
if(is.null(opt$BOARDS_DISABLED)) opt$BOARDS_DISABLED = NA

ENABLED  <- array(BOARDS %in% opt$BOARDS_ENABLED, dimnames=list(BOARDS))
DISABLED <- array(BOARDS %in% opt$BOARDS_DISABLED, dimnames=list(BOARDS))
ENABLED  <- ENABLED & !DISABLED
ENABLED

ENABLED[c("system","multi","corsa")] <- FALSE
if(0 && DEV && dir.exists("modulesx")) {
    ## Very early development modules/boards (ALWAYS SHOW FOR DEV)
    ##
    xboards <- dir("modulesx", pattern="Board.R$")
    xboards
    m=xboards[1]
    for(m in xboards) {
        message("[MAIN] loading DEVELOPMENT modules ",m)
        source(paste0("modulesx/",m))
    }
    ENABLED[] <- TRUE  ## enable all modules
    boards <- unique(c(boards, xboards))
}
ENABLED

## disable connectivity map if we have no signature database folder
has.sigdb <- length(dir(SIGDB.DIR,pattern="sigdb.*h5"))>0
has.sigdb
if(has.sigdb==FALSE) ENABLED["cmap"] <- FALSE

MAINTABS = c("DataView","Clustering","Expression","Enrichment",
             "Signature","CellProfiling","DEV")

main.init_time <- round(Sys.time() - main.start_time,digits=4)
main.init_time
message("[MAIN] main init time = ",main.init_time," ",attr(main.init_time,"units"))

## --------------------------------------------------------------------
## --------------------------- SERVER ---------------------------------
## --------------------------------------------------------------------

server = function(input, output, session) {

    message("\n========================================================")
    message("===================== SERVER ===========================")
    message("========================================================\n")
    dbg("[SERVER] 0: getwd = ",getwd())
    dbg("[SERVER] 0: HONCHO_URL = ",opt$HONCHO_URL)
    dbg("[SERVER] 0: SESSION = ",session$token)
    ##dbg("[SERVER] 0: names(session) = ",names(session))

    ## Logging of input/output events
    log.path <- "../logs/"
    log.path <- file.path(OPG,"logs")
    shinylogs::track_usage(storage_mode = shinylogs::store_rds(path = log.path))

    has.honcho <- Sys.getenv("HONCHO_TOKEN","")!="" &&
        !is.null(opt$HONCHO_URL) && opt$HONCHO_URL!=""
    if(1 && has.honcho) {
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

            ## not yet...
            ##if(!is.null(query[['pgxdir']])) {
            ##    pgx_dir <- query[['pgxdir']]
            ##}

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
    env[["user"]] <- shiny::callModule(UserBoard, "user", env)
    ##shinyjs::runjs("logout()")

    ## Modules needed after dataset is loaded (deferred)
    modules_loaded <- FALSE
    observeEvent( env[["load"]]$loaded(), {

        env.loaded <- env[["load"]]$loaded()
        message("[SERVER:env.loaded] env.loaded = ",env.loaded)

        ## on.exit({
        ##     message("[SERVER:env.loaded] on.exit::removing Modal")
        ##     Sys.sleep(4*modules_loaded)
        ##     shiny::removeModal()  ## remove modal from LoadingBoard
        ## })

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
            if(ENABLED["view"])   env[["view"]]   <- shiny::callModule( DataViewBoard, "view", env)
            if(ENABLED["clust"])  env[["clust"]]  <- shiny::callModule( ClusteringBoard, "clust", env)
            if(ENABLED["ftmap"])  env[["ftmap"]]  <- shiny::callModule( FeatureMapBoard, "ftmap", env)
            shiny::incProgress(0.2)
            if(ENABLED["expr"])   env[["expr"]]   <- shiny::callModule( ExpressionBoard, "expr", env)
            if(ENABLED["enrich"]) env[["enrich"]] <- shiny::callModule( EnrichmentBoard, "enrich", env)
            if(ENABLED["func"])   env[["func"]]   <- shiny::callModule( FunctionalBoard, "func", env)
            if(ENABLED["word"])   env[["word"]]   <- shiny::callModule( WordCloudBoard, "word", env)
            shiny::incProgress(0.4)
            if(ENABLED["drug"])   env[["drug"]]   <- shiny::callModule( DrugConnectivityBoard, "drug", env)
            if(ENABLED["isect"])  env[["isect"]]  <- shiny::callModule( IntersectionBoard, "isect", env)
            if(ENABLED["sig"])    env[["sig"]]    <- shiny::callModule( SignatureBoard, "sig", env)
            if(ENABLED["cor"])    env[["cor"]]    <- shiny::callModule( CorrelationBoard, "cor", env)
            shiny::incProgress(0.6)
            if(ENABLED["bio"])    env[["bio"]]    <- shiny::callModule( BiomarkerBoard, "bio", env)
            if(ENABLED["cmap"])   env[["cmap"]]   <- shiny::callModule( ConnectivityBoard, "cmap", env)
            if(ENABLED["scell"])  env[["scell"]]  <- shiny::callModule( SingleCellBoard, "scell", env)
            shiny::incProgress(0.8)
            if(ENABLED["tcga"])   env[["tcga"]]   <- shiny::callModule( TcgaBoard, "tcga", env)
            if(ENABLED["wgcna"])  env[["wgcna"]]  <- shiny::callModule( WgcnaBoard, "wgcna", env)
            if(ENABLED["comp"])   env[["comp"]]   <- shiny::callModule( CompareBoard, "comp", env)
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
                js.cb = "logout()"                
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
            referrer = env[["load"]]$auth$name(),
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
        if(0) {
            ## Return non-zero value so docker swarm can catch and restart
            ## the container upon on-failure
            dbg("[SERVER:quit] force stopping App... ")
            stopApp(99)
        }
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



    ## log(NULL)  ## force crash!!


}

## --------------------------------------------------------------------
## ------------------------------ UI ----------------------------------
## --------------------------------------------------------------------

TABVIEWS <- list(
    "load"   = tabView("Home",LoadingInputs("load"),LoadingUI("load")),
    "view"   = tabView("DataView",DataViewInputs("view"),DataViewUI("view")),
    "clust"  = tabView("Cluster samples",ClusteringInputs("clust"),ClusteringUI("clust")),
    "ftmap"  = tabView("Cluster features",FeatureMapInputs("ftmap"),FeatureMapUI("ftmap")),
    "wgcna"  = tabView("WGCNA (beta)",WgcnaInputs("wgcna"),WgcnaUI("wgcna")),
    "expr"   = tabView("Differential expression",ExpressionInputs("expr"),ExpressionUI("expr")),
    "cor"    = tabView("Correlation analysis", CorrelationInputs("cor"), CorrelationUI("cor")),
    "enrich" = tabView("Geneset enrichment",EnrichmentInputs("enrich"), EnrichmentUI("enrich")),
    "func"   = tabView("Pathway analysis", FunctionalInputs("func"), FunctionalUI("func")),
    "word"   = tabView("Word cloud", WordCloudInputs("word"), WordCloudUI("word")),
    "drug"   = tabView("Drug connectivity", DrugConnectivityInputs("drug"), DrugConnectivityUI("drug")),
    "isect"  = tabView("Compare signatures", IntersectionInputs("isect"), IntersectionUI("isect")),
    "sig"    = tabView("Test signatures", SignatureInputs("sig"), SignatureUI("sig")),
    "bio"    = tabView("Find biomarkers", BiomarkerInputs("bio"), BiomarkerUI("bio")),
    "cmap"   = tabView("Find similar experiments", ConnectivityInputs("cmap"), ConnectivityUI("cmap")),
    "scell"  = tabView("CellProfiling", SingleCellInputs("scell"), SingleCellUI("scell")),
    "tcga"   = tabView("TCGA survival (beta)", TcgaInputs("tcga"), TcgaUI("tcga")),
    "comp"   = tabView("Compare datasets (beta)", CompareInputs("comp"), CompareUI("comp"))
)

if(DEV) {
    if(ENABLED["corsa"]) TABVIEWS$corsa = tabView("CORSA (dev)",CorsaInputs("corsa"),
                                                  CorsaUI("corsa"))
    if(ENABLED["system"]) TABVIEWS$system = tabView("Systems analysis (dev)",
                                                    SystemInputs("system"),SystemUI("system"))
    if(ENABLED["multi"]) TABVIEWS$multi = tabView("Multi-level (dev)", MultiLevelInputs("multi"),
                                                  MultiLevelUI("multi"))
}

names(TABVIEWS)
TABVIEWS <- TABVIEWS[names(TABVIEWS) %in% names(which(ENABLED))]
names(TABVIEWS)

#-------------------------------------------------------
## Build USERMENU
#-------------------------------------------------------
user.tab <-  tabView(title = "Settings", id="user", UserInputs("user"), UserUI("user"))
##title = shiny::HTML("<span class='label label-info' id='authentication-user'></span>"),
logout.tab  <- shiny::tabPanel(shiny::HTML("<a onClick='logout()' id='authentication-logout'>Logout</a>"))

## conditionally add if firebase authentication is enabled
stop.tab    <- shiny::tabPanel(shiny::HTML("<a onClick='logout();quit();'>Quit</a>"))
if(opt$AUTHENTICATION == "shinyproxy") {
    ## For ShinyProxy we need to redirect to /logout for clean session
    ## logout. Then we need a redirect to the /login page.
    logout.tab  <- shiny::tabPanel(shiny::HTML("<a href='/login' onClick='shinyproxy_logout();' id='authentication-logout'>Logout</a>"))
}

upgrade.tab <- NULL
if(opt$AUTHENTICATION == "firebase") {
    upgrade.tab <- shiny::tabPanel(shiny::HTML("<a onClick='show_plans()' style='font-weight:bold;color:#2a9d8f;cursor:pointer;' id='authentication-upgrade'>Upgrade</a>"))
}

user.menu <- shiny::navbarMenu(
    ##title="User",
    title=icon("user-circle","fa"),
    user.tab,
    upgrade.tab,
    "----",
    shiny::tabPanel(title=shiny::HTML("<a href='https://omicsplayground.readthedocs.io' target='_blank'>Documentation</a>")),
    shiny::tabPanel(title=shiny::HTML("<a href='https://www.youtube.com/watch?v=_Q2LJmb2ihU&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-' target='_blank'>Video tutorials</a>")),
    shiny::tabPanel(title=shiny::HTML("<a href='https://groups.google.com/d/forum/omicsplayground' target='_blank'>Community Forum</a>")),
    shiny::tabPanel(title=shiny::HTML("<a href='https://github.com/bigomics/omicsplayground' target='_blank'>GitHub</a>")),
    "----",
    logout.tab,
    stop.tab
)

createUI <- function(tabs)
{
    message("\n======================================================")
    message("======================= UI ===========================")
    message("======================================================\n")

    version <- scan("../VERSION", character())[1]
    TITLE = paste(opt$TITLE,version)
    LOGO = shiny::div(shiny::img(src="bigomics-logo-white-48px.png", height="48px"),
               TITLE, id="navbar-logo", style="margin-top:-13px;")
    title = shiny::tagList(LOGO)
    windowTitle = TITLE
    theme = shinythemes::shinytheme("cerulean")
    id = "maintabs"

    ## Add Google Tag manager header code    
    gtag <- NULL
    if(Sys.getenv("OMICS_GOOGLE_TAG")!="") {
        gtag.html <- htmltools::includeHTML("www/google-tags.html")
        gtag.html <- sub("GTM-0000000",Sys.getenv("OMICS_GOOGLE_TAG"),gtag.html)
        gtag <- shiny::tags$head(gtag.html)
    }

    header = shiny::tagList(
        shiny::tags$head(shiny::tags$script(src="temp.js")),
        shiny::tags$head(shiny::tags$script(src="bigomics-extra.js")),  ## chatra,clarity
        gtag,   ## Google Tags???
        shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "styles.min.css")),
        shiny::tags$head(shiny::tags$link(rel="shortcut icon", href="favicon.ico")),
        shinyjs::useShinyjs(),
        sever::useSever(),
        shinylogs::use_tracking(),
        ##shinyalert::useShinyalert(),  # Set up shinyalert
        firebase::useFirebase(firestore = TRUE),
        ##TAGS.JSSCRIPT,  ## window size
        shiny::tags$script(async=NA, src="https://platform.twitter.com/widgets.js"),
        shiny::div(shiny::textOutput("current_dataset"), class='current-data'),
        shiny::div(class='label label-info current-user',id='authentication-user')
        ##QuestionBoard_UI("qa")
    )
    names(header) <- NULL

    footer.gif = shiny::tagList(
        shinybusy::busy_start_up(
            text = "\nPrepping your Omics Playground...", mode = "auto",
            background="#2780e3", color="#ffffff",
            ##loader = shiny::img(src=base64enc::dataURI(file="www/ready.png"))
            loader = shiny::img(src=base64enc::dataURI(file="www/monster-hi.png"))
        )
    )
    footer = footer.gif

    ##-------------------------------------
    ## create TAB list
    ##-------------------------------------
    createNavbarMenu <- function(title, tabs, icon=NULL) {
        tablist <- TABVIEWS[tabs]
        names(tablist) <- NULL
        do.call( navbarMenu, c(tablist, title=title, icon=icon) )
    }
    ##tablist <- TABVIEWS[tabs]
    tablist <- list()
    i=1
    for(i in 1:length(tabs)) {
        itab <- tabs[[i]]
        itab <- itab[which(ENABLED[itab])] ## only enabled
        if(length(itab)>1) {
            message("[MAIN] creating menu items for: ",paste(itab,collapse=" "))
            m <- createNavbarMenu( names(tabs)[i], itab )
            tablist[[i]] <- m
        } else if(length(itab)==1) {
            message("[MAIN] creating menu item for: ",itab)
            tablist[[i]] <- TABVIEWS[[itab]]
        } else {

        }
    }
    tablist <- tablist[!sapply(tablist,is.null)]

    ## add user menu (profile, help + logout)
    tablist[["usermenu"]] <- user.menu

    ##-------------------------------------
    ## create navbarPage
    ##-------------------------------------
    selected = "Home"
    names(tablist) <- NULL
    do.call( navbarPage, c(tablist,
                           title=title, id=id,
                           selected=selected,
                           windowTitle = windowTitle,
                           header = shiny::tagList(header),
                           footer = shiny::tagList(footer),
                           theme = theme))
}

tabs = list(
    "Home" = c("load"),
    "DataView" = "view",
    "Clustering" = c("clust","ftmap","wgcna"),
    "Expression" = c("expr","cor"),
    "Enrichment" = c("enrich","func","word","drug"),
    "Signature" = c("isect","sig","bio","cmap","comp","tcga"),
    "CellProfiling" = "scell",
    "DEV" = c("corsa","system","multi")
)


## Add Google Tag manager body code
gtag2 <- NULL
if(Sys.getenv("OMICS_GOOGLE_TAG")!="") {
    gtag2 <- htmltools::includeHTML("www/google-tags-noscript.html")
    gtag2 <- sub("GTM-0000000",Sys.getenv("OMICS_GOOGLE_TAG"),gtag2)
}

ui = tagList(
    gtag2,
    createUI(tabs)
)


## --------------------------------------------------------------------
## ------------------------------ RUN ---------------------------------
## --------------------------------------------------------------------


onStop(function() {
    message("*************** onStop:: doing application cleanup ****************")
    message("[APP] App died... ")

    ## Fill me...
})

onStart.FUN <- function() {
    message("*************** onStart:: preparing application ****************")
    message("[APP] App start!")
}

# shiny::shinyApp(ui, server, onStart=onStart.FUN)
run_application(ui, server)

##pkgs <- c( sessionInfo()[["basePkgs"]], names(sessionInfo()[["otherPkgs"]]),
##          names(sessionInfo()[["loadedOnly"]]) )

## --------------------------------------------------------------------
## ------------------------------ EOF----------------------------------
## --------------------------------------------------------------------
