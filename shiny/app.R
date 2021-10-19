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

##Sys.setenv("SHINYPROXY_USERNAME"="Test Person")
main.start_time <- Sys.time()

message("***********************************************")
message("*********** LOADING INITIAL LIBS **************")
message("***********************************************")

## some libraries that we often need and load fast
library(shiny)
library(shinyBS)
library(pryr)
library(grid)

## we need all these datasets that actually aren't datasets
## and so cannot be imported by data() function...
##library(org.Hs.eg.db) ## better use require inside?
##library(org.Mm.eg.db) ## better use require inside?

message("***********************************************")
message("***** RUNTIME ENVIRONMENT VARIABLES ***********")
message("***********************************************")

envcat <- function(var) message(var," = ",Sys.getenv(var))
envcat("SHINYPROXY_USERNAME")
envcat("SHINYPROXY_USERGROUPS")
envcat("PLAYGROUND_AUTHENTICATION")
envcat("PLAYGROUND_USERID")
envcat("PLAYGROUND_EXPIRY")
envcat("PLAYGROUND_QUOTA")
envcat("PLAYGROUND_LEVEL")
envcat("PLAYGROUND_HELLO")

## --------------------------------------------------------------------
## -------------------------- INIT ------------------------------------
## --------------------------------------------------------------------

message("\n")
message("***********************************************")
message("*********** SETTING GLOBAL VARIABLES **********")
message("***********************************************")

source("global.R")  ## global variable
message("OPG =",OPG)
message("RDIR =",RDIR)
message("FILES =",FILES)
message("FILESX =",FILESX)
message("PGX.DIR =",PGX.DIR)
message("SHINYPROXY = ",SHINYPROXY)

src.local=TRUE  ## local or not-local, that's the question...
src.local=FALSE ## local or not-local, that's the question...
source(file.path(RDIR,"pgx-include.R"),local=src.local)    ## lots of libraries and source()
source(file.path(RDIR,"pgx-functions.R"), local=src.local) ## functions...
source(file.path(RDIR,"pgx-files.R"), local=src.local)     ## file functions
source(file.path(RDIR,"pgx-init.R"),local=src.local)
source(file.path(RDIR,"auth.R"),local=src.local)

message("\n")
message("************************************************")
message("************* parsing OPTIONS file *************")
message("************************************************")

if(!file.exists("OPTIONS")) stop("FATAL ERROR: cannot find OPTIONS file")
opt <- pgx.readOptions(file="OPTIONS")

## over-ride options (for DEBUGGING)
##opt$AUTHENTICATION = "none"
##opt$AUTHENTICATION = "password"
##opt$AUTHENTICATION = "register"
opt$AUTHENTICATION = "firebase"

if(Sys.getenv("PLAYGROUND_AUTHENTICATION")!="") {
    auth <- Sys.getenv("PLAYGROUND_AUTHENTICATION")
    message("[ENV] overriding PLAYGROUND_AUTHENTICATION = ",auth)
    opt$AUTHENTICATION = auth
}

## copy to global environment
SHOW_QUESTIONS = FALSE
AUTHENTICATION = opt$AUTHENTICATION
USER_MODE      = opt$USER_MODE
WATERMARK      = opt$WATERMARK

DEV = (DEV && dir.exists("modulesx")) 
##DEV = FALSE
if(DEV) {
    message('******************************************************')
    message('****************** DEVELOPER MODE ********************')
    message('******************************************************')    
}

## show options
message("\n",paste(paste(names(opt),"\t= ",sapply(opt,paste,collapse=" ")),collapse="\n"),"\n")

## --------------------------------------------------------------------
## ------------------------ READ FUNCTIONS ----------------------------
## --------------------------------------------------------------------

source("app-init.R", local=FALSE)
message('>>>> Initializing data folder')
##pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=TRUE)
pgx.initDatasetFolder(PGX.DIR, force=FALSE, verbose=TRUE)

if(0) {    
    ##PGX.DIR="../test/"
    pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)    
    load("../data/geiger2016-arginine.pgx")
    load("../data/GSE10846-dlbcl-nc.pgx")
    load("../data/GSE157905-lenvatinib-bc.pgx")
    load("../data/gtex-aging-n40svaNnm.pgx")
    ngs = pgx.initialize(ngs)
}

## --------------------------------------------------------------------
## ----------------- READ MODULES/BOARDS ------------------------------
## --------------------------------------------------------------------

source("modules/AuthenticationModule.R",local=src.local)
source("modules/ComputePgxModule.R",local=src.local)
source("modules/MakeContrastModule.R",local=src.local)
source("modules/NormalizeCountsModule.R",local=src.local)
source("modules/BatchCorrectModule.R",local=src.local)
source("modules/UploadModule.R",local=src.local)
##source("modules/UsersMapModule.R_")

BOARDS <- c("load","view","clust","expr","enrich","isect","func",
            "word","drug","sig","scell","cor","bio","cmap","ftmap",
            "wgcna", "tcga","multi","system","qa","corsa","comp")
if(is.null(opt$BOARDS_ENABLED))  opt$BOARDS_ENABLED = BOARDS
if(is.null(opt$BOARDS_DISABLED)) opt$BOARDS_DISABLED = NA

ENABLED  <- array(BOARDS %in% opt$BOARDS_ENABLED, dimnames=list(BOARDS))
DISABLED <- array(BOARDS %in% opt$BOARDS_DISABLED, dimnames=list(BOARDS))
ENABLED  <- ENABLED & !DISABLED
ENABLED

boards <- dir("boards", pattern="Board.R$")
boards
for(m in boards) {
    message("[MAIN] loading board ",m)
    source(paste0("boards/",m), local=src.local)
    ##source(paste0("boards/",m), local=FALSE)
}

##ENABLED[c("wgcna","system","multi")] <- FALSE
ENABLED[c("system","multi","corsa")] <- FALSE 
if(0 && DEV && dir.exists("modulesx")) {
    ## Very early development modules/boards (ALWAYS SHOW FOR DEV)
    ##
    xboards <- dir("modulesx", pattern="Board.R$")
    xboards
    m=xboards[1]
    for(m in xboards) {
        message("[MAIN] loading DEVELOPMENT modules ",m)
        source(paste0("modulesx/",m), local=src.local)
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
    
    message("[SERVER] USER_MODE = ", USER_MODE)
    server.start_time <- Sys.time()
    
    limits <- c("samples" = opt$MAX_SAMPLES,
                "comparisons" = opt$MAX_COMPARISONS,
                "genes" = opt$MAX_GENES,
                "genesets" = opt$MAX_GENESETS,
                "datasets" = opt$MAX_DATASETS)
    env <- list()  ## communication "environment"
    env[["load"]]  <- shiny::callModule(
                                 LoadingBoard, "load",
                                 limits = limits,
                                 authentication = AUTHENTICATION,
                                 enable_upload = opt$ENABLE_UPLOAD,
                                 enable_delete = opt$ENABLE_DELETE,                                 
                                 enable_save = opt$ENABLE_SAVE)   
    

    already_loaded <- FALSE
    observeEvent( env[["load"]]$loaded(), {

        env.loaded <- env[["load"]]$loaded()
        message("[SERVER] env.loaded = ",env.loaded)                                    
        
        if(!env[["load"]]$loaded()){
            message("[SERVER] env.loaded = FALSE")                                    
            return(NULL)
        }

        message("[SERVER] env.loaded : 2 ")

        on.exit({
            message("[SERVER] on.exit::removing Modal")                        
            shiny::removeModal()
        })

        message("[SERVER] env.loaded : 2 ")

        if(already_loaded) {
            message("[SERVER] modules already loaded!")            
            return(NULL)
        }

        message("[SERVER] env.loaded : 3 ")
        
        already_loaded <<- TRUE

        message("[SERVER] env.loaded : 4 ")
        
        ## load other modules if
        message("[SERVER] --------- calling shiny modules ----------")
        if(ENABLED["view"])   env[["view"]]   <- shiny::callModule( DataViewBoard, "view", env)
        if(ENABLED["clust"])  env[["clust"]]  <- shiny::callModule( ClusteringBoard, "clust", env)
        if(ENABLED["ftmap"])  env[["ftmap"]]  <- shiny::callModule( FeatureMapBoard, "ftmap", env)    
        if(ENABLED["expr"])   env[["expr"]]   <- shiny::callModule( ExpressionBoard, "expr", env)
        if(ENABLED["enrich"]) env[["enrich"]] <- shiny::callModule( EnrichmentBoard, "enrich", env)
        if(ENABLED["func"])   env[["func"]]   <- shiny::callModule( FunctionalBoard, "func", env)
        if(ENABLED["word"])   env[["word"]]   <- shiny::callModule( WordCloudBoard, "word", env)
        if(ENABLED["drug"])   env[["drug"]]   <- shiny::callModule( DrugConnectivityBoard, "drug", env)
        if(ENABLED["isect"])  env[["isect"]]  <- shiny::callModule( IntersectionBoard, "isect", env)
        if(ENABLED["sig"])    env[["sig"]]    <- shiny::callModule( SignatureBoard, "sig", env)
        if(ENABLED["cor"])    env[["cor"]]    <- shiny::callModule( CorrelationBoard, "cor", env)
        if(ENABLED["bio"])    env[["bio"]]    <- shiny::callModule( BiomarkerBoard, "bio", env)
        if(ENABLED["cmap"])   env[["cmap"]]   <- shiny::callModule( ConnectivityBoard, "cmap", env)
        if(ENABLED["scell"])  env[["scell"]]  <- shiny::callModule( SingleCellBoard, "scell", env)
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
    
    ## message("[SERVER] all boards called:",paste(names(env),collapse=" "))
    message("[SERVER] boards enabled:",paste(names(which(ENABLED)),collapse=" "))
    
    output$current_dataset <- shiny::renderText({
        pgx <- env[["load"]][["inputData"]]()
        name <- gsub(".*\\/|[.]pgx$","",pgx$name)
        if(length(name)==0) name = "(no data)"
        name
    })
    
    ## Dynamically hide/show certain sections depending on USERMODE/object
    shiny::observe({
        pgx <- env[["load"]][["inputData"]]() ## trigger on change dataset
        
        ## hide all main tabs until we have an object
        if(is.null(pgx)) {
            lapply(MAINTABS, function(m) shiny::hideTab("maintabs",m))
            if(!opt$ENABLE_UPLOAD)  shiny::hideTab("load-tabs","Upload data")
            if(is.null(ACCESS.LOG)) shiny::hideTab("load-tabs","Visitors map")            
            return(NULL)
        }

        message("[SERVER] dataset changed. reconfiguring menu...")
        ## show all main tabs
        lapply(MAINTABS, function(m) shiny::showTab("maintabs",m))
        
        if(USER_MODE == "basic") {
            shiny::hideTab("maintabs","CellProfiling")
            shiny::hideTab("clust-tabs2","Feature ranking")
            shiny::hideTab("expr-tabs1","Volcano (methods)")
            shiny::hideTab("expr-tabs2","FDR table")
            shiny::hideTab("enrich-tabs1","Volcano (methods)")
            shiny::hideTab("enrich-tabs2","FDR table")
            shiny::hideTab("cor-tabs","Functional")    ## too slow
            shiny::hideTab("cor-tabs","Differential")  ## too complex
        }
        
        ## shiny::hideTab("cor-tabs","Functional")       
        if(USER_MODE == "dev" || DEV) {
            shiny::showTab("maintabs","DEV")
            shiny::showTab("view-tabs","Resource info")
            shiny::showTab("scell-tabs1","CNV")  ## DEV only
            shiny::showTab("scell-tabs1","Monocle") ## DEV only
            shiny::showTab("cor-tabs","Functional")
        } else {
            shiny::hideTab("maintabs","DEV")
            shiny::hideTab("view-tabs","Resource info")
            shiny::hideTab("scell-tabs1","CNV")  ## DEV only
            shiny::hideTab("scell-tabs1","Monocle") ## DEV only
            shiny::hideTab("cor-tabs","Functional")            
        }
        
        ## Dynamically show upon availability in pgx object
        if(opt$ENABLE_UPLOAD) shiny::showTab("load-tabs","Upload data")            
        tabRequire(pgx, "connectivity", "maintabs", "Similar experiments")
        tabRequire(pgx, "drugs", "maintabs", "Drug connectivity")
        tabRequire(pgx, "wordcloud", "maintabs", "Word cloud")
        tabRequire(pgx, "deconv", "maintabs", "CellProfiling")
        fileRequire("tcga_matrix.h5", "maintabs", "TCGA survival (beta)")         
        if(!is.null(ACCESS.LOG)) shiny::showTab("load-tabs","Visitors map")                    

        message("[SERVER] reconfiguring menu done.")        
    })

    server.init_time <- round(Sys.time() - server.start_time, digits=4)    
    message("[SERVER] server.init_time = ",server.init_time," ",attr(server.init_time,"units"))
    total.lapse_time <- round(Sys.time() - main.start_time,digits=4)
    message("[SERVER] total lapse time = ",total.lapse_time," ",attr(total.lapse_time,"units"))

}

## --------------------------------------------------------------------
## ------------------------------ UI ----------------------------------
## --------------------------------------------------------------------

help.tabs <- shiny::navbarMenu(
    "Help",
    shiny::tabPanel(title=shiny::HTML("<a href='https://omicsplayground.readthedocs.io' target='_blank'>Documentation")),
    shiny::tabPanel(title=shiny::HTML("<a href='https://www.youtube.com/watch?v=_Q2LJmb2ihU&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-' target='_blank'>Video tutorials</a>")),
    shiny::tabPanel(title=shiny::HTML("<a href='https://github.com/bigomics/omicsplayground' target='_blank'>GitHub")),
    shiny::tabPanel(title=shiny::HTML("<a href='https://hub.docker.com/r/bigomics/omicsplayground' target='_blank'>Docker")),
    shiny::tabPanel(title=shiny::HTML("<a href='https://groups.google.com/d/forum/omicsplayground' target='_blank'>Google groups"))
)

TABVIEWS <- list(
    "load"   = tabView("Home",LoadingInputs("load"),LoadingUI("load")),
    "view"   = tabView("DataView",DataViewInputs("view"),DataViewUI("view")),
    "clust"  = tabView("Cluster samples",ClusteringInputs("clust"),ClusteringUI("clust")),
    "ftmap"  = tabView("Feature maps (beta)",FeatureMapInputs("ftmap"),FeatureMapUI("ftmap")),    
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
    "cmap"   = tabView("Similar experiments", ConnectivityInputs("cmap"), ConnectivityUI("cmap")),
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

logout.tab <- shiny::tabPanel(title=shiny::HTML("<a id='logout' href='/logout'>Logout"))


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
    ##selected = "Home"    
    header = shiny::tagList(
        shiny::tags$head(shiny::tags$script(src="temp.js")),
        shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "playground.css")),
        shiny::tags$head(shiny::tags$link(rel="shortcut icon", href="favicon.ico")),
        shinyjs::useShinyjs(),
        firebase::useFirebase(),
        TAGS.JSSCRIPT,
        shiny::tags$script(async=NA, src="https://platform.twitter.com/widgets.js"),
        shiny::div(
            shiny::textOutput("current_dataset"),
            class='current-data'
        )
        ##QuestionBoard_UI("qa")
    )
    names(header) <- NULL
    
    footer.gif = shiny::tagList(
        shinybusy::busy_start_up(
            text = "\nPrepping your Omics Playground...", mode = "auto",
            background="#2780e3", color="#ffffff",
            loader = shiny::img(src=base64enc::dataURI(file="www/ready.png"))
        )
    )
    ## if(runif(1) < 0.1)
    footer = footer.gif ## every now and then show easter egg..
    
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
   
    ## add help menu
    tablist[["helpmenu"]] <- help.tabs
    if(SHINYPROXY) {
        tablist[["logout"]] <- logout.tab
    }
    
    # conditionally add if firebase authentication is enabled
    if(opt$AUTHENTICATION == "firebase") {
        login.tabs <- shiny::navbarMenu(
            "User",
            shiny::tabPanel(
                title = shiny::HTML(
                    "<span class='label label-info' id='authentication-user'></span>"
                )
            ),
            shiny::tabPanel(
                title = shiny::HTML(
                    "<a onClick='logout()' id='authentication-logout'>Logout</a>"
                )
            ),
            shiny::tabPanel(
                title = shiny::HTML(
                    "<a onClick='upgrade()' style='font-weight:bold;color:darkgreen;' id='authentication-upgrade'>Upgrade</a>"
                )
            )
        )
        tablist[["login"]] <- login.tabs
    }
    
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

ui = createUI(tabs)
##ui = navbarPage("Hello World!")

## --------------------------------------------------------------------
## ------------------------------ RUN ---------------------------------
## --------------------------------------------------------------------

shiny::shinyApp(ui, server)

##pkgs <- c( sessionInfo()[["basePkgs"]], names(sessionInfo()[["otherPkgs"]]),
##          names(sessionInfo()[["loadedOnly"]]) )

## --------------------------------------------------------------------
## ------------------------------ EOF----------------------------------
## --------------------------------------------------------------------

