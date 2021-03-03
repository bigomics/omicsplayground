##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

#########################################################################
##                                                                     ##
##              Main application for Omics Playground                  ##
##                                                                     ##
#########################################################################

library(shiny)
library(shinyjs)
library(shinyWidgets)
library(waiter)
library(plotly)

message("\n\n")
message("###############################################################")
message("##################### OMICS PLAYGROUND ########################")
message("###############################################################")

message("\n")
message("************************************************")
message("********* RUNTIME ENVIRONMENT VARIABLES ********")
message("************************************************")

Sys.setlocale("LC_CTYPE","en_US.UTF-8") 
Sys.setlocale("LC_TIME","en_US.UTF-8")
##Sys.setlocale("LC_ALL", "C")  ## really??
##Sys.setenv("SHINYPROXY_USERNAME"="Test Person")

envcat <- function(var) message(var," = ",Sys.getenv(var))
envcat("SHINYPROXY_USERNAME")
envcat("SHINYPROXY_USERGROUPS")
envcat("PLAYGROUND_AUTHENTICATION")
envcat("PLAYGROUND_USERID")
envcat("PLAYGROUND_EXPIRY")
envcat("PLAYGROUND_LEVEL")
envcat("PLAYGROUND_HELLO")

## --------------------------------------------------------------------
## -------------------------- INIT ------------------------------------
## --------------------------------------------------------------------

message("\n")
message("***********************************************")
message("*********** SETTING GLOBAL VARIABLES **********")
message("***********************************************")

source("global.R")

message("OPG =",OPG)
message("RDIR =",RDIR)
message("FILES =",FILES)
message("FILESX =",FILESX)
message("PGX.DIR =",PGX.DIR)
message("DEBUG = ",DEBUG)
message("WATERMARK = ",WATERMARK)
message("SHINYPROXY = ",SHINYPROXY)

src.local=TRUE  ## need???
src.local=FALSE  ## need???
source(file.path(RDIR,"pgx-include.R"),local=src.local)  ## pass local vars
source(file.path(RDIR,"pgx-functions.R"), local=src.local)  ## pass local vars
source(file.path(RDIR,"pgx-files.R"), local=src.local)  ## pass local vars
source(file.path(RDIR,"pgx-init.R"),local=src.local)     ## pass local vars

message("\n")
message("************************************************")
message("************* parsing OPTIONS file *************")
message("************************************************")

options(shiny.maxRequestSize = 999*1024^2)  ##max 999Mb upload
if(!file.exists("OPTIONS")) stop("FATAL ERROR: cannot find OPTIONS file")
opt <- pgx.readOptions(file="OPTIONS")

## over-ride options (for DEBUGGING)
## opt$AUTHENTICATION = "none"
## opt$AUTHENTICATION = "password"
## opt$AUTHENTICATION = "register"
## opt$AUTHENTICATION = "firebase"


if(Sys.getenv("PLAYGROUND_AUTHENTICATION")!="") {
    auth <- Sys.getenv("PLAYGROUND_AUTHENTICATION")
    message("[ENV] overriding PLAYGROUND_AUTHENTICATION = ",auth)
    opt$AUTHENTICATION = auth
}

## copy to global environment
SHOW_QUESTIONS = FALSE
##WATERMARK      = opt$WATERMARK
##USER_MODE      = opt$USER_MODE
AUTHENTICATION = opt$AUTHENTICATION
DEV = (USER_MODE=="dev" && dir.exists("../../omicsplayground-dev"))

## show options
message("\n",paste(paste(names(opt),"\t= ",sapply(opt,paste,collapse=" ")),collapse="\n"),"\n")

## --------------------------------------------------------------------
## ------------------------ READ FUNCTIONS ----------------------------
## --------------------------------------------------------------------

source("init.R", local=FALSE)
source("modules/AuthenticationModule.R",local=src.local)
source("modules/ComputePgxModule.R",local=src.local)
source("modules/MakeContrastModule.R",local=src.local)
source("modules/NormalizeCountsModule.R",local=src.local)
source("modules/SuperBatchCorrectModule.R",local=src.local)
source("modules/UploadModule.R",local=src.local)
##source("modules/UsersMapModule.R_")

##pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)
pgx.initDatasetFolder(PGX.DIR, force=FALSE, verbose=1)


if(0) {    
    ##PGX.DIR="../test/"
    ##pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)    
    load("../data/geiger2016-arginine.pgx")
    load("../data/GSE10846-dlbcl-nc.pgx")
    ngs = pgx.initialize(ngs)
}



## --------------------------------------------------------------------
## ------------------------ READ BOARDS -------------------------------
## --------------------------------------------------------------------

BOARDS <- c("load","view","clust","expr","enrich","isect","func",
            "word","drug","sig","scell","cor","bio","cmap",
            "wgcna", "tcga","system","multi","qa")
if(is.null(opt$BOARDS_ENABLED)) opt$BOARDS_ENABLED = BOARDS
if(is.null(opt$BOARDS_DISABLED)) opt$BOARDS_DISABLED = NA
ENABLED  <- array(BOARDS %in% opt$BOARDS_ENABLED, dimnames=list(BOARDS))
DISABLED <- array(BOARDS %in% opt$BOARDS_DISABLED, dimnames=list(BOARDS))
ENABLED  <- ENABLED & !DISABLED
ENABLED

boards <- dir("boards", pattern="Board.R$")
for(m in boards) {
    message("[MAIN] loading board ",m)
    source(paste0("boards/",m), local=src.local)
    ##source(paste0("boards/",m), local=FALSE)
}

##ENABLED[c("wgcna","system","multi")] <- FALSE
ENABLED[c("system","multi")] <- FALSE
if(1 && DEV && dir.exists("../../omicsplayground-dev")) {
    xboards <- dir("../../omicsplayground-dev/modulesx", pattern="Board.R$")
    ##xboards <- dir("../../omicsplayground-dev/modulesx", pattern="TcgaBoard.R$")
    xboards
    m=xboards[1]
    for(m in xboards) {
        message("[MAIN] loading extra board ",m)
        source(paste0("../../omicsplayground-dev/modulesx/",m), local=src.local)
    }
    ##ENABLED[c("system","multi")] <- TRUE
    boards <- unique(c(boards, xboards))
}
ENABLED

## disable connectivity map if we have no signature database folder
has.sigdb <- length(dir(FILESX,pattern="sigdb.*h5")>0)
has.sigdb
if(has.sigdb==FALSE) ENABLED["cmap"] <- FALSE

MAINTABS = c("DataView","Clustering","Expression","Enrichment",
             "Signature","CellProfiling","Dev")

if(0) {
    save.image(file="../cache/image.RData")
    system.time( load(file="../cache/image.RData") )
}

## --------------------------------------------------------------------
## --------------------------- SERVER ---------------------------------
## --------------------------------------------------------------------

server = function(input, output, session) {
    
    message("\n========================================================")
    message("===================== SERVER ===========================")
    message("========================================================\n")

    message("[MAIN] calling boards...")
    message("[MAIN] USER_MODE = ", USER_MODE)
    
    library(firebase)
    firebase=firebase2=NULL
    if(AUTHENTICATION=="firebase") {
        firebase  <- FirebaseEmailPassword$new()
        firebase2 <- FirebaseSocial$new()
    }
    
    ## firebase <- NULL
    max.limits <- c("samples" = opt$MAX_SAMPLES,
                    "comparisons" = opt$MAX_COMPARISONS,
                    "genes" = opt$MAX_GENES,
                    "genesets" = opt$MAX_GENESETS)
    env <- list()  ## communication environment
    env[["load"]]   <- callModule(
        LoadingBoard, "load", max.limits = max.limits,
        authentication = AUTHENTICATION, enable_delete = opt$ENABLE_DELETE,
        enable_save = opt$ENABLE_SAVE,
        firebase=firebase, firebase2=firebase2)
    env[["view"]]   <- callModule( DataViewBoard, "view", env)
    env[["clust"]]  <- callModule( ClusteringBoard, "clust", env)
    env[["expr"]]   <- callModule( ExpressionBoard, "expr", env)
    env[["enrich"]] <- callModule( EnrichmentBoard, "enrich", env)
    if(ENABLED["func"])   env[["func"]]   <- callModule( FunctionalBoard, "func", env)
    if(ENABLED["word"])   env[["word"]]   <- callModule( WordCloudBoard, "word", env)
    if(ENABLED["drug"])   env[["drug"]]   <- callModule( DrugConnectivityBoard, "drug", env)
    if(ENABLED["isect"])  env[["isect"]]  <- callModule( IntersectionBoard, "isect", env)
    if(ENABLED["sig"])    env[["sig"]]    <- callModule( SignatureBoard, "sig", env)
    if(ENABLED["scell"])  env[["scell"]]  <- callModule( SingleCellBoard, "scell", env)
    if(ENABLED["cor"])    env[["cor"]]    <- callModule( CorrelationBoard, "cor", env)
    if(ENABLED["bio"])    env[["bio"]]    <- callModule( BiomarkerBoard, "bio", env)
    if(ENABLED["cmap"])   env[["cmap"]]   <- callModule( ConnectivityBoard, "cmap", env)
    if(ENABLED["tcga"])   env[["tcga"]]   <- callModule( TcgaBoard, "tcga", env)
    if(ENABLED["wgcna"])  env[["wgcna"]]  <- callModule( WgcnaBoard, "wgcna", env)
    if(1 && DEV) {
        if(ENABLED["system"]) env[["system"]] <- callModule( SystemBoard, "system", env)
        if(ENABLED["multi"])  env[["multi"]]  <- callModule( MultiLevelBoard, "multi", env)
    }
    env[["qa"]] <- callModule( QuestionBoard, "qa", lapse = -1)
    
    ## message("[MAIN] all boards called:",paste(names(env),collapse=" "))
    message("[MAIN] boards enabled:",paste(names(which(ENABLED)),collapse=" "))
    ## outputOptions(output, "clust", suspendWhenHidden=FALSE) ## important!!!
    
    output$current_dataset <- renderText({
        pgx <- env[["load"]][["inputData"]]()
        name <- gsub(".*\\/|[.]pgx$","",pgx$name)
        if(length(name)==0) name = "(no data)"
        name
    })
      
    ## Dynamicall hide/show certain sections depending on USERMODE/object
    observe({
        pgx <- env[["load"]][["inputData"]]() ## trigger on change dataset

        ## hide all main tabs until we have an object
        if(is.null(pgx)) {
            lapply(MAINTABS, function(m) hideTab("maintabs",m))
            if(!opt$ENABLE_UPLOAD)  hideTab("load-tabs","Upload data")
            if(is.null(ACCESS.LOG)) hideTab("load-tabs","Visitors map")            
            return(NULL)
        }

        message("[MAIN] dataset changed. reconfiguring menu...")
        ## show all main tabs
        lapply(MAINTABS, function(m) showTab("maintabs",m))
        
        if(USER_MODE == "basic") {
            hideTab("maintabs","CellProfiling")
            hideTab("enrich-tabs1","GeneMap")
            hideTab("clust-tabs2","Feature ranking")
            hideTab("expr-tabs1","Volcano (methods)")
            hideTab("expr-tabs2","FDR table")
            hideTab("enrich-tabs1","Volcano (methods)")
            hideTab("enrich-tabs2","FDR table")
            ## hideTab("cor-tabs","Functional")
        }

        ## hideTab("cor-tabs","Functional")
        
        if(USER_MODE == "dev" || DEV) {
            showTab("maintabs","Dev")
            showTab("view-tabs","Resource info")
            showTab("enrich-tabs1","GeneMap")
            showTab("scell-tabs1","CNV")  ## DEV only
            showTab("scell-tabs1","Monocle") ## DEV only
            showTab("cor-tabs","Functional")
        } else {
            hideTab("view-tabs","Resource info")
            hideTab("maintabs","Dev")
            hideTab("scell-tabs1","CNV")  ## DEV only
            hideTab("scell-tabs1","Monocle") ## DEV only       
        }
        
        ## Dynamically show upon availability in pgx object
        if(opt$ENABLE_UPLOAD) showTab("load-tabs","Upload data")            
        tabRequire(pgx, "connectivity", "maintabs", "Similar experiments")
        tabRequire(pgx, "drugs", "maintabs", "Drug connectivity")
        tabRequire(pgx, "wordcloud", "maintabs", "Word cloud")
        tabRequire(pgx, "deconv", "maintabs", "CellProfiling")
        fileRequire("tcga_matrix.h5", "maintabs", "TCGA survival (beta)")         
        if(!is.null(ACCESS.LOG)) showTab("load-tabs","Visitors map")                    

        message("[MAIN] reconfiguring menu done.")        
    })
    waiter_hide()
}

## --------------------------------------------------------------------
## ------------------------------ UI ----------------------------------
## --------------------------------------------------------------------

help.tabs <- navbarMenu(
    "Help",
    tabPanel(title=HTML("<a href='https://omicsplayground.readthedocs.io' target='_blank'>Documentation")),
    tabPanel(title=HTML("<a href='https://www.youtube.com/watch?v=_Q2LJmb2ihU&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-' target='_blank'>Video tutorials</a>")),
    tabPanel(title=HTML("<a href='https://github.com/bigomics/omicsplayground' target='_blank'>GitHub")),
    tabPanel(title=HTML("<a href='https://hub.docker.com/r/bigomics/omicsplayground' target='_blank'>Docker")),
    tabPanel(title=HTML("<a href='https://groups.google.com/d/forum/omicsplayground' target='_blank'>Google groups"))
)

TABVIEWS <- list(
    "load" = tabView("Home",LoadingInputs("load"),LoadingUI("load")),
    "view" = tabView("DataView",DataViewInputs("view"),DataViewUI("view")),
    "clust" = tabView("Unsupervised clustering",ClusteringInputs("clust"),ClusteringUI("clust")),
    "wgcna" = tabView("WGCNA (beta)",WgcnaInputs("wgcna"),WgcnaUI("wgcna")),
    "expr" = tabView("Differential expression",ExpressionInputs("expr"),ExpressionUI("expr")),
    "cor"  = tabView("Correlation analysis", CorrelationInputs("cor"), CorrelationUI("cor")),
    "enrich" = tabView("Geneset enrichment",EnrichmentInputs("enrich"), EnrichmentUI("enrich")),
    "func" = tabView("Pathway analysis", FunctionalInputs("func"), FunctionalUI("func")),
    "word" = tabView("Word cloud", WordCloudInputs("word"), WordCloudUI("word")),
    "drug" = tabView("Drug connectivity", DrugConnectivityInputs("drug"), DrugConnectivityUI("drug")),
    "isect" = tabView("Compare signatures", IntersectionInputs("isect"), IntersectionUI("isect")),
    "sig" = tabView("Test signatures", SignatureInputs("sig"), SignatureUI("sig")),
    "bio" = tabView("Find biomarkers", BiomarkerInputs("bio"), BiomarkerUI("bio")),
    "cmap" = tabView("Similar experiments", ConnectivityInputs("cmap"), ConnectivityUI("cmap")),
    "scell" = tabView("CellProfiling", SingleCellInputs("scell"), SingleCellUI("scell")),
    "tcga" = tabView("TCGA survival (beta)", TcgaInputs("tcga"), TcgaUI("tcga"))
    ##"system" = tabView("Systems analysis", SystemInputs("system"), SystemUI("system")),
    ##"multi" = tabView("Multi-level", MultiLevelInputs("multi"), MultiLevelUI("multi"))
)

names(TABVIEWS)
TABVIEWS <- TABVIEWS[names(TABVIEWS) %in% names(which(ENABLED))]
names(TABVIEWS)

TABVIEWS[["logout"]] <- tabPanel(title=HTML("<a id='logout' href='/logout'>Logout"))
ENABLED["logout"] <- ifelse(SHINYPROXY, TRUE, FALSE)
## ENABLED["logout"]=TRUE

tabs = list( home=c("load","logout"))
createUI <- function(tabs)
{
    message("\n======================================================")
    message("======================= UI ===========================")
    message("======================================================\n")
    
    version <- scan("../VERSION", character())[1]
    TITLE = paste(opt$TITLE,version)
    LOGO = div(img(src="bigomics-logo-white-48px.png", height="48px"),
               TITLE, id="navbar-logo", style="margin-top:-13px;")    
    title = tagList(LOGO)
    windowTitle = TITLE
    theme = shinythemes::shinytheme("cerulean")
    id = "maintabs"
    ##selected = "Home"    
    header = tagList(
        tags$head(tags$link(rel = "stylesheet", href = "playground.css")),
        tags$head(tags$link(rel="shortcut icon", href="favicon.ico")),
        shinyjs::useShinyjs(),
        firebase::useFirebase(),
        TAGS.JSSCRIPT,
        tags$script(async=NA, src="https://platform.twitter.com/widgets.js"),
        use_waiter(),
        div(textOutput("current_dataset"),class='current-data')
        ##QuestionBoard_UI("qa")
    )
    names(header) <- NULL
   
    footer = tagList(
        ##social_buttons(),
        waiter_show_on_load(
            html = spin_wave(),
            ##color = "#2780e3",
            color = "#1967be"            
            ## logo = "ready.png"
        ) # place at the bottom
    )

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
       
    ##-------------------------------------
    ## create navbarPage
    ##-------------------------------------
    selected = "Home"    
    names(tablist) <- NULL
    do.call( navbarPage, c(tablist,
                           title=title, id=id,
                           selected=selected,
                           windowTitle = windowTitle,
                           header = tagList(header),
                           footer = tagList(footer),
                           theme = theme))
}

tabs = list(
    "logout",
    "Home" = c("load"),
    "DataView" = "view",
    "Clustering" = c("clust","wgcna"),
    "Expression" = c("expr","cor"),
    "Enrichment" = c("enrich","func","word","drug"),
    "Signature" = c("isect","sig","bio","cmap","tcga"),
    "CellProfiling" = "scell",
    "Dev" = c("system","multi")
)
ui = createUI(tabs)


## --------------------------------------------------------------------
## ------------------------------ RUN ---------------------------------
## --------------------------------------------------------------------

shiny::shinyApp(ui, server)


##pkgs <- c( sessionInfo()[["basePkgs"]], names(sessionInfo()[["otherPkgs"]]),
##          names(sessionInfo()[["loadedOnly"]]) )

## --------------------------------------------------------------------
## ------------------------------ EOF----------------------------------
## --------------------------------------------------------------------
