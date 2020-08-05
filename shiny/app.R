#########################################################################
##                                                                     ##
##              Main application for Omics Playground                  ##
##                                                                     ##
#########################################################################

message("\n\n")
message("#################################################################")
message("##################### OMICS PLAYGROUND ##########################")
message("#################################################################")
message("\n\n")

DEBUG = FALSE
DEBUG = TRUE
##.setSmtpServer("~/bigomics/server-conf/smtp_server/infomaniak.env")

library(shiny)
library(shinyjs)
library(shinyWidgets)
library(waiter)
library(plotly)

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
message("*******************************************")
message("******* SETTING GLOBAL VARIABLES **********")
message("*******************************************")

OPG = ".."
RDIR = file.path(OPG,"R")
FILES = file.path(OPG,"lib")
FILESX = file.path(OPG,"libx")
PGX.DIR = file.path(OPG,"data")
##PGX.DIR = c(PGX.DIR,file.path(OPG,"data-extra"))

source(file.path(RDIR,"pgx-include.R"),local=TRUE)  ## pass local vars
source("global.R", local=TRUE)
source("modules/AuthenticationDialog.R")

message("\n")
message("*****************************************")
message("******** parsing OPTIONS file ***********")
message("*****************************************")

source(file.path(RDIR,"pgx-files.R"), local=TRUE)  ## pass local vars
options(shiny.maxRequestSize = 999*1024^2)  ##max 999Mb upload
if(!file.exists("OPTIONS")) stop("FATAL ERROR: cannot find OPTIONS file")
opt <- pgx.readOptions(file="OPTIONS")

## over-ride options
## opt$AUTHENTICATION = "none"
## opt$AUTHENTICATION = "password"
## opt$AUTHENTICATION = "register"
if(Sys.getenv("PLAYGROUND_AUTHENTICATION")!="") {
    auth <- Sys.getenv("PLAYGROUND_AUTHENTICATION")
    message("[ENV] overriding PLAYGROUND_AUTHENTICATION = ",auth)
    opt$AUTHENTICATION = auth
}


WATERMARK      = opt$WATERMARK
SHOW_QUESTIONS = FALSE
DEV.VERSION    = opt$DEV_VERSION && dir.exists("../../omicsplayground-dev")
USER_MODE      = opt$USER_MODE
AUTHENTICATION = opt$AUTHENTICATION

## show options
message("\n",paste(paste(names(opt),"\t= ",sapply(opt,paste,collapse=" ")),collapse="\n"),"\n")

## --------------------------------------------------------------------
## ------------------------ READ FUNCTIONS ----------------------------
## --------------------------------------------------------------------

## pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)
pgx.initDatasetFolder(PGX.DIR, force=FALSE, verbose=1)
source("../R/pgx-init.R", local=TRUE)  ## pass local vars
##source("../R/pgx-functions.R", local=TRUE)  ## pass local vars

if(0) {
    load("../data/geiger2016-arginine.pgx")
    load("../data/GSE10846-dlbcl.pgx")
    load("../data/GSE10846-xgreta.pgx")
    load("../data/GSE72056-scmelanoma.pgx")
    load("../data/tcga-brca_pub.pgx")
    load("../data/GSE22886-immune.pgx")
    load("../data/pgx-grassi2020-10x-S.pgx")
    ngs = pgx.initialize(ngs)
}

## --------------------------------------------------------------------
## ------------------------ READ MODULES ------------------------------
## --------------------------------------------------------------------

MODULES <- c("load","view","clust","expr","enrich","isect","func",
             "word","drug","sig","scell","cor","bio","cmap",
             "tcga","bc","multi","qa")
if(is.null(opt$MODULES_ENABLED)) opt$MODULES_ENABLED = MODULES
if(is.null(opt$MODULES_DISNABLED)) opt$MODULES_DISNABLED = NA
ENABLED  <- array(MODULES %in% opt$MODULES_ENABLED, dimnames=list(MODULES))
DISABLED <- array(MODULES %in% opt$MODULES_DISABLED, dimnames=list(MODULES))
ENABLED  <- ENABLED & !DISABLED
ENABLED

modules <- dir("modules", pattern="Module.R$")
for(m in modules) {
    message("[MAIN] loading module ",m)
    source(paste0("modules/",m), local=TRUE)
    ##source(paste0("modules/",m), local=FALSE)
}

if(DEV.VERSION && dir.exists("../../omicsplayground-dev")) {
    xmodules <- dir("../../omicsplayground-dev/shiny/modules", pattern="Module..R$")
    for(m in xmodules) {
        message("[MAIN] loading module ",m)
        source(paste0("../../omicsplayground-dev/shiny/modules/",m), local=TRUE)
    }
    ENABLED[c("tcga","bc","multi")] <- TRUE
    modules <- c(modules, xmodules)
} else {
    ENABLED[c("tcga","bc","multi")] <- FALSE
}
ENABLED

## disable connectivity map if we have no signature database folder
has.sigdb <- length(dir(FILESX,pattern="sigdb.*h5")>0)
has.sigdb
if(has.sigdb==FALSE) ENABLED["cmap"] <- FALSE

MAINTABS = c("DataView","Clustering","Expression","Enrichment",
             "Signature","CellProfiling","Development")

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
    message("[MAIN] calling modules...")
    
    max.limits <- c("samples" = opt$MAX_SAMPLES,
                    "comparisons" = opt$MAX_COMPARISONS,
                    "genes" = opt$MAX_GENES)
    env <- list()  ## communication environment
    env[["load"]]   <- callModule(
        LoadingModule, "load", max.limits = max.limits,
        authentication = AUTHENTICATION)
    env[["view"]]   <- callModule( DataViewModule, "view", env)
    env[["clust"]]  <- callModule( ClusteringModule, "clust", env)
    env[["expr"]]   <- callModule( ExpressionModule, "expr", env)
    env[["enrich"]] <- callModule( EnrichmentModule, "enrich", env)
    if(ENABLED["func"])   env[["func"]]   <- callModule( FunctionalModule, "func", env)
    if(ENABLED["word"])   env[["word"]]   <- callModule( WordCloudModule, "word", env)
    if(ENABLED["drug"])   env[["drug"]]   <- callModule( DrugConnectivityModule, "drug", env)
    if(ENABLED["isect"])  env[["isect"]]  <- callModule( IntersectionModule, "isect", env)
    if(ENABLED["sig"])    env[["sig"]]    <- callModule( SignatureModule, "sig", env)
    if(ENABLED["scell"])  env[["scell"]]  <- callModule( SingleCellModule, "scell", env)
    if(ENABLED["cor"])    env[["cor"]]    <- callModule( CorrelationModule, "cor", env)
    if(ENABLED["bio"])    env[["bio"]]    <- callModule( BiomarkerModule, "bio", env)
    if(ENABLED["cmap"])   env[["cmap"]]   <- callModule( ConnectivityModule, "cmap", env)
    if(ENABLED["tcga"])   env[["tcga"]]   <- callModule( TcgaModule, "tcga", env)
    if(ENABLED["bc"])     env[["bc"]]     <- callModule( BatchCorrectModule, "bc", env)
    if(ENABLED["multi"])  env[["multi"]]  <- callModule( MultiLevelModule, "multi", env)
    env[["qa"]] <- callModule( QuestionModule, "qa", lapse = -1)
    
    ## message("[MAIN] all modules called:",paste(names(env),collapse=" "))
    message("[MAIN] modules enabled:",paste(names(which(ENABLED)),collapse=" "))
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
        
        hideTab("view-tabs","Resource info")
        hideTab("maintabs","Development")

        if(opt$USER_MODE == "BASIC") {
            hideTab("maintabs","CellProfiling")
            hideTab("enrich-tabs1","GeneMap")
            hideTab("clust-tabs2","Feature ranking")
            hideTab("expr-tabs1","Volcano (methods)")
            hideTab("expr-tabs2","FDR table")
            hideTab("enrich-tabs1","Volcano (methods)")
            hideTab("enrich-tabs2","FDR table")
        }

        if(!opt$USER_MODE == "DEV" || DEV.VERSION) {
            showTab("maintabs","Development")
            showTab("view-tabs","Resource info")
            showTab("enrich-tabs1","GeneMap")
            hideTab("scell-tabs1","CNV")  ## DEV only
            hideTab("scell-tabs1","Monocle") ## DEV only       
        }

        ## Dynamically show upon availability in pgx object
        if(opt$ENABLE_UPLOAD) showTab("load-tabs","Upload data")            
        showHideTab(pgx, "connectivity", "maintabs", "Similar experiments")
        showHideTab(pgx, "drugs", "maintabs", "Drug connectivity")
        showHideTab(pgx, "wordcloud", "maintabs", "Word cloud")
        showHideTab(pgx, "deconv", "maintabs", "CellProfiling")         
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
    "clust" = tabView("Clustering",ClusteringInputs("clust"),ClusteringUI("clust")),
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
    "scell" = tabView("CellProfiling", SingleCellInputs("scell"), SingleCellUI("scell"))    
)

if(DEV.VERSION) {
    TABVIEWS <- c(
        TABVIEWS,
        list(
            "bc" = tabView("Batch-effects analysis", BatchCorrectInputs("bc"), BatchCorrectUI("bc")),
            "tcga" = tabView("TCGA survival", TcgaInputs("tcga"), TcgaUI("tcga")),
            "multi" = tabView("Multi-level", MultiLevelInputs("multi"), MultiLevelUI("multi"))
        )
    )
}

tabs = "load"
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
        shinyjs::useShinyjs(),
        TAGS.JSSCRIPT,
        tags$script(async=NA, src="https://platform.twitter.com/widgets.js"),
        use_waiter(),
        div(textOutput("current_dataset"),class='current-data')
        ##QuestionModule_UI("qa")
    )
    names(header) <- NULL
    
    footer = tagList(
        ##social_buttons(),
        waiter_show_on_load(
            html = spin_wave(),
            ##color = "#2780e3",
            color = "#1967be",            
            logo = "ready.png"
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
    for(i in 1:length(tabs)) {
        itab <- tabs[[i]]
        itab <- itab[which(ENABLED[itab])] ## only enabled
        if(length(itab)>1) {
            m <- createNavbarMenu( names(tabs)[i], itab )
            tablist[[i]] <- m
        } else if(length(itab)==1) {
            tablist[[i]] <- TABVIEWS[[itab]] 
        } else {
        }
    }
    tablist <- tablist[!sapply(tablist,is.null)]
    
    ## add help menu
    tablist[["helpmenu"]] <- help.tabs
    names(tablist) <- NULL

    ##-------------------------------------
    ## create navbarPage
    ##-------------------------------------
    selected = "Home"    
    do.call( navbarPage, c(tablist,
                           title=title, id=id,
                           selected=selected,
                           windowTitle = windowTitle,
                           header = tagList(header),
                           footer = tagList(footer),
                           theme = theme))
}

ui = createUI(
    tabs = list(
        "Home" = "load",
        "DataView" = "view",
        "Clustering" = "clust",
        "Expression" = c("expr","cor"),
        "Enrichment" = c("enrich","func","word","drug"),
        "Signature" = c("isect","sig","bio","cmap"),
        "CellProfiling" = "scell",
        "Development" = c("bc","tcga","multi")
    )
)

## --------------------------------------------------------------------
## ------------------------------ RUN ---------------------------------
## --------------------------------------------------------------------

shiny::shinyApp(ui, server)

##pkgs <- c( sessionInfo()[["basePkgs"]], names(sessionInfo()[["otherPkgs"]]),
##          names(sessionInfo()[["loadedOnly"]]) )


## --------------------------------------------------------------------
## ------------------------------ EOF----------------------------------
## --------------------------------------------------------------------
