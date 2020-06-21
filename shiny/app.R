#########################################################################
##                                                                     ##
##              Main application for Omics Playground                  ##
##                                                                     ##
#########################################################################


## --------------------------------------------------------------------
## ------------------------ CHECKS ------------------------------------
## --------------------------------------------------------------------

res <- system("ping orca-server",intern=TRUE)
orca.ok <- (length(res)>0)
orca.ok
if(!orca.ok) {
    cat("###############################################################\n")
    cat("##### ERROR:: ORCA server not running. please start ORCA. #####\n")
    cat("###############################################################\n")
    stop()
}


## --------------------------------------------------------------------
## ------------------------- INIT -------------------------------------
## --------------------------------------------------------------------

cat("===================== INIT =======================\n")

library(shiny)
library(shinyjs)
library(shinyWidgets)
library(waiter)

RDIR = "../R"
FILES = "../lib"
FILESX = "../libx"
PGX.DIR = c("../data","../data-extra")
PGX.DIR = "../data"
dir.exists(PGX.DIR)

## --------------------------------------------------------------------
## ----------------------- READ OPTIONS -------------------------------
## --------------------------------------------------------------------

source("../R/pgx-files.R", local=TRUE)  ## pass local vars
options(shiny.maxRequestSize = 999*1024^2)  ##max 999Mb upload
if(!file.exists("OPTIONS")) stop("FATAL ERROR: cannot find OPTIONS file")
opt <- pgx.readOptions(file="OPTIONS")

WATERMARK = opt$WATERMARK
SHOW_QUESTIONS = FALSE
DEV.VERSION = opt$DEV_VERSION && dir.exists("../../omicsplayground-dev")
if(opt$USER_MODE=="BASIC") DEV.VERSION = FALSE

## show options
cat(paste(paste(names(opt), "\t= ", sapply(opt,paste,collapse=" ")),collapse="\n"),"\n")

## --------------------------------------------------------------------
## ------------------------ READ FUNCTIONS ----------------------------
## --------------------------------------------------------------------

source("../R/pgx-include.R", local=TRUE)  ## pass local vars
## pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)
pgx.initDatasetFolder(PGX.DIR, force=FALSE, verbose=1)
source("../R/pgx-init.R", local=TRUE)  ## pass local vars
source("global.R", local=TRUE)

if(0) {
    load("../data/geiger2016-arginine.pgx")
    load("../data/GSE10846-dlbcl.pgx")
    load("../data/GSE102908-ibetX.pgx")
    load("../data/tcga-brca_pub.pgx")
    load("../data/GSE22886-immune.pgx")   
    ngs = pgx.initialize(ngs)
}

## --------------------------------------------------------------------
## ------------------------ READ MODULES ------------------------------
## --------------------------------------------------------------------

MODULES <- c("load","view","clust","expr","enrich","isect","func",
             "word","drug","sig","scell","cor","bio","cmap",
             "tcga","bc","multi","qa")
if(is.null(opt$MODULES_ENABLED)) opt$MODULES_ENABLED = MODULES
ENABLED  <- array(MODULES %in% opt$MODULES_ENABLED, dimnames=list(MODULES))
DISABLED <- array(MODULES %in% opt$MODULES_DISABLED, dimnames=list(MODULES))
ENABLED  <- ENABLED & !DISABLED
ENABLED

cat("[MAIN] sourcing modules...\n")

source("modules/LoadingModule.R", local=TRUE)
source("modules/DataViewModule.R", local=TRUE)
source("modules/ClusteringModule.R", local=TRUE)
source("modules/ExpressionModule.R", local=TRUE)
source("modules/EnrichmentModule.R", local=TRUE)
source("modules/IntersectionModule.R", local=TRUE)
source("modules/FunctionalModule.R", local=TRUE)
source("modules/WordCloudModule.R", local=TRUE)
source("modules/DrugConnectivityModule.R", local=TRUE)
source("modules/SignatureModule.R", local=TRUE)
source("modules/SingleCellModule.R", local=TRUE)
source("modules/CorrelationModule.R", local=TRUE)
source("modules/BiomarkerModule.R", local=TRUE)
source("modules/QuestionModule.R", local=TRUE)
source("modules/ConnectivityModule.R", local=TRUE)
##source("modules/UsersMapModule.R", local=TRUE)

if(DEV.VERSION && dir.exists("../../omicsplayground-dev")) {
    source("../../omicsplayground-dev/shiny/modules/TcgaModule.R", local=TRUE)
    source("../../omicsplayground-dev/shiny/modules/BatchCorrectModule.R", local=TRUE)
    source("../../omicsplayground-dev/shiny/modules/MultiLevelModule.R", local=TRUE)
    ENABLED[c("tcga","bc","multi")] <- TRUE
} else {
    ENABLED[c("tcga","bc","multi")] <- FALSE
}
ENABLED

cat("[MAIN] starting server and UI...\n")

MAINTABS = c("DataView","Clustering","Expression","Enrichment",
             "Signature","SingleCell","Development")

server = function(input, output, session) {

    ## useShinyjs()
    ## useShinydashboard()
    ## useShinydashboardPlus()

    cat("===================== SERVER =======================\n")
    cat("[MAIN] calling modules...\n")

    max.limits <- c("samples" = opt$MAX_SAMPLES,
                    "comparisons" = opt$MAX_COMPARISONS,
                    "genes" = opt$MAX_GENES)
    env <- list()  ## communication environment
    env[["load"]]   <- callModule(
        LoadingModule, "load", hideModeButton = opt$HIDE_MODEBUTTON,
        max.limits = max.limits, defaultMode = opt$USER_MODE )
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
    env[["qa"]]     <- callModule( QuestionModule, "qa", lapse = -1)
    
    cat("[MAIN] all modules called\n")
    
    output$current_dataset <- renderText({
        pgx <- env[["load"]][["inputData"]]()
        name <- gsub(".*\\/|[.]pgx$","",pgx$name)
        if(length(name)==0) name = "(no data)"
        name
    })
   
    ## Hide/show certain sections depending on USER MODE
    observe({
        pgx <- env[["load"]][["inputData"]]() ## trigger on change dataset

        ## hide all main tabs
        if(is.null(pgx)) {
            lapply(MAINTABS, function(m) hideTab("maintabs",m))
            if(!opt$ENABLE_UPLOAD)  hideTab("load-tabs","Upload data")
            if(is.null(ACCESS.LOG)) hideTab("load-tabs","Visitors map")            
            return(NULL)
        }

        ## show all main tabs
        lapply(MAINTABS, function(m) showTab("maintabs",m))
                
        ## show single-cell module??
        show.cc <- ( (opt$SINGLE_CELL == "AUTO" && ncol(pgx$counts) >= 500) ||
                     (opt$SINGLE_CELL == "AUTO" && grepl("^scRNA",pgx$datatype)) ||
                     opt$SINGLE_CELL == "TRUE")
        if(is.null(show.cc) || is.na(show.cc) || length(show.cc)==0) show.cc <- FALSE
        show.cc <- show.cc && "deconv" %in% names(pgx)
        
        hideTab("view-tabs","Resource info")
        hideTab("maintabs","Development")
        hideTab("maintabs","SingleCell")
        
        hideTab("enrich-tabs1","GeneMap")
        hideTab("clust-tabs2","Feature ranking")
        hideTab("expr-tabs1","Volcano (methods)")
        hideTab("expr-tabs2","FDR table")
        hideTab("enrich-tabs1","Volcano (methods)")
        hideTab("enrich-tabs2","FDR table")

        hideTab("maintabs","SingleCell")
        hideTab("scell-tabs1","CNV")
        hideTab("scell-tabs1","Monocle")        
        if(show.cc) showTab("maintabs","SingleCell")

        if(opt$USER_MODE == "PRO") {
            showTab("clust-tabs2","Feature ranking")
            showTab("expr-tabs1","Volcano (methods)")
            showTab("expr-tabs2","FDR table")
            showTab("enrich-tabs1","Volcano (methods)")
            showTab("enrich-tabs2","FDR table")
        }

        if(DEV.VERSION) {
            showTab("maintabs","Development")
            showTab("view-tabs","Resource info")
            showTab("enrich-tabs1","GeneMap")
            showTab("scell-tabs1","CNV")
            showTab("scell-tabs1","Monocle")
        }

        ## Dynamically show upon availability
        if(opt$ENABLE_UPLOAD) showTab("load-tabs","Upload data")            
        showHideTab(pgx, "connectivity", "maintabs", "Similar experiments")
        showHideTab(pgx, "drugs", "maintabs", "Drug connectivity")
        showHideTab(pgx, "wordcloud", "maintabs", "Word cloud") 
        if(!is.null(ACCESS.LOG)) showTab("load-tabs","Visitors map")            
        
    })

    waiter_hide()
}

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
    "scell" = tabView("SingleCell", SingleCellInputs("scell"), SingleCellUI("scell"))    
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
createNavbarPage <- function(tabs)
{
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
        tags$head(tags$link(rel = "stylesheet", href = "navbar.css")),
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
        waiter_show_on_load(spin_fading_circles()) # place at the bottom
    )

    ## create TAB list
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

    selected = "Home"    
    do.call( navbarPage, c(tablist,
                           title=title, id=id,
                           selected=selected,
                           windowTitle = windowTitle,
                           header = tagList(header),
                           footer = tagList(footer),
                           theme = theme) )

}

ui = createNavbarPage(
    tabs = list(
        "Home" = "load",
        "DataView" = "view",
        "Clustering" = "clust",
        "Expression" = c("expr","cor"),
        "Enrichment" = c("enrich","func","word","drug"),
        "Signature" = c("isect","sig","bio","cmap"),
        "SingleCell" = "scell",
        "Development" = c("bc","tcga","multi")
    )
)

shiny::shinyApp(ui, server)


##pkgs <- c( sessionInfo()[["basePkgs"]], names(sessionInfo()[["otherPkgs"]]),
##          names(sessionInfo()[["loadedOnly"]]) )
