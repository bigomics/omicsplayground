##################################################################################
##                                                                              ##
##                 Main application for Omics Playground                        ##  
##                                                                              ##
##################################################################################

library(shiny)
library(shinyjs)
library(devtools)
require(shinyWidgets)
library(waiter)

cat("===================== INIT =======================\n")

RDIR = "../R"
FILES = "../lib"
##PGX.DIR = "../data"
##PGX.DIR = "/data/PublicData/archs4data/gse25k"
PGX.DIR = c("../data","../data-extra")
dir.exists(PGX.DIR)

source("../R/pgx-include.R", local=TRUE)  ## pass local vars
## pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)
pgx.initDatasetFolder(PGX.DIR, force=FALSE, verbose=1)
source("../R/pgx-init.R", local=TRUE)  ## pass local vars

options(shiny.maxRequestSize = 999*1024^2)  ##max 999Mb upload
if(!file.exists("OPTIONS")) stop("FATAL ERROR: cannot find OPTIONS file")
opt <- pgx.readOptions(file="OPTIONS")

DEV.VERSION = FALSE
if(dir.exists("../../omicsplayground-dev")) DEV.VERSION = TRUE

if(opt$USER_MODE=="basic") {
    cat("********************* BASIC MODE **********************\n")
    DEV.VERSION = FALSE
}

if(0) {
    load("../data/geiger2016-arginine.pgx")
    load("../data/GSE10846-dlbcl.pgx")
    load("../data/am2019-jev.pgx")
    load("../data/tcga-brca_pub-gx.pgx")
    load("../../omicsplayground-dev/data/CCLE-drugSX2.pgx")
    ngs = pgx.initialize(ngs)
}

source("global.R", local=TRUE)
source("../R/pgx-modules.R", local=TRUE)
source("modules/LoadingModule.R", local=TRUE)
source("modules/DataViewModule.R", local=TRUE)
source("modules/ClusteringModule.R", local=TRUE)
source("modules/ExpressionModule.R", local=TRUE)
source("modules/EnrichmentModule.R", local=TRUE)
source("modules/IntersectionModule.R", local=TRUE)
source("modules/FunctionalModule.R", local=TRUE)
source("modules/DrugConnectivityModule.R", local=TRUE)
source("modules/SignatureModule.R", local=TRUE)
source("modules/ProfilingModule.R", local=TRUE)
source("modules/CorrelationModule.R", local=TRUE)
source("modules/BiomarkerModule.R", local=TRUE)

if(DEV.VERSION && dir.exists("../../omicsplayground-dev")) {
    source("../../omicsplayground-dev/shiny/modules/ConnectivityModule.R", local=TRUE)
    source("../../omicsplayground-dev/shiny/modules/TcgaModule.R", local=TRUE)
    source("../../omicsplayground-dev/shiny/modules/BatchCorrectModule.R", local=TRUE)
    source("../../omicsplayground-dev/shiny/modules/MultiLevelModule.R", local=TRUE)
}

server = function(input, output, session) {

    useShinyjs()
    ##useShinydashboard()
    ##useShinydashboardPlus()
    cat("===================== SERVER =======================\n")
    cat("calling modules... ")

    max.limits <- c("samples" = opt$MAX_SAMPLES,
                    "comparisons" = opt$MAX_COMPARISONS,
                    "genes" = opt$MAX_GENES)
    
    env <- list()  ## communication environment
    ## env[["load"]][["inputData"]] <- reactive({ ngs })
    env[["load"]]   <- callModule(
        LoadingModule, "load", hideModeButton=opt$HIDE_MODEBUTTON,
        max.limits = max.limits, defaultMode=opt$USER_MODE )
    env[["view"]]   <- callModule( DataViewModule, "view", env)
    env[["clust"]]  <- callModule( ClusteringModule, "clust", env)
    env[["expr"]]   <- callModule( ExpressionModule, "expr", env)
    env[["enrich"]] <- callModule( EnrichmentModule, "enrich", env)
    env[["isect"]]  <- callModule( IntersectionModule, "isect", env)
    env[["func"]]   <- callModule( FunctionalModule, "func", env)
    env[["drug"]]   <- callModule( DrugConnectivityModule, "drug", env)
    env[["sig"]]    <- callModule( SignatureModule, "sig", env)
    env[["prof"]]   <- callModule( ProfilingModule, "prof", env)
    env[["cor"]]    <- callModule( CorrelationModule, "cor", env)
    env[["bio"]]    <- callModule( BiomarkerModule, "bio", env)
    if(DEV.VERSION) {
        env[["cmap"]] <- callModule( ConnectivityModule, "cmap", env)
        env[["tcga"]] <- callModule( TcgaModule, "tcga", env)
        env[["bc"]]   <- callModule( BatchCorrectModule, "bc", env)
        env[["multi"]]   <- callModule( MultiLevelModule, "multi", env)
    }

    cat("[OK]\n")

    output$current_dataset <- renderText({
        pgx <- env[["load"]][["inputData"]]()
        name <- gsub(".*\\/|[.]pgx$","",pgx$name)
        if(length(name)==0) name = "(no data)"
        name
    })

    ## Timed UI messages...
    nwarn = 0
    observe({
        usermode <- env[["load"]][["usermode"]]()
        if(opt$USER_MODE=="basic") usermode <- "BASIC" ## override
        if(usermode=="BASIC") {
            shinyjs::hide(selector = "div.download-button")
            shinyjs::hide(selector = "div.modebar")
            shinyjs::hide(selector = "div.pro-feature")
            ## if(nwarn==3) sendSweetAlert( session=session, title="", text="download is disabled")
        }
        invalidateLater(1000*30)  ## every 30 seconds check...
        nwarn <<- nwarn + 1
    })

    ## Hide/show certain sections depending on USER MODE
    observe({
        pgx <- env[["load"]][["inputData"]]() ## trigger on change dataset
        usermode <- env[["load"]][["usermode"]]()  ## trigger on button
        if(length(usermode)==0) usermode <- "BASIC"
        dbg("usermode = ",usermode)
        
        hideTab("view-tabs","Resource info")
        hideTab("maintabs","Development")
        hideTab("maintabs","Biomarker analysis")
        hideTab("maintabs","Drug connectivity")
        hideTab("maintabs","SingleCell")
        shinyjs::hide(selector = "div.download-button")
        shinyjs::hide(selector = "div.modebar")
        shinyjs::hide(selector = "div.pro-feature")

        hideTab("enrich-tabs1","GeneMap")
        hideTab("clust-tabs2","Feature ranking")
        hideTab("expr-tabs1","Volcano (methods)")
        hideTab("expr-tabs2","FDR table")
        hideTab("enrich-tabs1","Volcano (methods)")
        hideTab("enrich-tabs2","FDR table")
        hideTab("prof-tabs1","Monocle")

        if(toupper(opt$ENABLE_UPLOAD) %in% c("NO","FALSE")) {
            hideTab("load-tabs","Upload data")            
        }
        
        if(usermode != "BASIC") {
            showTab("maintabs","SingleCell")
            showTab("maintabs","Biomarker analysis")
            showTab("maintabs","Drug connectivity")

            showTab("clust-tabs2","Feature ranking")
            showTab("expr-tabs1","Volcano (methods)")
            showTab("expr-tabs2","FDR table")
            showTab("enrich-tabs1","Volcano (methods)")
            showTab("enrich-tabs2","FDR table")
            shinyjs::show(selector = "div.download-button")
            shinyjs::show(selector = "div.modebar")
            shinyjs::show(selector = "div.pro-feature")
        }

        if(DEV.VERSION) {
            showTab("maintabs","Development")
            showTab("view-tabs","Resource info")
            showTab("enrich-tabs1","GeneMap")
            showTab("prof-tabs1","Monocle")
        }
        
    })

    waiter_hide()
}

version <- scan("../VERSION", character())[1]
TITLE = paste(opt$TITLE,version)
## TITLE = "Omics PlayCloud"
logo = div(img(src="bigomics-logo-white-48px.png", height="48px"),
           TITLE, id="navbar-logo", style="margin-top:-13px;")

dev.tabs <- NULL
if(DEV.VERSION) {
    dev.tabs <- navbarMenu(
        "Development",
        tabView("Batch-effects analysis", BatchCorrectInputs("bc"), BatchCorrectUI("bc")),
        tabView("Connectivity mapping", ConnectivityInputs("cmap"), ConnectivityUI("cmap")),
        tabView("TCGA survival", TcgaInputs("tcga"), TcgaUI("tcga")),
        tabView("Multi-level", MultiLevelInputs("multi"), MultiLevelUI("multi"))
    )
}

help.tabs <- navbarMenu(
    "Help",
    tabPanel(title=HTML("<a href='https://omicsplayground.readthedocs.io' target='_blank'>Documentation")),
    tabPanel(title=HTML("<a href='https://www.youtube.com/watch?v=_Q2LJmb2ihU&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-' target='_blank'>Video tutorials</a>")),
    tabPanel(title=HTML("<a href='https://github.com/bigomics/omicsplayground' target='_blank'>GitHub")),
    tabPanel(title=HTML("<a href='https://hub.docker.com/r/bigomics/omicsplayground' target='_blank'>Docker")),
    tabPanel(title=HTML("<a href='https://groups.google.com/d/forum/omicsplayground' target='_blank'>Google groups"))
)


ui = navbarPage(
    title = logo, windowTitle = TITLE,
    theme = shinythemes::shinytheme("cerulean"),
    ##includeCSS("www/navbar.css"),
    id = "maintabs",
    selected = "Home",
    header = tagList(
        tags$head(tags$link(rel = "stylesheet", href = "navbar.css")),
        TAGS.JSSCRIPT,
        shinyjs::useShinyjs(),
        use_waiter(),
        div(textOutput("current_dataset"),class='current-data')
    ),
    tabView("Home",LoadingInputs("load"),LoadingUI("load")),
    tabView("DataView",DataViewInputs("view"),DataViewUI("view")),
    tabView("Clustering",ClusteringInputs("clust"),ClusteringUI("clust")),
    navbarMenu(
        "Expression",
        tabView("Differential expression",ExpressionInputs("expr"),ExpressionUI("expr")),
        tabView("Correlation analysis", CorrelationInputs("cor"), CorrelationUI("cor"))
    ),
    navbarMenu(
        "Enrichment",
        tabView("Geneset enrichment",EnrichmentInputs("enrich"),EnrichmentUI("enrich")),
        tabView("Pathway analysis", FunctionalInputs("func"), FunctionalUI("func")),
        tabView("Drug connectivity", DrugConnectivityInputs("drug"), DrugConnectivityUI("drug"))
    ),
    navbarMenu(
        "Signature",
        tabView("Intersection analysis", IntersectionInputs("isect"), IntersectionUI("isect")),
        tabView("Signature analysis", SignatureInputs("sig"), SignatureUI("sig")),
        tabView("Biomarker analysis", BiomarkerInputs("bio"), BiomarkerUI("bio"))
    ),
    tabView("SingleCell", ProfilingInputs("prof"), ProfilingUI("prof")),
    help.tabs,
    dev.tabs,
    footer = tagList(
        ## social_buttons(),
        waiter_show_on_load(spin_fading_circles()) # place at the bottom
    )
)

shiny::shinyApp(ui, server)


