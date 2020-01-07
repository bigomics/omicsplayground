## install.packages("shinydashboardPlus")
## remotes::install_github("JohnCoene/waiter")

library(shiny)
library(shinyjs)
library(devtools)
require(shinyWidgets)
library(waiter)

cat("===================== INIT =======================\n")

RDIR="../R"
FILES="../lib"
PGX.DIR="~/bigomics/data/archs4data/gse"
PGX.DIR="~/Projects/Data/archs4data/gse"
PGX.DIR="../data"
dir.exists(PGX.DIR)

source("../R/pgx-functions.R", local=TRUE)  ## pass local vars
source("../R/pgx-files.R", local=TRUE)  ## pass local vars
##pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)
pgx.initDatasetFolder(PGX.DIR, force=FALSE, verbose=1)
source("../R/pgx-init.R", local=TRUE)  ## pass local vars
options(shiny.maxRequestSize = 200*1024^2)  ##max 200Mb upload

## DEV.VERSION = TRUE

if(0) {
    load("../data/geiger2016-arginine.pgx")
    load("../../alex/alex2019-data.pgx")
    load("../data/GSE10846-dlbcl.pgx")
    load("../data/GSE101766.pgx")
    load("../data-ext/vogel2019-tcells.pgx")
    load("../../exampledata/mouse/mouse.pgx")    
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
source("modules/SignatureModule.R", local=TRUE)
source("modules/ProfilingModule.R", local=TRUE)
source("modules/CorrelationModule.R", local=TRUE)
if(DEV.VERSION) source("../../omicsplayground-dev/shiny/modules/BiomarkerModule.R", local=TRUE)
if(DEV.VERSION) source("../../omicsplayground-dev/shiny/modules/MetaModule.R", local=TRUE)
if(DEV.VERSION) source("../../omicsplayground-dev/shiny/modules/BatchCorrectModule.R", local=TRUE)

server = function(input, output, session) {

    useShinyjs()
    ##useShinydashboard()
    ##useShinydashboardPlus()
    cat("===================== SERVER =======================\n")
    cat("calling modules... ")
    
    env <- list()  ## communication environment 
    ## env[["load"]][["inputData"]] <- reactive({ ngs })    
    env[["load"]]   <- callModule( LoadingModule, "load", hideUserMode=FALSE)
    env[["view"]]   <- callModule( DataViewModule, "view", env)
    env[["clust"]]  <- callModule( ClusteringModule, "clust", env)
    env[["expr"]]   <- callModule( ExpressionModule, "expr", env)
    env[["enrich"]] <- callModule( EnrichmentModule, "enrich", env)
    env[["isect"]]  <- callModule( IntersectionModule, "isect", env)
    env[["func"]]   <- callModule( FunctionalModule, "func", env)
    env[["sig"]]    <- callModule( SignatureModule, "sig", env)
    env[["prof"]]   <- callModule( ProfilingModule, "prof", env)
    env[["cor"]]    <- callModule( CorrelationModule, "cor", env)
    if(DEV.VERSION) env[["bio"]]  <- callModule( BiomarkerModule, "bio", env)
    if(DEV.VERSION) env[["meta"]] <- callModule( MetaModule, "meta", env)
    if(DEV.VERSION) env[["bc"]]   <- callModule( BatchCorrectModule, "bc", env)

    cat("[OK]\n")

    output$current_dataset <- renderText({
        pgx <- env[["load"]][["inputData"]]()
        name <- gsub(".*\\/|[.]pgx$","",pgx$name)
        if(length(name)==0) name = "(no data)"
        name
    })

    ## Hide/show certain section
    observe({
        usermode <- env[["load"]][["usermode"]]()
        if(length(usermode)==0) usermode <- "BASIC"
        dbg("usermode = ",usermode)
        hideTab("view-tabs","Resource info")
        hideTab("enrich-tabs1","GeneMap")

        hideTab("maintabs","Development")
        ## hideTab("maintabs","BatchCorrect")
        ## hideTab("maintabs","Biomarker analysis")
        ## hideTab("maintabs","MetaAnalysis")
        ## hideTab("bio-tab1","Multi-level")

        if(usermode=="BASIC") {
            hideTab("maintabs","scProfiling")
            hideTab("clust-tabs2","Feature ranking")
            hideTab("expr-tabs1","Volcano (methods)")                        
            hideTab("expr-tabs2","FDR table")            
            hideTab("enrich-tabs1","Volcano (methods)")                        
            hideTab("enrich-tabs2","FDR table")            
        } else {
            showTab("maintabs","scProfiling")
            showTab("clust-tabs2","Feature ranking")
            showTab("expr-tabs1","Volcano (methods)")                        
            showTab("expr-tabs2","FDR table")            
            showTab("enrich-tabs1","Volcano (methods)")                        
            showTab("enrich-tabs2","FDR table")            
        }
        if(DEV.VERSION) {
            showTab("maintabs","Development")            
            ## showTab("maintabs","Biomarker analysis")
            ## showTab("maintabs","BatchCorrect")
            ## showTab("maintabs","MetaAnalysis")
            showTab("view-tabs","Resource info")
            showTab("enrich-tabs1","GeneMap")                        
            showTab("bio-tab1","Multi-level")
        }
    })
    
    hide_waiter()    
}


version <- scan("../VERSION", character())[1]
TITLE = paste("Omics Playground",version)
## TITLE = "Omics PlayCloud"
logo = div(img(src="bigomics-logo-white-48px.png", height="48px"),
           TITLE, id="navbar-logo", style="margin-top:-13px;")

dev.tabs <- NULL
if(DEV.VERSION) {
    dev.tabs <- navbarMenu(
        "Development",
        tabView("Biomarker analysis", BiomarkerInputs("bio"), BiomarkerUI("bio")),        
        tabView("BatchCorrect", BatchCorrectInputs("bc"), BatchCorrectUI("bc")),
        tabView("MetaAnalysis", MetaInputs("meta"), MetaUI("meta"))
    )
}


ui = navbarPage( 
    title = logo, windowTitle = TITLE,
    theme = shinythemes::shinytheme("cerulean"),
    ##includeCSS("www/navbar.css"),
    id = "maintabs",
    header = tagList(
        tags$head(tags$link(rel = "stylesheet", href = "navbar.css")),
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
        tabView("Functional analysis", FunctionalInputs("func"), FunctionalUI("func"))
    ),
    navbarMenu(
        "Signature",
        tabView("Intersection analysis", IntersectionInputs("isect"), IntersectionUI("isect")),
        tabView("Signature analysis", SignatureInputs("sig"), SignatureUI("sig"))
    ),
    tabView("scProfiling", ProfilingInputs("prof"), ProfilingUI("prof")),
    dev.tabs,
    footer = tagList(
        social_buttons(),
        show_waiter_on_load(spin_fading_circles()) # place at the bottom
    )
)

shiny::shinyApp(ui, server)

