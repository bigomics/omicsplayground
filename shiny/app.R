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
PGX.DIR="../data"
##PGX.DIR="~/bigomics/data/archs4data/gse1"

source("../R/pgx-init.R", local=TRUE)  ## pass local vars
DEV.VERSION = TRUE
options(shiny.maxRequestSize = 200*1024^2)  ##max 200Mb upload

if(0) {
    load("../data/geiger2016-arginine.pgx")
    load("../../alex/alex2019-data.pgx")
    load("../data/GSE10846-dlbcl.pgx")
    load("../data/GSE101766.pgx")    
    ngs=pgx.initialize(ngs)    
}


source("global.R")
source("../R/pgx-modules.R")
source("modules/DataViewModule.R", local=TRUE)
source("modules/ClusteringModule.R", local=TRUE)
source("modules/ExpressionModule.R", local=TRUE)
source("modules/EnrichmentModule.R", local=TRUE)
source("modules/IntersectionModule.R", local=TRUE)
source("modules/FunctionalModule.R", local=TRUE)
source("modules/SignatureModule.R", local=TRUE)
source("modules/LoadingModule.R", local=TRUE)
source("modules/BiomarkerModule.R", local=TRUE)
source("modules/ProfilingModule.R", local=TRUE)

server = function(input, output, session) {

    useShinyjs()
    ##useShinydashboard()
    ##useShinydashboardPlus()
    cat("===================== SERVER =======================\n")

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
    env[["bio"]]    <- callModule( BiomarkerModule, "bio", env)
    env[["prof"]]   <- callModule( ProfilingModule, "prof", env)

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
        hideTab("bio-tab1","Multi-level")

        if(usermode=="BASIC") {
            hideTab("maintabs","Biomarker")
            hideTab("maintabs","scProfiling")
            hideTab("clust-tabs2","Feature ranking")
            hideTab("expr-tabs1","Volcano (methods)")                        
            hideTab("expr-tabs2","FDR table")            
            hideTab("enrich-tabs1","Volcano (methods)")                        
            hideTab("enrich-tabs2","FDR table")            
        } else {
            showTab("maintabs","Biomarker")
            showTab("maintabs","scProfiling")
            showTab("clust-tabs2","Feature ranking")
            showTab("expr-tabs1","Volcano (methods)")                        
            showTab("expr-tabs2","FDR table")            
            showTab("enrich-tabs1","Volcano (methods)")                        
            showTab("enrich-tabs2","FDR table")            
        }
        if(DEV.VERSION) {
            showTab("view-tabs","Resource info")
            showTab("enrich-tabs1","GeneMap")                        
            showTab("bio-tab1","Multi-level")
        }
    })
    
    hide_waiter()    
}

tabView <- function(title, tab.inputs, tab.ui) {
    tabPanel(title,
             sidebarLayout(
                 sidebarPanel( width=2, tab.inputs, id="sidebar"),
                 mainPanel( width=10, tab.ui)
             ))
}

TITLE = "Omics PlayCloud"
TITLE = "Omics Playground v2"
logo = div(img(src="bigomics-logo-white-48px.png", height="48px"),
           TITLE, id="navbar-logo", style="margin-top:-13px;")

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
    tabView("Expression",ExpressionInputs("expr"),ExpressionUI("expr")),
    tabView("Enrichment",EnrichmentInputs("enrich"),EnrichmentUI("enrich")),
    tabView("Intersection", IntersectionInputs("isect"), IntersectionUI("isect")),
    tabView("Functional", FunctionalInputs("func"), FunctionalUI("func")),
    tabView("Signature", SignatureInputs("sig"), SignatureUI("sig")),
    tabView("Biomarker", BiomarkerInputs("bio"), BiomarkerUI("bio")),
    tabView("scProfiling", ProfilingInputs("prof"), ProfilingUI("prof")),

    footer = tagList(
        social_buttons(),
        show_waiter_on_load(spin_fading_circles()) # place at the bottom
    )
)

shiny::shinyApp(ui, server)
