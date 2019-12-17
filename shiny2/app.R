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
source("../R/pgx-init.R", local=TRUE)  ## pass local vars

##load("../data/geiger2016-arginine.pgx"); ngs=pgx.initialize(ngs)    

USERMODE = reactiveVal("PRO")
DEV.VERSION = TRUE

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

    ##useShinydashboard()
    ##useShinydashboardPlus()
    cat("===================== SERVER =======================\n")

    env <- list()  ## communication environment 
    ## env[["load"]][["inputData"]] <- reactive({ ngs })    
    env[["load"]]   <- callModule( LoadingModule, "load")
    env[["expr"]]   <- callModule( ExpressionModule, "expr", env)
    env[["view"]]   <- callModule( DataViewModule, "view", env)
    env[["clust"]]  <- callModule( ClusteringModule, "clust", env)
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
        HTML("<div class='current-data'>",name,"</div>")
    })

    ## Hide/show certain section
    observe({
        usermode <- env[["load"]][["usermode"]]()
        if(length(usermode)==0) return(NULL)
        dbg("usermode = ",usermode)
        if(usermode=="BASIC") {
            hideTab("maintabs","Biomarker")
            hideTab("maintabs","scProfiling")
        } else {
            showTab("maintabs","Biomarker")
            showTab("maintabs","scProfiling")
        }
    })
    
    waiter_hide()    
}

tabView <- function(title, tab.inputs, tab.ui) {
    tabPanel(title,
             sidebarLayout(
                 sidebarPanel( width=2, tab.inputs, id="sidebar"),
                 mainPanel( width=10, tab.ui)
             ))
}

title = div(img(src="bigomics-logo-white-48px.png", width="48px"),
            "Omics Playground v2", id="navbar-logo", style="margin-top:-13px;")

ui = navbarPage( 
    title = title, windowTitle="Omics Playground v2",
    theme = shinythemes::shinytheme("cerulean"),
    ##includeCSS("www/navbar.css"),
    id = "maintabs",
    header = tagList(
        tags$head( tags$link(rel = "stylesheet", href = "navbar.css")),
        shinyjs::useShinyjs(),        
        use_waiter(include_js = FALSE),
        htmlOutput("current_dataset")
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
        waiter_show_on_load(spin_fading_circles()) # place at the bottom
    )
)

shiny::shinyApp(ui, server)
