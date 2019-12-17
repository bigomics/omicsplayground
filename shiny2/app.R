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
## DEV.VERSION = TRUE

source("../R/pgx-modules.R")
source("modules/DataViewModule.R", local=TRUE)
source("modules/ClusteringModule.R", local=TRUE)
source("modules/ExpressionModule.R", local=TRUE)
source("modules/EnrichmentModule.R", local=TRUE)
source("modules/IntersectionModule.R", local=TRUE)
source("modules/FunctionalModule.R", local=TRUE)
source("modules/SignatureModule.R", local=TRUE)
source("modules/LoadingModule.R", local=TRUE)

server = function(input, output, session) {

    ##useShinydashboard()
    ##useShinydashboardPlus()
    cat("===================== SERVER =======================\n")

    ## inputData <- reactive({ ngs })
    output$main_usermode <- renderText("PRO")
    
    inputData <- callModule(LoadingModule, "home1")
    callModule( DataViewModule, "dataview1", inputData)
    callModule( ClusteringModule, "clust1", inputData)
    callModule( ExpressionModule, "expr1", inputData)
    callModule( EnrichmentModule, "enrich1", inputData)
    callModule( IntersectionModule, "isect1", inputData)
    callModule( FunctionalModule, "func1", inputData)
    callModule( SignatureModule, "sig1", inputData)

    output$current_dataset <- renderText({        
        HTML("<div class='current-data'>",inputData()$name,"</div>")
    })
    
    hide_waiter()    
}

tabView <- function(title, tab.inputs, tab.ui) {
    tabPanel(title, sidebarLayout(
                        sidebarPanel( width = 2, tab.inputs ),
                        mainPanel( width = 10, tab.ui)
                    ))
}

title = div(img(src="bigomics-logo-white-48px.png", width="48px"),
            "Omics Playground v2", id="navbar-logo", style="margin-top:-13px;")

ui = navbarPage( 
    title = title, windowTitle="Omics Playground v2",
    theme = shinythemes::shinytheme("cerulean"),
    ##includeCSS("www/navbar.css"),
    header = tagList(
        tags$head( tags$link(rel = "stylesheet", href = "navbar.css")),
        shinyjs::useShinyjs(),        
        use_waiter(include_js = FALSE),
        htmlOutput("current_dataset")
    ),
    tabView("Home",LoadingInputs("home1"),LoadingUI("home1")),
    tabView("DataView",DataViewInputs("dataview1"),DataViewUI("dataview1")),
    tabView("Clustering",ClusteringInputs("clust1"),ClusteringUI("clust1")),
    tabView("Expression",ExpressionInputs("expr1"),ExpressionUI("expr1")),
    tabView("Enrichment",EnrichmentInputs("enrich1"),EnrichmentUI("enrich1")),
    tabView("Intersection", IntersectionInputs("isect1"), IntersectionUI("isect1")),
    tabView("Functional", FunctionalInputs("func1"), FunctionalUI("func1")),
    tabView("Signature", SignatureInputs("sig1"), SignatureUI("sig1")),
    footer = tagList(
        show_waiter_on_load(spin_fading_circles()) # place at the bottom
    )
)

shiny::shinyApp(ui, server)
