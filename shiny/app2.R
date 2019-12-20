## install.packages("shinydashboardPlus")
## remotes::install_github("JohnCoene/waiter")

library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinydashboardPlus)
library(devtools)
require(shinyWidgets)
library(waiter)
 
cat("===================== INIT =======================\n")

RDIR="../R"
FILES="../lib"
PGX.DIR="../data"
source("../R/pgx-init.R", local=TRUE)  ## pass local vars
DEV.VERSION = TRUE

##load("../data/geiger2016-arginine.pgx"); ngs=pgx.initialize(ngs)    

source("../R/pgx-modules.R")
source("modules/LoadingModule.R", local=TRUE)
source("modules/DataViewModule.R", local=TRUE)
source("modules/ClusteringModule.R", local=TRUE)
source("modules/ExpressionModule.R", local=TRUE)
source("modules/EnrichmentModule.R", local=TRUE)
source("modules/IntersectionModule.R", local=TRUE)
source("modules/FunctionalModule.R", local=TRUE)
source("modules/SignatureModule.R", local=TRUE)
source("modules/BiomarkerModule.R", local=TRUE)
source("modules/ProfilingModule.R", local=TRUE)


server = function(input, output, session) {

    ##useShinydashboard()
    useShinydashboardPlus()
    js$showControlSidebar()
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
    ##env[["bio"]]    <- callModule( BiomarkerModule, "bio", env)
    env[["prof"]]   <- callModule( ProfilingModule, "prof", env)

    ## How to hide/show the right (control)sidebar
    observeEvent(input$menuitem, {
        ##js$hideControlSidebar()
        state <- isolate(input$sidebarCollapsed)
        print(paste("sidebarCollapsed=",state))
        print(paste("sidebarItemExpanded=",isolate(input$sidebarItemExpanded)))
        cat("names(input)=",names(input),"\n")
        if(state) js$hideControlSidebar()
        if(!state) js$showControlSidebar()
    })    
    
    output$sidebarmenu <- renderMenu({
        ## Build sidebar menu
        menuitems <-  tagList(
            menuItem("Home", tabName = "load", icon=icon("home")),
            menuItem("DataView", tabName = "view", icon=icon("table")),            
            menuItem("Clustering", tabName = "clust", icon=icon("project-diagram")),
            menuItem("Expression", tabName = "expr", icon=icon("chart-bar")),
            menuItem("Enrichment", tabName = "enrich", icon=icon("expand-arrows-alt")),
            menuItem("Intersection", tabName = "isect", icon=icon("crosshairs")),
            menuItem("Functional", tabName = "func", icon=icon("square-root-alt")),
            menuItem("Signature", tabName = "sig", icon=icon("signature")),
            ##menuItem("Biomarker", tabName = "bio", icon=icon("marker")),
            menuItem("scProfiling", tabName = "prof", icon=icon("dot-circle"))
        )
        sidebarMenu( id="menuitem", .list = menuitems)
    })

    output$current_dataset <- renderText({
        pgx <- env[["load"]][["inputData"]]()
        name <- gsub(".*\\/|[.]pgx$","",pgx$name)
        if(length(name)==0) name = "(no data)"
        HTML("<div class='current-data'>",name,"</div>")
    })
    
    hide_waiter()    
}

footer.txt = "Powered by <a href='http://bigomics.ch'>BigOmics Analytics</a>."

logo = tagList(span(class="logo-lg", "Omics Playground v2"), 
               tags$img(src="bigomics-logo-white-32px.png"))

logo = tagList(div(img(src="company-logo.svg", height="36px"), id="navbar-logo", class="logo-lg", style="margin-top:-4px;"), div(img(src="company-logo.svg", height="36px"), id="navbar-logo", style="margin-top:-4px; margin-left:-4px;"))

ui = dashboardPagePlus( 
    title = "Omics Playground v2",
    ##title = "BC | Platforms",    
    skin = "blue",
    header = dashboardHeaderPlus(
        title = logo,
        enable_rightsidebar = TRUE,
        rightSidebarIcon = "ellipsis-v",
        ## ---------- items in the top menu aligned left:
        left_menu = tagList(
            dropdownButton(
                inputId="file-menu", ## 
                label="File", icon=NULL, circle=FALSE,
                actionLink(inputId="file_upload", label="Upload"),
                actionLink(inputId="file_load", label="Load"),
                actionLink(inputId="file_save", label="Save")
            ),
            dropdownButton(
                inputId="options-menu", ## 
                label="Options", icon=NULL, circle=FALSE,
                actionLink(inputId="option1", label="Option1"),
                actionLink(inputId="option2", label="Option2")
            )
        )
        ## dropdownMenu(
        ##     type = "notifications",
        ##     notificationItem("Today is a good day!"),
        ##     notificationItem("Sale!")
        ## )
    ),
    sidebar = dashboardSidebar(
        sidebarMenuOutput("sidebarmenu")
    ),
    rightsidebar = rightSidebar(
        background = "light",
        .items = tagList(
            conditionalPanel("input.menuitem=='load'", LoadingInputs("load")),
            conditionalPanel("input.menuitem=='view'", DataViewInputs("view")),
            conditionalPanel("input.menuitem=='clust'", ClusteringInputs("clust")),
            conditionalPanel("input.menuitem=='expr'", ExpressionInputs("expr")),
            conditionalPanel("input.menuitem=='enrich'", EnrichmentInputs("enrich")),
            conditionalPanel("input.menuitem=='isect'", IntersectionInputs("isect")),
            conditionalPanel("input.menuitem=='func'", FunctionalInputs("func")),
            conditionalPanel("input.menuitem=='sig'", SignatureInputs("sig")),
            conditionalPanel("input.menuitem=='bio'", BiomarkerInputs("bio")),
            conditionalPanel("input.menuitem=='prof'", ProfilingInputs("prof"))
        )
    ),
    body = dashboardBody(
        shinyjs::useShinyjs(),        
        ##useShinydashboard(),
        useShinydashboardPlus(),        
        use_waiter(include_js = FALSE), # do not include js
        includeCSS("www/playground.css"),
        extendShinyjs(text = 'shinyjs.hideControlSidebar = function(params) { $("body").removeClass("control-sidebar-open"); $(window).trigger("resize"); }'),
        extendShinyjs(text='shinyjs.showControlSidebar = function(params) { $("body").addClass("control-sidebar-open"); $(window).trigger("resize"); }'),
        htmlOutput("current_dataset"),
        tabItems(
            tabItem(tabName = "load", LoadingUI("load")),
            tabItem(tabName = "view", DataViewUI("view")),
            tabItem(tabName = "clust", ClusteringUI("clust")),
            tabItem(tabName = "expr", ExpressionUI("expr")),
            tabItem(tabName = "enrich", EnrichmentUI("enrich")),
            tabItem(tabName = "isect", IntersectionUI("isect")),
            tabItem(tabName = "func", FunctionalUI("func")),
            tabItem(tabName = "sig", SignatureUI("sig")),
            tabItem(tabName = "bio", BiomarkerUI("bio")),
            tabItem(tabName = "prof", ProfilingUI("prof"))
        )
    ),    
    footer = tagList(
        dashboardFooter(left_text = HTML(footer.txt)),
        ## social_buttons(),
        show_waiter_on_load(spin_fading_circles()) # place at the bottom
    )
)

shiny::shinyApp(ui, server)
