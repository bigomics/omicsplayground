## install.packages("shinydashboardPlus")
## remotes::install_github("JohnCoene/waiter")

library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjqui)
library(devtools)
require(shinyWidgets)
require(shinyBS)
library(waiter)
 
cat("===================== INIT =======================\n")

RDIR="../R/"
FILES="../lib/"

##PGX.DIR="../../playground/pgx/"
source("../R/pgx-init.R", local=TRUE)  ## pass local vars
load("../data/geiger2016-arginine.pgx")
ngs <- pgx.initialize(ngs)    

source("../R/pgx-modules.R")
source("modules/TestBoard.R")
source("modules/WordcloudBoard.R")
source("modules/ClusteringModule.R", local=TRUE)
source("modules/DataViewModule.R", local=TRUE)

USERMODE <- reactiveVal("PRO")
DEV.VERSION=1

server = function(input, output, session) {

    ##useShinydashboard()
    useShinydashboardPlus()
    js$showControlSidebar()
    cat("===================== SERVER =======================\n")
    
    inputData <- reactive({
        input$btn_load
        return(ngs)
    })
    callModule( TestBoard, "test1", inputData)
    callModule( DataViewModule, "dataview1", inputData)
    ##callModule( WordcloudBoard, "wordcloud1", inputData, ui_size)
    callModule( ClusteringModule, "clust1", inputData)

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
            menuItem("DataView", tabName = "dataview", icon=icon("table")),
            menuItem("Clustering", tabName = "clustering", icon=icon("project-diagram")),
            ##menuItem("Expression", tabName = "expression", icon=icon("bars")),
            ##menuItem("Enrichment", tabName = "enrichment", icon=icon("bars")),
            menuItem("Analysis", icon=icon("bullseye"), ##tabName = "functional", 
                     menuSubItem("KEGG", tabName = "kegg"),
                     menuSubItem("WordCloud", tabName = "wordcloud")
                     ),
            menuItem("TestBoard", tabName = "test", icon=icon("table"))
        )
        if(input$ui_usermode=="expert") {
            moreitems <- tagList(
                menuItem("Signature", tabName = "signature", icon=icon("bars")),
                menuItem("Biomarker", tabName = "biomarker", icon=icon("bars")),
                menuItem("scProfiling", tabName = "profiling", icon=icon("bars"))
            )
            menuitems <- c(menuitems, moreitems)
        }
        sidebarMenu( id="menuitem", .list = menuitems)
    })

    
    observeEvent( input$btn_load, {    
    })

    hide_waiter()
    
}

ui = dashboardPagePlus( 
    title = "Omics Playground v2",
    footer = dashboardFooter(
        left_text = HTML("Proudly designed and created by <a href='http://bigomics.ch'>BigOmics Analytics</a>."),
        ),
    header = dashboardHeaderPlus(
        title = tagList(
            span(class = "logo-lg", "Omics Playground v2"), 
            tags$img(src="bigomics-logo-white-32px.png")),
        ## enable_rightsidebar = TRUE,
        rightSidebarIcon = "gears",
        
        left_menu = tagList(
            dropdownButton(##id="menu_file", title="File",
                size="xs", label="File", icon=NULL, circle=FALSE,
                actionLink(inputId="btn_load", label="Load"),
                actionLink(inputId="btn_save", label="Save")
            ),
            dropdownBlock(
                id="menu_options", title="Options",
                radioButtons(inputId="ui_usermode", label="User mode",
                             c("basic","expert"), selected="basic",
                             inline=TRUE)
            )
        ),
        dropdownMenu(
            type = "notifications",
            notificationItem("Today is a good day!"),
            notificationItem("Sale!")
        )
        
        ),
    sidebar = dashboardSidebar(
        sidebarMenuOutput("sidebarmenu")
    ),
    rightsidebar = rightSidebar(
        background = "light",
        .items = tagList(
            conditionalPanel("input.menuitem=='functional'",
                             selectInput("input1","Input1",1:10),
                             selectInput("input2","Input2",1:10)                             
                             ),
            conditionalPanel("input.menuitem=='clustering'", ClusteringInputs("clust1")),
            conditionalPanel("input.menuitem=='dataview'", DataViewInputs("dataview1"))
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
        
        tabItems(
            tabItem(tabName = "test", TestBoardUI("test1")),
            tabItem(tabName = "dataview", DataViewUI("dataview1")),
            tabItem(tabName = "clustering", ClusteringUI("clust1")),
            ##tabItem(tabName = "tsne", uiOutput("tsneUI")),
            ##tabItem(tabName = "tsne", ClusteringBoard_tsneUI("clust1")),            
            tabItem(tabName = "wordcloud", WordcloudBoardUI("wordcloud1"))
        ),
        show_waiter_on_load(spin_fading_circles()) # place at the bottom
    )

    
)

shiny::shinyApp(ui, server)
