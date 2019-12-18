##install.packages("V8")
##install.packages("node-dev")

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(plotly)
library(shinyjs)
library(shinyBS)
library(V8)

ui <- dashboardPage(
    dashboardHeader(),
    dashboardSidebar(),
    dashboardBody(
        useShinyjs(),
        extendShinyjs(text = 'shinyjs.hideSidebar = function(params) { $("body").addClass("sidebar-collapse"); $(window).trigger("resize"); }'),
        extendShinyjs(text='shinyjs.showSidebar = function(params) { $("body").removeClass("sidebar-collapse"); $(window).trigger("resize"); }'),
        bsButton("showpanel", "Show/Hide sidebar",icon = icon("toggle-off"), type = "toggle",style = "info", value = TRUE),
        fluidRow(tabsetPanel(id='tabs',
                             tabPanel(value=1,title="Tab1"),
                             tabPanel(value=2,title="Tab2"),
                             tabPanel(value=3, title="Plot",
                                      fluidRow(
                                          column(12,
                                                 plotlyOutput('plot', height=800))))
                             )
                 )))


server <- function(input, output, session) { 
    output$plot <- renderPlotly({
        plot_ly(data = iris, x = ~Sepal.Length, y = ~Petal.Length)
    })
    observe({
        if(input$showpanel == TRUE) {
            js$showSidebar()
        }
        else {
            js$hideSidebar()
        }
    })
}

shinyApp(ui, server)


  

ui <- dashboardPagePlus(
    header = dashboardHeaderPlus(
        enable_rightsidebar = TRUE
    ),
    dashboardSidebar(),
    dashboardBody(
        useShinyjs(),
        ##extendShinyjs(text = 'shinyjs.hideSidebar = function(params) { $("body").addClass("sidebar-collapse"); $(window).trigger("resize"); }'),
        ##extendShinyjs(text='shinyjs.showSidebar = function(params) { $("body").removeClass("sidebar-collapse"); $(window).trigger("resize"); }'),
        extendShinyjs(text = 'shinyjs.hideSidebar = function(params) { $("body").removeClass("control-sidebar-open"); $(window).trigger("resize"); }'),
        extendShinyjs(text='shinyjs.showSidebar = function(params) { $("body").addClass("control-sidebar-open"); $(window).trigger("resize"); }'),
        bsButton("showpanel", "Show/Hide sidebar",icon = icon("toggle-off"), type = "toggle",style = "info", value = TRUE),
        fluidRow(tabsetPanel(id='tabs',
                             tabPanel(value=1,title="Tab1"),
                             tabPanel(value=2,title="Tab2"),
                             tabPanel(value=3, title="Plot",
                                      fluidRow(
                                          column(12,
                                                 plotlyOutput('plot', height=800))))
                             )
                 )),
    rightsidebar = rightSidebar()
)


server <- function(input, output, session) { 
    output$plot <- renderPlotly({
        plot_ly(data = iris, x = ~Sepal.Length, y = ~Petal.Length)
    })
    observe({
        if(input$showpanel == TRUE) {
            js$showSidebar()
        }
        else {
            js$hideSidebar()
        }
    })
}

shinyApp(ui, server)

