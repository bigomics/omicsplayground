##library(shiny)
##library(shinyjs)
##library(shinyalert)

TimerUI <- function(id) {
  ns <- shiny::NS(id)
  tagList(
      shinyjs::useShinyjs(),      
      shinyalert::useShinyalert(),  # Set up shinyalert
      shiny::actionButton( ns("reset"), "Press me to reset the timer"),
      ##checkboxInput( ns("pause"), "Click me to pause the timer"),      
      br(),
      shiny::textOutput(ns("currentTime"), container=span)
  )
}

TimerModule <- function(input, output, session, callbackR,
                        exit.time = 10, exit.title = "Timed out", exit.text = "", 
                        warning.time = 5, warning.title = "Warning", warning.text = "", 
                        poll = 1)
{
    
    startTime = Sys.time()
    message("startTime = ",startTime)
    has.warned <- FALSE
    timer.on <- reactiveVal(TRUE)
    
    resetTimer <- function() {
        has.warned <<- FALSE                    
        startTime <<- Sys.time()        
        timer.on(TRUE)
    }
    
    observeEvent(input$reset, {
        resetTimer()
    })
    
    ##getTime <- reactive({ Sys.time() - startTime })    
    output$currentTime <- renderText({
        shiny::invalidateLater(poll*500)
        format(Sys.time() - startTime)
    })    

    observe({

        if(timer.on()) shiny::invalidateLater(poll*1000)
        seconds.passed <- as.numeric(Sys.time() - startTime)
        message("seconds.passed = ",seconds.passed)
        
        if(seconds.passed > exit.time && timer.on()) {
            message("10 seconds passed!")                        
            shinyalert::closeAlert()
            FUN.callbackR <- function(x) {
                resetTimer()
                callbackR()
                timer.on(TRUE)
            }
            shinyalert::shinyalert(
                            title= exit.title,
                            text = exit.text,
                            timer = 0,
                            callbackR = FUN.callbackR,
                            type = "warning")
            ##startTime  <<- Sys.time()
            has.warned <<- FALSE
            timer.on(FALSE)
            return()
        }

        if(seconds.passed > warning.time && !has.warned && timer.on()) {
            message("5 seconds passed!")
            has.warned <<- TRUE
            shinyalert::closeAlert()
            shinyalert::shinyalert(
                            title = warning.title,
                            text = warning.text,
                            ##timer = 2000,                            
                            type="warning")
            return()            
        }
    })
}



if(0) {

    ui <- fluidPage(
        shinyjs::useShinyjs(),
        TimerUI("test")
    )

    server <- function(input, output, session) {
        callbackR <- function(x) {
            message("time out! callbackR called!!")
            stopApp()
        }
        callModule(TimerModule, "test", callbackR,
                   exit.time = 20,
                   exit.title = "Oh No!",
                   exit.text = "Your FREE session has timed out.",                
                   warning.time = 10,
                   warning.title = "Warning",
                   warning.text = "Your FREE session will expire in 10 seconds.",
                   poll = 1) 
        
    }
    shinyApp(ui = ui, server = server)

}
