##library(shiny)
##library(shinyjs)
##library(shinyalert)

TimerModuleUI <- function(id) {
    ns <- shiny::NS(id)
}

TimerModule <- function(id,
                        timeout,
                        warn_before = 0,
                        max_warn = 3,
                        reset = reactive(0),
                        run = reactive(TRUE),                        
                        poll = Inf)
{
    moduleServer(id, function(input, output, session)
    {

        start_time <- shiny::reactiveVal(Sys.time())
        ##start_time <- Sys.time()
        
        nwarn = 0
        warn_start = timeout - max_warn * warn_before
        message("[timer_module] warn_start = ",warn_start)

        if(warn_start < 0) {
            stop("invalid and max_warn and warn_before")
        }
        
        observeEvent( reset(), {
            if(reset()==0) return()
            message("[TimerModule] reset at",Sys.time())
            start_time(Sys.time())
            nwarn <<- 0
        })

        lapse_time <- function() {
            difftime(Sys.time(), start_time(), units = "secs")
        }

        timer <- reactive({
            if(!run()) {            
                message("[TimerModule:timer] stopped by control")
                shiny::invalidateLater(Inf)                
            } else if(lapse_time() <= timeout) {
                shiny::invalidateLater(poll*1000)
            } else {
                message("[TimerModule] timer stopped at",Sys.time())
                shiny::invalidateLater(Inf)
            }
            lapse_time()
        })    

        timeout.trigger <- reactive({
            is_lapsed <- (lapse_time() > timeout)
            if(timeout < 0) return()
            if(!run()) {
                message("[TimerModule:timeout.trigger] stopped by control")
                shiny::invalidateLater(Inf)
            } else if(is_lapsed) {
                message("**** TIME OUT! time out at = ",lapse_time())                
                shiny::invalidateLater(Inf)
            } else {
                shiny::invalidateLater(timeout*1000)                
            }
            is_lapsed
        })

        warn.trigger <- reactive({
            is_warned <- (lapse_time() > warn_start)
            is_lapsed <- (lapse_time() > timeout)            
            if(warn_before==0) return()
            if(!run()) {
                message("[TimerModule:warn.trigger] stopped by control")                
                shiny::invalidateLater(Inf)
            } else if(is_warned && nwarn < max_warn) {
                message("**** giving warn! lapse_time = ",lapse_time())
                message("**** nwarn = ",nwarn)                                
                shiny::invalidateLater(warn_before*1000)
                nwarn <<- nwarn + 1
            } else if(is_warned && nwarn == max_warn) {
                message("sorry no warn left")                
                shiny::invalidateLater(Inf)
            } else {
                shiny::invalidateLater(warn_start*1000)
            }
            (nwarn * !is_lapsed)
        })
                
        list(
            timer    = timer,
            lapse_time = lapse_time,
            warn     = warn.trigger,
            timeout  = timeout.trigger
        )
    })
}


if(FALSE) {
    
    require(shiny)
    shinyApp(
        ui = fluidPage(
            actionButton("reset","reset"),
            checkboxInput("run","run",value=TRUE),            
            textOutput("timer"),            
            textOutput("status1"),            
            textOutput("status2"),
            uiOutput("warning"),
            uiOutput("timeout")
        ),
        server = function(input, output, session) {
            
            TimerModule(
                "timer", 
                timeout = 15,
                warn_before = 5,
                max_warn = 2,
                poll = Inf,
                reset = reactive(input$reset),
                run = reactive(input$run)                
            ) -> tm
            
            output$timeout <- renderUI({
                if(!tm$timeout()) return(NULL)
                showModal( modalDialog("Sorry, time's up friend!"))
            })

            output$warning <- renderUI({
                  if(!tm$warn()) return(NULL)
                 showModal(modalDialog("Warning time out soon..."))
            })
            
            output$timer <- renderText({
                tm$timer()               
            })

            output$status1 <- renderText({
                message("status1: reacted!!!!")
                paste("timeout = ",tm$timeout())                
            })
          
            output$status2 <- renderText({
                message("status2: reacted!!!!")
                paste("warn = ",tm$warn())
            })

        }
    )
    

}
