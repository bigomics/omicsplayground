##library(shiny)
##library(shinyjs)
##library(shinyalert)

TimerModuleUI <- function(id) {
    ns <- shiny::NS(id)
}

TimerModule <- function(id,
                        timeout,
                        grace_time = 0,
                        max_grace = 3,
                        reset = reactive(0),
                        timeout.callback = NULL,
                        grace.callback = NULL,                         
                        poll = 1)
{
    moduleServer(id, function(input, output, session)
    {

        start_time <- Sys.time()
        message("[timer_module] start_time = ",start_time)
        
        ngrace = 0
        warn.time = timeout - max_grace * grace_time
        message("[timer_module] warn.time = ",warn.time)
        if(warn.time < 0) {
            stop("invalid and max_grace and grace_time")
        }
        
        observeEvent( reset(), {
            if(reset()==0) return()            
            start_time <<- Sys.time()        
        })

        lapse_time <- function() {
            lapse_time <- difftime(Sys.time(), start_time, units = "secs")
            lapse_time
        }

        timer <- reactive({
            shiny::invalidateLater(poll*1000)
            lapse_time()
        })    

        timeout.trigger <- reactive({
            is_lapsed <- (lapse_time() > timeout)
            if(is_lapsed) {
                message("**** TIME OUT! time out at = ",lapse_time())                
                if(!is.null(timeout.callback)) timeout.callback()                
                shiny::invalidateLater(Inf)
            } else {
                ##shiny::invalidateLater(poll*1000)
                shiny::invalidateLater(timeout*1000)                
            }
            is_lapsed
        })

        grace.trigger <- reactive({
            is_warned <- (lapse_time() > warn.time)
            is_lapsed <- (lapse_time() > timeout)            
            if(is_warned && ngrace < max_grace) {
                message("**** giving grace! lapse_time = ",lapse_time())
                message("**** ngrace = ",ngrace)                                
                shiny::invalidateLater(grace_time*1000)
                if(!is.null(grace.callback)) grace.callback()
                ngrace <<- ngrace + 1
            } else if(is_warned && ngrace == max_grace) {
                message("sorry no grace left")                
                shiny::invalidateLater(Inf)
            } else {
                shiny::invalidateLater(warn.time*1000)
            }
            (ngrace * !is_lapsed)
        })
                
        list(
            timer    = timer,
            grace    = grace.trigger,
            timeout  = timeout.trigger
        )
    })
}


if(FALSE) {
    
    require(shiny)
    shinyApp(
        ui = fluidPage(
            textOutput("timer"),
            textOutput("status1"),            
            textOutput("status2"),
            uiOutput("warning"),
            uiOutput("timeout")
        ),
        server = function(input, output, session) {
            
            TimerModule(
                "timer", 
                timeout = 10,
                grace_time = 3,
                max_grace = 0,
                poll = 1,
                timeout.callback = function() {message("Callback: TIME OUT!!!!!!")},
                grace.callback = function() {message("Callback: Gracing you...")}                 
            ) -> tm

            social_modal_server("social", r.show = tm$grace)
            
            output$timeout <- renderUI({
                if(!tm$timeout()) return(NULL)
                showModal( modalDialog("Sorry, time's up friend!"))
            })
            output$warning <- renderUI({
                if(!tm$grace()) return(NULL)
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
                paste("grace = ",tm$grace())
            })
        }
    )
    

}
