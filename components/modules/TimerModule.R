##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#
#
#

TimerModuleUI <- function(id) {
  ns <- shiny::NS(id)
  ## empty
}

TimerModule <- function(id,
                        condition,
                        timeout,
                        warn_before = 0,
                        max_warn = 3,
                        warn_callback = NULL,
                        timeout_callback = NULL
                        ) {
  moduleServer(id, function(input, output, session) {

    rv <- reactiveValues(reset=0, run=TRUE)
    nwarn <- 0
    start_time <- shiny::reactiveVal(Sys.time())
    warn_start <- timeout - max_warn * warn_before

    if (warn_start < 0) {
      stop("invalid and max_warn and warn_before")
    }

    observeEvent( condition(), {
      if(timeout>0 && condition()) {
        reset_timer()
        rv$run <- TRUE
      } else {
        rv$run <- FALSE
      }
    })

    reset_timer <- function() {
      start_time(Sys.time())
      nwarn <<- 0
    }
    
    observeEvent(rv$reset, {
      if (rv$reset == 0) return(NULL)
      reset_timer()
    })

    lapse_time <- function() {
      now <- Sys.time()
      difftime(now, start_time(), units = "secs")
    }

    timeout_event <- reactive({
      is_lapsed <- (lapse_time() > timeout)
      if (timeout < 0) {
        return(NULL)
      }
      delta <- NULL
      if (!rv$run) {
        delta <- Inf
      } else if (is_lapsed) {
        delta <- Inf
      } else { 
        delta <- 0.2 * timeout * 1000   ## trigger exactly at timout    
      }
      shiny::invalidateLater(delta)      
      is_lapsed
    })

    warn_event <- reactive({
      is_warned <- (lapse_time() > warn_start)
      is_lapsed <- (lapse_time() > timeout)
      if (warn_before == 0 || !rv$run ) {
        shiny::invalidateLater(Inf)        
        return(0)
      }
      if (is_warned && nwarn < max_warn) {
        shiny::invalidateLater(warn_before * 1000)
        nwarn <<- nwarn + 1
      } else if (is_warned && nwarn == max_warn) {
        shiny::invalidateLater(Inf)
      } else {
        ## warn_start is exactly at first warn event
        shiny::invalidateLater(warn_start * 1000)  
      }
      (nwarn * !is_lapsed)
    })

    ## ---------- built-in observers ------------------
    observeEvent(warn_event(), {
      if(warn_event()==0) return(NULL)
      if(is.null(warn_callback)) return(NULL)      
      warn_callback()
    })

    observeEvent(timeout_event(), {
      if(timeout_event()==0) return(NULL)
      if(is.null(timeout_callback)) return(NULL)      
      timeout_callback()
    })
    
    ## ---------- exported 'public functions' ---------------
    reset <- function() {
      rv$reset <- rv$reset + 1
    }
    run <- function(state=TRUE) {
      rv$run <- state
    }
    
    list(
      lapse_time = lapse_time,
      warn_event = warn_event,
      timeout_event = timeout_event,
      reset = reset,  ## function!
      run = run       ## function!
    )
  })
}

