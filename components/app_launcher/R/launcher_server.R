##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


#' The application server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
launcher_server <- function(id, parent, load_example = NULL,
                            app_launchers = NULL, pgx=NULL ) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    coming_soon_alert <- function() {      
      shinyalert::shinyalert("Coming soon!",
        text = "Sorry. This feature is not available yet",
        type = "info")
    }

    dev_alert <- function() {      
      shinyalert::shinyalert("Developers only!",
        text = "Warning. Early development. Only use for testing.",
        type = "warning")
    }

    shiny::observeEvent(input$logo_click, {
      ui.showAboutModal()
    })

    ## If no dataset is loaded (X is NULL/empty), show a helpful popup only
    ## after the Qsee module itself is visible.  In the main application Qsee
    ## is a hidden parent nav panel during startup.
    ##
    ## Depend on is_visible() alone, not on input$nav: nav only changes once
    ## bigdash.selectTab() above completes its own client round trip, and
    ## chaining the popup off *that* instead of the visibility probe just
    ## adds a redundant extra hop before the popup can appear.
    check_pgx_loaded <- function() {
      X <- pgx$X
      has.pgx <- !(is.null(X) || length(X) == 0 || nrow(X) == 0 || ncol(X) == 0)
      if(has.pgx) {
        message("[launcher_server] PGX is loaded")
        return(TRUE)
      }
      shinyalert::shinyalert("No dataset",
        text = "You need first load a dataset to use this module.",
        type = "warning")
      return(FALSE)
    }
    
    ## ---------------- dashboards ------------------
    observeEvent(input$launch_playground, {
      ## playground is not yet fully a shiny module and is "preloaded"
      ## for now. Just select navtab
      bslib::nav_select("app-sidebar", "Dashboard", session=parent)
    })
    
    observeEvent(input$launch_qsee, {
      if(!isTRUE(opt$DEVMODE)) {
        coming_soon_alert()
        return(NULL)
      }
      dev_alert()
      run_app <- app_launchers[["qsee"]]
      if (!is.null(run_app) && is.function(run_app)) {
        run_app()
      } 
    })
    
    observeEvent(input$launch_across, {
      if(!isTRUE(opt$DEVMODE)) {
        coming_soon_alert()
        return(NULL)
      }
      dev_alert()
      run_app <- app_launchers[["across"]]
      if (!is.null(run_app) && is.function(run_app)) {
        run_app()
      } 
    })

    observeEvent( input$launch_mythril, {
      if(!isTRUE(opt$DEVMODE)) {
        coming_soon_alert()
        return(NULL)
      }      
      shinyalert::shinyalert("Empty!", "Please fill this stub.")
      ## run_app <- app_launchers[["mythril"]]
      ## if (!is.null(run_app) && is.function(run_app)) {
      ##   run_app()
      ## } 
    })

    ## ---------------- mini apps ------------------
    observeEvent(input$launch_prism, {
      dev_alert()
      if(!isTRUE(opt$DEVMODE)) {
        return(NULL)
      }
      bslib::nav_select("app-sidebar", "Prism", session=parent)
    })

    observeEvent(input$launch_idconvert, {
      dev_alert() 
      if(!isTRUE(opt$DEVMODE)) {
        return(NULL)
      }
      bslib::nav_select("app-sidebar", "IDconvert", session=parent)
    })

    ## Quick action: load the example dataset and jump to the Dashboard
    observeEvent(
      list(input$load_example),
    {
      if (is.null(load_example)) {
        warning("[launcher_server] !!! no load_example trigger available")
        return()
      }
      if (is.null(load_example())) {
        load_example(1)
      } else {
        load_example(load_example() + 1)
      }
      ##bslib::nav_select("app-sidebar", "Dashboard", session = parent)
    })

    ## Quick action: go to the Upload panel
    observeEvent(input$upload_new_data, {
      bslib::nav_select("app-sidebar", "Upload", session = parent)
    })

    ## Quick action: go to the Obi panel
    observeEvent(input$chat_with_obi, {
      if(check_pgx_loaded()) {
        bslib::nav_select("app-sidebar", "Copilot", session = parent)
      }
    })

    ## Quick action: go to AI Studio
    observeEvent(input$launch_studio, {
      if(check_pgx_loaded()) {
        bslib::nav_select("app-sidebar", "Studio", session = parent)
      }
    })
    
  })
}
