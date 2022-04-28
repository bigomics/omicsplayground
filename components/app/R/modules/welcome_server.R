##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

WelcomeBoard <- function(id,
                      auth
                      )
{
  moduleServer(id, function(input, output, session) 
  {
    ns <- session$ns ## NAMESPACE

    output$welcome <- shiny::renderText({
        name <- auth$name()
        dbg("[HomeBoard] name = ",name)        
        if(name %in% c("",NA,NULL)) {
          welcome <- "Welcome back..."
        } else {          
          welcome <- paste0("Welcome back",name,"...")
        }
        welcome
    })

    observeEvent( input$load_data, {
      ## shinyjs::click("")
    })

    observeEvent( input$upload_new, {
      ## shinyjs::click("")
    })
    

    
  })
}
