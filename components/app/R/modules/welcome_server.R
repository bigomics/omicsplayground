##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

WelcomeBoard <- function(id, auth )
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
          first.name <- strsplit("ivo kwee",split="[@ .]")[[1]][1]
          first.name <- paste0(toupper(substring(first.name,1,1)),
                               substring(first.name,2,nchar(first.name)))
          welcome <- paste0("Welcome back ",first.name,"...")
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
