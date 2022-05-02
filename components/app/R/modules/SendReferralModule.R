##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

if(0) {
  ## install.packages("shinyFeedback")

  shiny::shinyApp(
    ui = shiny::fluidPage(
      shiny::actionButton("show","show"),
      send_referral_ui("referral") 
    ),
    server = function(input, output) {
      send_referral_server(
        id = "referral",
        r.show = reactive(input$show)
      )
    }
  )


}

SendReferralModuleUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("modal"))
}

SendReferralModule <- function(id, r.show=reactive(0), close.callback="")
{
  shiny::moduleServer(id, function(input, output, session)
  {
    ns <- session$ns
    rv <- reactiveValues(
        success = 1,
        emails  = c()
    )
    
    ## JS logout callback
    js.cb = "function(x){logout();}"

    ## -------------------- modal UI --------------------------
    email_modal <- eventReactive( r.show(), {
      
      if(r.show()==0) return()
      
      shiny::showModal(
        shiny::modalDialog(
          title = "Session Expired",
          size = "l",
          shinyFeedback::useShinyFeedback(),
          p("Please enter three email addresses."),
          fluidRow(
            column(
              4,
              textInput(
                ns("name1"),
                "Name"
              ),
              textInput(
                ns("email1"),
                "Email"
              )
            ),
            column(
              4,
              textInput(
                ns("name2"),
                "Name"
              ),
              textInput(
                ns("email2"),
                "Email"
              )
            ),
            column(
              4,
              textInput(
                ns("name3"),
                "Name"
              ),
              textInput(
                ns("email3"),
                "Email"
              )
            )
          ),
          p(
            class = "text-danger text-center",
            id = "referral-global-error"
          ),
          footer = fillRow(
            flex = c(NA,1,NA),
            actionButton(
              ns("close"),
              "Maybe later...",
              class = "btn btn-warning"
              ##icon = icon("times")
            ),
            br(),
            actionButton(
              ns("sendRefs"),
              "Send emails",
              class = "btn btn-primary",
              icon = icon("paper-plane")
            )
#            tags$a(
#              class = "btn btn-danger",
#              icon("times"),
#              "Close",
#              onClick = HTML(close.callback)
#            )
          )
        )
      )  ## end of showModal
    })

    observeEvent(input$close,{
      removeModal()
      rv$success <- 0
      rv$emails <- c()
    })
        
    output$modal <- shiny::renderUI({
      email_modal()
    }) ## end of output$modal
    
    ## react if send button if pressed ----------------------------    
    shiny::observeEvent( input$sendRefs,
    {
      
      message("[observeEvent:input$sendRefs] reacted!")

      message("[observeEvent:input$sendRefs] name1 = ",input$name1)
      message("[observeEvent:input$sendRefs] name2 = ",input$name2)
      message("[observeEvent:input$sendRefs] name3 = ",input$name3)            

      message("[observeEvent:input$sendRefs] email1 = ",input$email1)
      message("[observeEvent:input$sendRefs] email2 = ",input$email2)
      message("[observeEvent:input$sendRefs] email3 = ",input$email3)            
      
      ## check inputs
      input_errors <- FALSE
      
      emails <- trimws(
        c(
          email1 = input$email1,
          email2 = input$email2,
          email3 = input$email3
        )
      )
      
      ## check emails
      check_email <- function(e) {
        k <- match(e,names(emails))
        shinyFeedback::hideFeedback(e)
        if(emails[e] == "") {
          shinyFeedback::showFeedbackWarning(inputId = e ,"Missing email")
          input_errors <<- TRUE        
        } else if(!grepl("@",emails[e])) {
          shinyFeedback::showFeedbackWarning(inputId = e ,"Invalid email")
          input_errors <<- TRUE        
        } else if(duplicated(emails)[k]) {
          shinyFeedback::showFeedbackWarning(inputId = e ,"Duplicated email")
          input_errors <<- TRUE        
        }
      }
      
      check_email("email1")
      check_email("email2")
      check_email("email3")
            
      # check names
      check_name <- function(e) {
        shinyFeedback::hideFeedback(e)
        if(input[[e]] == "") {
          shinyFeedback::showFeedbackWarning(inputId = e ,"Missing name")
          input_errors <<- TRUE        
        }
      }
      
      check_name("name1")
      check_name("name2")
      check_name("name3")
      
      if(input_errors)
        return()

      ## send emails
      body <- list(
        referrer = "The user",
        referrals = list(
          list(
            name = input$name1,
            email = input$email1
          ),
          list(
            name = input$name2,
            email = input$email2
          ),
          list(
            name = input$name3,
            email = input$email3
          )
        )
      )
      all_good = TRUE

      if(1) {
        token <- Sys.getenv("HONCHO_TOKEN", "")
        uri <- sprintf("%s/referral?token=%s", opt$HONCHO_URL, token)
        response <- httr::POST(
          uri,
          body = body,
          encode = "json"
        )
        
        # check response
        content <- httr::content(response)
        all_good <- lapply(content, function(ref) {
          return(ref$success)
        }) %>% 
        unlist() %>% 
          all()
      }
      
      if(!all_good) {
        session$sendCustomMessage(
          "referral-global-error", 
          list(
            message = "One or more of these email address was erroneous"
          )
        )
        return()
      }

      ## remove modal
      removeModal()
      rv$success <- rv$success + 1
      rv$emails <- emails
    }) ## observer sendRefs event


    ## return object
    list(
        success = reactive(rv$success),
        emails  = reactive(rv$emails)
    )

  }) ## moduleServer
  
}
