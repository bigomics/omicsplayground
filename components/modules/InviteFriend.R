##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

InviteFriendUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::actionButton(
    ns("action"), "Invite!",
    width = "auto", class = "quick-button"
  )
}

InviteFriendModule <- function(
  id,
  auth,
  callbackR = NULL
  )
{
  moduleServer(id, function(input, output, session)
  {

    ns <- session$ns ## NAMESPACE    
    
    iv <- shinyvalidate::InputValidator$new()
    iv$add_rule("email", shinyvalidate::sv_required())
    iv$add_rule("email", shinyvalidate::sv_email())
    iv$enable()

    showModal <- function() {
      body <- tagList(
        HTML("<center><h3><b>and earn some swag!</b></h3><p><p>Invite your friends to Omics Playground and earn some Bigomics swag like cool stickers, our 'Friendly Monster' T-Shirt or one of our awesome sustainable Dopper water bottles. Read more about it <a href='https://bigomics.ch/omics-playground' target='_blank'><u>here</u></a>.<br><br>"),
        HTML("<img src='https://i0.wp.com/bigomics.ch/wp-content/uploads/2024/03/BigOmics-T-Shirt.webp?resize=76%2C86&ssl=1'><img src='https://i0.wp.com/bigomics.ch/wp-content/uploads/2024/03/BigOmics-Water-bottle.webp?resize=76%2C86&ssl=1'>"),
        HTML("<br><br><p>Enter your friend's email:"),
        div(
          shiny::textInput(
                   inputId = ns("email"),
                   label = "", placeholder = "Email address..."
                 ),
          style = "margin-top: -30px;"
        ),
        br(),
        shiny::actionButton(ns("invite"), "Invite!", class = "btn btn-primary")
      )
    
      modal <- shiny::modalDialog(
        title = NULL,
        bsutils::modalHeader(
          div(class = "modal-title", "Share the Love. Invite A Friend."),
          style = "background-color: #f0f9fd;"
        ),
        body,
        footer = NULL,
        size = "l",
        easyClose = TRUE,
        tags$style(".modal-dialog {width: 720px;}"),
        tags$style(".modal-content {background-color: #f0f9fd;}"),
        tags$style(".modal-header {padding: 0px;}")
      )

      shiny::showModal(modal)
    }

    r_click <- shiny::reactiveVal(0)
    ext_click <- function() {
      r_click(r_click() + 1)
    }

    click <- shiny::reactive({
      r_click() + input$action
    })
    
    shiny::observeEvent( click(), {
      showModal() 
    })
    
    shiny::observeEvent( input$invite, {
      invite_email <- input$email
      if(!checkValidEmailFormat(invite_email)) {
        ##shinyalert::shinyalert(text="Not a valid email")
        dbg("[observeInviteFriendButton] Not a valid email")
        return(NULL)
      }
      
      message("sending invite email to", invite_email, "\n")
      user_name <- auth$username
      user_email <- auth$email
      friend_email <- invite_email
      gmail_creds <- file.path(ETC, "gmail_creds")
      sendEmail(user_email, user_name, friend_email, path_to_creds = gmail_creds)

      shiny::removeModal()
      shiny::removeModal()      

      ## thank you modal
      if(!is.null(callbackR)) {
        dbg("[sendInviteModule] callbackR called")
        callbackR()
      } else {
        shinyalert::shinyalert(
          text = "Your friend has been invited. Thank you!",
          timer = 3000
        )
      }
      
    })

    sendEmail <- function(user_email, user_name, friend_email,
                          path_to_creds = "gmail_creds") {
      dbg("[sendInviteEmail] reacted! : path_to_creds =", path_to_creds)
      
      if (!file.exists(path_to_creds)) {
        message("[sendInviteEmail] WARNING : mail not sent. cannot get mail creds =", path_to_creds)
        return(NULL)
      }
      
      dbg("[sendInviteEmail] user_name =", user_name)
      dbg("[sendInviteEmail] user_email =", user_email)
      dbg("[sendInviteEmail] friend_email =", friend_email)
      
      user_email <- trimws(user_email)
      user_name <- trimws(user_name)
      friend_email <- trimws(friend_email)
      
      body_msg <- "Hi. I always thought omics analysis was so difficult, but now I am using BigOmics Playground to analyze my own omics data. No coding required. It's so easy and fun! You should really try it!"
      
      if (is.null(user_name) || is.na(user_name) || user_name == "") user_name <- user_email
      
      blastula::smtp_send(
        blastula::compose_email(
          body = blastula::md(
            glue::glue(
         "Hello,

          You've been invited by your friend <strong>{user_name}</strong> who
          thinks you'll like Omics Playground.

          You can read more about how Omics Playground can improve your RNAseq
          and proteomics analysis
          <a href='https://bigomics.ch/omics-playground'>here</a>.

          To create your free account, please visit the BigOmics website
          <strong>www.bigomics.ch</strong> and register.

          A note from your friend:

          \"{body_msg}\"

          Yours,

          BigOmics Team"
          )
          ),
          footer = blastula::md(
            glue::glue("Email sent on {blastula::add_readable_time()}.")
          )
        ),
        from = "bigomics.app@gmail.com",
        to = friend_email,
##        cc = "support@bigomics.ch",
        subject = paste("A friend invited you to Omics Playground"),
        credentials = blastula::creds_file(path_to_creds)
      )

      ## record the invite
      invite.file <- file.path(ETC, "INVITES.log")
      do.append <- file.exists(invite.file)
      timestamp <- as.character(Sys.time())
      invite_data <- list("2024-03-25 14:06:18", "from.me@test.com", "to.friend@test.com")
      invite_data <- list(timestamp, user_email, friend_email)
      data.table::fwrite(invite_data, file = invite.file, quote = TRUE, append = do.append)

      
    }  ## end of sendEmail


    ## return
    list(
      click = ext_click ## function!
    )
    
  }) ## end of moduleServer
}

