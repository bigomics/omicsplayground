##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##


ui.inviteModal <- function(id) {

  body <- tagList(
    HTML("<center><h3><b>and earn some swag!</b></h3><p><p>Invite your friends to Omics Playground and earn some Bigomics swag including cool stickers, a 'Friendly Monster' T-Shirt or one of our awesome sustainable Dopper water bottles. Read more <a href='https://bigomics.ch/swag'><u>here</u></a>."),
    HTML("<br><br><p>Enter your friend's email:"),
    div( shiny::textInput(
        inputId = paste0(id,"_email"),
        label = "", placeholder = "Email address..."),
        style = "margin-top: -30px;"),
    br(),
    shiny::actionButton(paste0(id,"_button"),"Invite!", class="btn btn-primary" )
  )          

  modal <- shiny::modalDialog(
    title = NULL,
    bsutils::modalHeader(
      div(class="modal-title", "Share the Love. Invite A Friend."),
      style = "background-color: #f0f9fd;"),
    body,
    footer = NULL,
    size = "l",
    easyClose = TRUE,
    tags$style(".modal-dialog {width: 720px;}"),
    tags$style(".modal-content {background-color: #f0f9fd;}"),    
    tags$style(".modal-header {padding: 0px;}")    
  )
  ##return(div(id=id, modal))
  return(shiny::showModal(modal))
}


sendInviteEmail <- function(user_email, user_name, friend_email, 
                            path_to_creds = "gmail_creds")
{
  dbg("[sendInviteEmail] reacted! : path_to_creds =", path_to_creds)
  
  if (!file.exists(path_to_creds)) {
    message("[sendInviteEmail] WARNING : mail not sent. cannot get mail creds =", path_to_creds)
    return(NULL)
  }

  dbg("[sendInviteEmail] user_name =", user_name)
  dbg("[sendInviteEmail] user_email =", user_email)
  dbg("[sendInviteEmail] friend_email =", friend_email)
  
  user_email   <- trimws(user_email)
  user_name    <- trimws(user_name)
  friend_email <- trimws(friend_email)

  body_msg <- "Hi. I always thought omics analysis was so difficult, but now I am using BigOmics Playground to analyze my own omics data. No coding required. It's so easy and fun! You should really try it!"
  
  if(is.null(user_name) || is.na(user_name) || user_name == "") user_name <- user_email

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
    cc = "support@bigomics.ch",    
    subject = paste("A friend invited you to Omics Playground"),
    credentials = blastula::creds_file(path_to_creds)
  )

  ## record the invite
  invite.file <- file.path(ETC, "INVITES.log") 
  do.append <- file.exists(invite.file)
  timestamp <- as.character(Sys.time())
  invite_data <- list( "2024-03-25 14:06:18", "from.me@test.com", "to.friend@test.com")
  invite_data <- list(timestamp, user_email, friend_email)
  data.table::fwrite(invite_data, file = invite.file, quote = TRUE, append = do.append)

  
}
