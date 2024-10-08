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
    callbackR = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    ## email text input validator
    iv <- shinyvalidate::InputValidator$new()
    iv$enable()

    observeEvent(input$email, {
      iv$add_rule("email", shinyvalidate::sv_required())
      iv$add_rule("email", shinyvalidate::sv_email())
    })

    showModal <- function() {
      body <- tagList(
        HTML("<center><h3><b>and earn some swag!</b></h3><p><p>Invite your friends to Omics Playground and earn some exclusive Bigomics swag like cool stickers, our 'Friendly Monster' T-Shirt or one of our awesome sustainable Dopper water bottles. Read more about it <a href='https://bigomics.ch/invite' target='_blank'><u>here</u></a>.<br><br>"),
        HTML("<img src='https://i0.wp.com/bigomics.ch/wp-content/uploads/2024/03/Stickers-03.webp?resize=76%2C86&ssl=1'>&nbsp;<img src='https://i0.wp.com/bigomics.ch/wp-content/uploads/2024/03/BigOmics-T-Shirt.webp?resize=76%2C86&ssl=1'>&nbsp;<img src='https://i0.wp.com/bigomics.ch/wp-content/uploads/2024/03/BigOmics-Water-bottle.webp?resize=76%2C86&ssl=1'>"),
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

    shiny::observeEvent(
      {
        list(r_click(), input$action)
      },
      {
        if (r_click() || input$action) {
          showModal()
        }
      }
    )

    shiny::observeEvent(input$invite, {
      friend_email <- input$email

      if (!checkValidEmailFormat(friend_email)) {
        ## shinyalert::shinyalert(text="Not a valid email")
        dbg("[observeInviteFriend] error: Not a valid email")
        return(NULL)
      }

      ## check personal email
      is_personal_email <- checkPersonalEmail(friend_email)
      if (is_personal_email) {
        shinyalert::shinyalert(text = "Please use institutional or business email")
        return(NULL)
      }

      ## check own email
      own_email <- friend_email == auth$email
      if (own_email) {
        shinyalert::shinyalert(text = "Meh... You cannot invite yourself.")
        return(NULL)
      }

      ## check already registered
      if (checkExistUserFolder(friend_email)) {
        shinyalert::shinyalert(text = "Your friend is already on Omics Playground")
        return(NULL)
      }

      ## check already invited
      invite_file2 <- file.path(auth$user_dir, "INVITES.log")
      if (!is.null(invite_file2) && file.exists(invite_file2)) {
        invite_list <- data.table::fread(invite_file2)
        colnames(invite_list) <- c("time", "from", "to")
        already_invited <- sum(invite_list$to == friend_email & invite_list$from == auth$email)
        if (already_invited > 3) {
          shinyalert::shinyalert(text = "You've already invited your friend many times!")
          return(NULL)
        }
      }

      ## Send email
      user_name <- auth$username
      user_email <- auth$email
      gmail_creds <- file.path(ETC, "gmail_creds")

      if (!file.exists(gmail_creds)) {
        shinyalert::shinyalert(
          text = "Error. Cannot connect to mail server",
          timer = 4000
        )
        return(NULL)
      }

      message("sending invite email to ", friend_email, "\n")
      sendInvitationEmail(user_email, user_name, friend_email,
        path_to_creds = gmail_creds
      )

      ## record the invite
      invite.file <- file.path(ETC, "INVITES.log")
      invite.file2 <- file.path(auth$user_dir, "INVITES.log")
      do.append <- file.exists(invite.file)
      timestamp <- as.character(Sys.time())
      invite_list <- list(timestamp, user_email, friend_email)
      data.table::fwrite(invite_list, file = invite.file, quote = TRUE, append = do.append)
      data.table::fwrite(invite_list, file = invite.file2, quote = TRUE, append = do.append)

      ## send confirmation
      sendConfirmationEmail(user_email, user_name, friend_email,
        path_to_creds = gmail_creds
      )

      ## remove modals
      shiny::removeModal()
      shiny::removeModal()

      ## thank you modal
      if (!is.null(callbackR)) {
        dbg("[sendInviteModule] callbackR called")
        callbackR()
      } else {
        shinyalert::shinyalert(
          text = "Your friend has been invited. Thank you!",
          timer = 4000
        )
      }
    })

    randomMotto <- function() {
      motto_list <- c(
        "Omics Playground. Never stop discovering.",
        "Omics Playground. Play.See.Discover.",
        "Omics Playground. Created with love by BigOmics Analytics.",
        "Omics Playground. Created in Ticino, the sunny side of Switzerland.",
        "Omics Playground. Easy but powerful.",
        "Omics Playground. Advanced omics analysis for everyone."
      )
      sample(motto_list, 1)
    }

    sendInvitationEmail <- function(user_email, user_name, friend_email,
                                    path_to_creds = "gmail_creds") {
      if (!file.exists(path_to_creds)) {
        message("[sendInviteEmail] WARNING : mail not sent. cannot get mail creds =", path_to_creds)
        return(NULL)
      }

      user_email <- trimws(user_email)
      user_name <- trimws(user_name)
      friend_email <- trimws(friend_email)

      if (is.null(user_name) || is.na(user_name) || user_name == "") user_name <- user_email

      blastula::smtp_send(
        blastula::compose_email(
          body = blastula::md(
            glue::glue(
              "Hi there,

Your friend {user_name} thinks you'd be a perfect fit for Omics Playground! We're thrilled to invite you to join our platform.

Omics Playground is an analysis and visualization cloud-based platform that is helping more than 1700 researchers worldwide to interactively explore RNA-Seq and proteomics data. By signing up, you'll gain access to more than 18 analysis modules to enhance your research.

To get started, simply click on the link below to create your account on Omics Playground:

     https://eu1.hubs.ly/H08lRhC0

Once registered, you'll be on your way to unlocking the full potential of Omics Playground and contributing to our growing community of researchers.

We can't wait to welcome you aboard!

Best,

The BigOmics Team
"
            )
          ),
          footer = blastula::md(randomMotto())
        ),
        from = "bigomics.app@gmail.com",
        to = friend_email,
        ##        cc = "support@bigomics.ch",
        subject = paste("You're invited! Join Omics Playground today"),
        credentials = blastula::creds_file(path_to_creds)
      )
    } ## end of sendInvitationEmail

    sendConfirmationEmail <- function(user_email, user_name, friend_email,
                                      path_to_creds = "gmail_creds") {
      if (!file.exists(path_to_creds)) {
        message("[sendConfirmationEmail] WARNING : mail not sent. cannot get mail creds =", path_to_creds)
        return(NULL)
      }
      user_email <- trimws(user_email)
      user_name <- trimws(user_name)
      friend_email <- trimws(friend_email)

      if (is.null(user_name) || is.na(user_name) || user_name == "") user_name <- user_email

      numref <- 1
      numsuccess <- 0
      ##      invite_log = file.path( auth$user_dir, "INVITES.log")""
      invite_file <- file.path(ETC, "INVITES.log")
      if (file.exists(invite_file)) {
        all_invites <- read.csv(invite_file, header = FALSE)
        sel <- which(all_invites[, 2] == user_email)
        all_refs <- unique(all_invites[sel, 3])
        numref <- length(all_refs)

        all_registered <- list.dirs(PGX.DIR, full.names = FALSE, recursive = FALSE)
        all_registered <- grep("@", all_registered, value = TRUE)
        numsuccess <- length(intersect(all_refs, all_registered))
      }

      blastula::smtp_send(
        blastula::compose_email(
          body = blastula::md(
            glue::glue(
              "
Dear {user_name},

Thank you for referring your friend {friend_email} to join Omics Playground!

As of now, you've referred {numref} number of colleagues of which {numsuccess} have successfully registered. You'll be one step closer to claiming your exclusive BigOmics swag!

Your referral helps our community grow and brings us closer to making omics data analysis accessible to everyone. We appreciate your support and enthusiasm for our platform!

Best,

The BigOmics Team
"
            )
          ),
          footer = blastula::md(randomMotto())
        ),
        from = "bigomics.app@gmail.com",
        to = user_email,
        ##      cc = "support@bigomics.ch",
        subject = paste("Your Invite has been sent!"),
        credentials = blastula::creds_file(path_to_creds)
      )
    } ## end of sendEmail

    ## return
    list(
      click = ext_click ## exported function!
    )
  }) ## end of moduleServer
}
