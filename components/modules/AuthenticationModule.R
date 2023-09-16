##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

AuthenticationUI <- function(id) {
  ns <- shiny::NS(id) ## namespace
}

NoAuthenticationModule <- function(id,
                                   show_modal = TRUE,
                                   username = "",
                                   email = "") {
  shiny::moduleServer(
    id, function(input, output, session) {
      message("[NoAuthenticationModule] >>>> no authentication <<<<")
      ns <- session$ns
      USER <- shiny::reactiveValues(
        method = "none",
        logged = FALSE,
        username = "",
        email = "",
        level = "",
        limit = "",
        options = opt, ## init from global
        user_dir = PGX.DIR ## global
      )

      m <- splashLoginModal(
        ns = ns,
        with.username = FALSE,
        with.email = FALSE,
        with.password = FALSE,
        title = "Sign in",
        subtitle = "Ready to explore your data?",
        button.text = "Sure I am!"
      )
      shiny::showModal(m)

      resetUSER <- function() {
        USER$logged <- FALSE
        USER$username <- ""
        USER$email <- ""
        USER$level <- ""
        USER$limit <- ""
        if (show_modal) {
          shiny::showModal(m)
        } else {
          USER$logged <- TRUE
        }
        USER$username <- username
        USER$email <- email
      }

      output$showLogin <- shiny::renderUI({
        resetUSER()
      })

      output$login_warning <- shiny::renderText("")

      shiny::observeEvent(
        {
          input$login_btn
        },
        {
          shiny::removeModal()
          USER$logged <- TRUE

          ## set options. NEED RETHINK!!! should we allow USERDIR???
          ## should we allow OPTIONS???
          USER$options <- read_user_options(PGX.DIR)

          if (!is.null(username) && username != "") {
            USER$username <- username
            USER$email <- email
            ## dbg("[NoAuthenticationModule] setting username=",username)
          }
        }
      )

      ## export 'public' function
      USER$resetUSER <- resetUSER

      return(USER)
    } ## end-of-server
  )
}

AuthenticationModuleApacheCookie <- function(id,
                                   show_modal = TRUE,
                                   username = "",
                                   email = "") {
  shiny::moduleServer(
    id, function(input, output, session) {
      message("[NoAuthenticationModule] >>>> no authentication -- reading user from apache cookie <<<<")
      email <- extract_cookie_value(session$request$HTTP_COOKIE, "user")
      ns <- session$ns
      USER <- shiny::reactiveValues(
        method = "none",
        logged = FALSE,
        username = "",
        email = email,
        level = "",
        limit = "",
        options = opt, ## init from global
        user_dir = PGX.DIR ## global
      )

      m <- splashLoginModal(
        ns = ns,
        with.username = FALSE,
        with.email = FALSE,
        with.password = FALSE,
        title = "Sign in",
        subtitle = "Ready to explore your data?",
        button.text = "Sure I am!"
      )
      shiny::showModal(m)

      resetUSER <- function() {
        USER$logged <- FALSE
        USER$username <- ""
        USER$email <- ""
        USER$level <- ""
        USER$limit <- ""
        if (show_modal) {
          shiny::showModal(m)
        } else {
          USER$logged <- TRUE
        }
        USER$username <- username
        USER$email <- email
      }

      output$showLogin <- shiny::renderUI({
        resetUSER()
      })

      output$login_warning <- shiny::renderText("")

      shiny::observeEvent(input$login_btn, {
        shiny::removeModal()
        USER$logged <- TRUE

        # set options
        USER$options <- read_user_options(PGX.DIR)
      })

      ## export 'public' function
      USER$resetUSER <- resetUSER

      return(USER)
    } ## end-of-server
  )
}

## ================================================================================
## FirebaseAuthenticationModule
## ================================================================================

FirebaseAuthenticationModule <- function(id,
                                         domain = NULL,
                                         credentials_file = NULL,
                                         firebase.rds = "firebase.rds",
                                         allow_personal = TRUE,
                                         allow_new_users = TRUE) {
  shiny::moduleServer(id, function(input, output, session) {
    message("[AuthenticationModule] >>>> using FireBase authentication <<<<")

    if (file.exists(firebase.rds)) {
      firebase_config <- firebase:::read_config(firebase.rds)
    } else {
      stop("[FATAL ERROR] no firebase.rds file found. please create.")
    }
    Sys.setenv(OMICS_GOOGLE_PROJECT = firebase_config$projectId)

    ns <- session$ns
    USER <- shiny::reactiveValues(
      method = "firebase",
      logged = FALSE,
      username = "",
      password = "",
      email = "",
      level = "",
      limit = "",
      token = NULL,
      uid = NULL,
      stripe_id = NULL,
      href = NULL,
      options = opt,
      user_dir = PGX.DIR
    )

    firebase <- firebase::FirebaseSocial$
      new(persistence = "local")

    firebase2 <- firebase::FirebaseEmailLink$
      new(persistence = "local")

    shinyjs::runjs("logout()") ## comment out for login persistence

    observeEvent(input$launchGoogle, {
      firebase$launch_google(flow = "popup")
    })

    resetUSER <- function() {
      USER$logged <- FALSE
      USER$username <- ""
      USER$password <- ""
      USER$email <- ""
      USER$level <- ""
      USER$limit <- ""
      USER$token <- ""

      ## sign out (THIS LOOSES PERSISTENCE!)
      firebase$sign_out()
      dbg("[FirebaseAuthenticationModule:resetUSER] *** signing out of firebase **** ")

      m <- splashLoginModal(
        ns = ns,
        with.email = FALSE,
        with.username = FALSE,
        with.password = FALSE,
        with.firebase = TRUE,
        with.firebase_emailonly = FALSE,
        title = HTML("Register <div style='font-size:0.5em;'>or</div> Sign in"),
        subtitle = "Enter your email and we'll send you a magic link for password-free access.",
        button.text = NULL
      )

      ## login modal
      shiny::showModal(m)
    }

    ### I don't really understand this???... (IK 10jul23)
    first_time <- TRUE
    observeEvent(USER$logged, {
      ## no need to show the modal if the user is logged this is due
      ## to persistence. But if it is the first time of the session
      ## we force reset/logout to delete sleeping (persistent?) logins.
      if (USER$logged && !first_time) {
        ## create user_dir and set correct user_dir path
        USER$user_dir <- file.path(PGX.DIR, USER$email)
        create_user_dir_if_needed(USER$user_dir, PGX.DIR)
        if (!opt$ENABLE_USERDIR) {
          USER$user_dir <- file.path(PGX.DIR)
        }
        # set options
        USER$options <- read_user_options(USER$user_dir)
        return()
      }
      first_time <<- FALSE
      resetUSER()
    })


    email_waiter <- waiter::Waiter$new(
      id = ns("emailSubmit"),
      html = div(waiter::spin_3(),
        style = "transform: scale(0.6);"
      ),
      color = waiter::transparent(.8)
    )

    observeEvent(input$emailSubmit, {
      shiny::req(input$emailInput)

      ## >>> We could check here for email validaty and intercept the
      ## login process for not authorized people with wrong domain
      ## or against a subscription list.
      email <- tolower(input$emailInput)
      check <- checkEmail(
        email = email,
        domain = domain,
        credentials_file = credentials_file,
        check.personal = !allow_personal,
        check.existing = !allow_new_users
      )

      if (!check$valid) {
        js.emailFeedbackMessage(session, check$msg, "error")
        shiny::updateTextInput(session, "emailInput", value = "")
        return(NULL)
      }

      ## >>> OK let's send auth request
      js.emailFeedbackMessage(session, "Email sent, check your inbox.", "success")
      email_waiter$show()
      firebase2$send_email(email)
      email_waiter$hide()
    })

    observeEvent(firebase$get_signed_in(), {
      response <- firebase$get_signed_in()

      if (!response$success) {
        info("[FirebaseAuthenticationModule] WARNING : Firebase sign in NOT succesful")
        resetUSER()
        return()
      }

      ## even if the response is fine, we still need to check against
      ## the allowed domain or CREDENTIALS list again, especially if
      ## the user used the social buttons to login
      user_email <- tolower(response$response$email)
      check2 <- checkEmail(
        email = user_email,
        domain = domain,
        credentials_file = credentials_file,
        check.personal = !allow_personal,
        check.existing = !allow_new_users
      )

      if (!check2$valid) {
        shinyalert::shinyalert(
          title = "",
          text = paste("Sorry.", check2$msg),
          size = "xs"
        )
        ## js.emailFeedbackMessage(session, check2$msg, "error")
        return(NULL)
      }

      on.exit({
        if (USER$logged) removeModal()
      })

      USER$logged <- TRUE
      USER$uid <- as.character(response$response$uid)
      USER$username <- response$response$displayName
      USER$email <- user_email

      if (!is.null(USER$username)) USER$username <- as.character(USER$username)
      if (!is.null(USER$email)) USER$email <- as.character(USER$email)
      if (is.null(USER$username)) USER$username <- ""
      if (is.null(USER$email)) USER$email <- ""

      # create user_dir (always), set path, and set options
      USER$user_dir <- file.path(PGX.DIR, USER$email)
      create_user_dir_if_needed(USER$user_dir, PGX.DIR)
      if (!opt$ENABLE_USERDIR) {
        USER$user_dir <- file.path(PGX.DIR)
      }
      USER$options <- read_user_options(USER$user_dir)
      session$sendCustomMessage("get-permissions", list(ns = ns(NULL)))
    })

    observeEvent(input$stripeId, {
      USER$stripe_id <- input$stripeId$id
      USER$href <- input$stripeId$href
    })

    observeEvent(input$permissions, {
      perm <- input$permissions

      USER$level <- "free"
      if (perm$success) {
        USER$level <- "premium"
      }

      session$sendCustomMessage(
        "set-user",
        list(
          user = USER$email,
          level = USER$level,
          pricing = Sys.getenv("OMICS_STRIPE_PREMIUM_PRICE")
        )
      )
    })

    observeEvent(input$firebaseUpgrade, {
      current.plan <- USER$level
      showModal(
        upgrade.dialog(ns, current.plan)
      )
    })

    observeEvent(input$manage, {
      response <- httr::POST(
        "https://api.stripe.com/v1/billing_portal/sessions",
        body = list(
          customer = USER$stripe_id,
          return_url = USER$href
        ),
        httr::authenticate(
          Sys.getenv("OMICS_STRIPE_KEY"),
          ""
        ),
        encode = "form"
      )

      httr::warn_for_status(response)
      content <- httr::content(response)
      session$sendCustomMessage("manage-sub", content$url)
    })

    observeEvent(input$launchFacebook, {
      shinyalert::shinyalert(
        title = "",
        text = "Sorry. Facebook login is not yet avaible.",
        size = "xs"
      )
    })

    observeEvent(input$launchApple, {
      shinyalert::shinyalert(
        title = "",
        text = "Sorry. Apple login is not yet avaible.",
        size = "xs"
      )
    })

    observeEvent(input$launchMicrosoft, {
      shinyalert::shinyalert(
        title = "",
        text = "Sorry. Microsoft login is not yet avaible.",
        size = "xs"
      )
    })

    observeEvent(input$launchTwitter, {
      shinyalert::shinyalert(
        title = "",
        text = "Sorry. Twitter login is not yet avaible.",
        size = "xs"
      )
    })

    ## export 'public' functions
    USER$resetUSER <- resetUSER

    return(USER)
  })
}

## ================================================================================
## EmailAuthenticationModule (using Firebase, no stripe!!!)
## ================================================================================

EmailLinkAuthenticationModule <- function(id,
                                          pgx_dir,
                                          domain = NULL,
                                          credentials_file = NULL,
                                          allow_new_users = TRUE,
                                          allow_personal = TRUE,
                                          firebase.rds = "firebase.rds") {
  shiny::moduleServer(id, function(input, output, session) {
    message("[EmailLinkAuthenticationModule] >>>> using email link (Firebase) authentication <<<<")

    if (file.exists(firebase.rds)) {
      firebase_config <- firebase:::read_config(firebase.rds)
    } else {
      stop("[EmailLinkAuthenticationModule] FATAL ERROR : no firebase.rds file")
    }
    Sys.setenv(OMICS_GOOGLE_PROJECT = firebase_config$projectId)
    if (!is.null(credentials_file) && credentials_file == FALSE) credentials_file <- NULL

    ns <- session$ns
    USER <- shiny::reactiveValues(
      method = "email-link",
      logged = FALSE,
      username = "",
      password = "",
      email = "",
      level = "",
      limit = "",
      token = "",
      uid = NA,
      stripe_id = NA,
      href = NA,
      options = opt,
      user_dir = PGX.DIR
    )

    firebase <- firebase::FirebaseSocial$
      new(persistence = "local")

    firebase2 <- firebase::FirebaseEmailLink$
      new(persistence = "local")

    observeEvent(input$launchGoogle, {
      firebase$launch_google(flow = "popup")
    })

    shinyjs::runjs("logout()") ## comment out for login persistence

    resetUSER <- function() {
      message("[EmailLinkAuthenticationModule] resetting USER... ")

      USER$logged <- FALSE
      USER$username <- ""
      USER$password <- ""
      USER$email <- ""
      USER$level <- ""
      USER$limit <- ""
      USER$token <- ""

      ## sign out (THIS LOOSES PERSISTENCE!)
      firebase$sign_out()
      dbg("[EmailLinkAuthenticationModule:resetUSER] *** signing out of firebase **** ")

      title <- HTML("Sign up <div style='font-size:0.4em;'>or</div> Log in")
      if (!is.null(credentials_file) && file.exists(credentials_file)) {
        title <- "Log in"
      }

      m <- splashLoginModal(
        ns = ns,
        with.email = FALSE,
        with.username = FALSE,
        with.password = FALSE,
        with.firebase = FALSE,
        with.firebase_emailonly = TRUE,
        title = title,
        subtitle = "Enter your email and we'll send you a magic link for password-free access.",
        button.text = "Send!"
      )

      ## login modal
      shiny::showModal(m)
    }

    first_time <- TRUE
    observeEvent(USER$logged, {
      ## no need to show the modal if the user is logged this is due
      ## to persistence. But if it is the first time of the session
      ## we force reset/logout to delete sleeping logins.
      if (USER$logged && !first_time) {
        return()
      }
      first_time <<- FALSE
      resetUSER()
    })

    email_waiter <- waiter::Waiter$new(
      id = ns("emailSubmit"),
      html = div(waiter::spin_3(),
        style = "transform: scale(0.6);"
      ),
      color = waiter::transparent(.8)
    )

    observeEvent(input$emailSubmit, {
      ## >>> We could check here for email validaty and intercept the
      ## login process for not authorized people with wrong domain
      email <- tolower(input$emailInput)
      check <- checkEmail(
        email = email,
        domain = domain,
        credentials_file = credentials_file,
        check.personal = !allow_personal,
        check.existing = !allow_new_users
      )

      if (!check$valid) {
        js.emailFeedbackMessage(session, check$msg, "error")
        shiny::updateTextInput(session, "emailInput", value = "")
        return(NULL)
      }

      ## >>> OK let's send auth request
      js.emailFeedbackMessage(session, "Email sent, check your inbox.", "success")
      email_waiter$show()
      firebase2$send_email(email)
      email_waiter$hide()
    })

    observeEvent(firebase$get_signed_in(), {
      response <- firebase$get_signed_in()

      if (!response$success) {
        info("[EmailLinkAuthenticationModule] WARNING : Firebase signin NOT succesful")
        resetUSER()
        return()
      }

      on.exit({
        if (USER$logged) removeModal()
      })

      USER$logged <- TRUE
      USER$uid <- as.character(response$response$uid)
      USER$username <- response$response$displayName
      USER$email <- response$response$email

      if (!is.null(USER$username)) USER$username <- as.character(USER$username)
      if (!is.null(USER$email)) USER$email <- as.character(USER$email)
      if (is.null(USER$username)) USER$username <- ""
      if (is.null(USER$email)) USER$email <- ""

      # create user_dir (always), set path, and set options
      USER$user_dir <- file.path(PGX.DIR, USER$email)
      create_user_dir_if_needed(USER$user_dir, PGX.DIR)
      if (!opt$ENABLE_USERDIR) {
        USER$user_dir <- file.path(PGX.DIR)
      }
      USER$options <- read_user_options(USER$user_dir)

      session$sendCustomMessage("set-user", list(user = USER$email))
      session$sendCustomMessage("get-permissions", list(ns = ns(NULL)))
    })

    ## export 'public' functions
    USER$resetUSER <- resetUSER

    return(USER)
  })
}

## ================================================================================
## PasswordAuthenticationModule (ask login.name + password)
## ================================================================================

PasswordAuthenticationModule <- function(id,
                                         credentials_file,
                                         allow_personal = TRUE,
                                         domain = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    message("[PasswordAuthenticationModule] >>>> using password authentication <<<<")

    ns <- session$ns
    if (!is.null(credentials_file) && credentials_file == FALSE) credentials_file <- NULL

    USER <- shiny::reactiveValues(
      method = "password",
      logged = FALSE,
      username = NA,
      email = NA,
      password = NA,
      level = "",
      limit = "",
      options = opt,
      user_dir = PGX.DIR
    )

    login_modal <- splashLoginModal(
      ns = ns,
      with.email = TRUE,
      with.username = FALSE,
      with.password = TRUE,
      title = "Log in",
      subtitle = "Ready to explore your data?",
      button.text = "Start!"
    )
    shiny::showModal(login_modal) ## need first time

    resetUSER <- function() {
      USER$logged <- FALSE
      USER$username <- NA
      USER$email <- NA
      USER$password <- NA
      USER$level <- ""
      USER$limit <- ""
      shiny::showModal(login_modal)
    }

    CREDENTIALS <- read.csv(credentials_file, colClasses = "character")

    output$showLogin <- shiny::renderUI({
      shiny::showModal(login_modal)
    })

    output$login_warning <- shiny::renderText("")

    shiny::observeEvent(input$login_btn, {
      login.OK <- FALSE
      valid.date <- FALSE
      valid.user <- FALSE

      ##      login_username <- input$login_username
      login_email <- input$login_email
      login_password <- input$login_password


      ## >>> We check here for email validaty and intercept the
      ## login process for not authorized people with wrong domain
      check <- checkEmail(
        email = login_email,
        domain = domain,
        credentials_file = credentials_file,
        check.personal = !allow_personal,
        check.existing = FALSE
      )

      if (!check$valid) {
        output$login_warning <- shiny::renderText(check$msg)
        shinyjs::delay(4000, {
          output$login_warning <- shiny::renderText("")
        })
        return(NULL)
      }

      if (is.null(login_password) || login_password == "") {
        output$login_warning <- shiny::renderText("missing password")
        shinyjs::delay(4000, {
          output$login_warning <- shiny::renderText("")
        })
        return(NULL)
      }

      sel <- which(CREDENTIALS$email == login_email)[1]
      valid.user <- isTRUE(length(sel) > 0)
      valid.pw <- isTRUE(CREDENTIALS[sel, "password"] == input$login_password)
      valid.date <- isTRUE(Sys.Date() < as.Date(CREDENTIALS[sel, "expiry"]))

      login.OK <- (valid.user && valid.pw && valid.date)

      if (login.OK) {
        message("[PasswordAuthenticationModule::login] PASSED : login OK! ")
        output$login_warning <- shiny::renderText("")
        shiny::removeModal()
        sel <- which(CREDENTIALS$email == login_email)[1]
        cred <- CREDENTIALS[sel, ]
        USER$username <- cred$username
        USER$email <- cred$email
        USER$level <- cred$level
        USER$limit <- cred$limit
        USER$logged <- TRUE

        # create user_dir (always), set path, and set options
        USER$user_dir <- file.path(PGX.DIR, USER$email)
        create_user_dir_if_needed(USER$user_dir, PGX.DIR)
        if (!opt$ENABLE_USERDIR) {
          USER$user_dir <- file.path(PGX.DIR)
        }
        USER$options <- read_user_options(USER$user_dir)

        ## need for JS hsq tracking
        session$sendCustomMessage("set-user", list(user = USER$username))
      } else {
        message("[PasswordAuthenticationModule::login] WARNING : login failed ")
        if (!valid.date) {
          output$login_warning <- shiny::renderText("Registration expired")
        }
        if (!valid.pw) {
          output$login_warning <- shiny::renderText("Invalid password")
        }
        if (!valid.user) {
          output$login_warning <- shiny::renderText("Invalid user")
        }
        shinyjs::delay(4000, {
          output$login_warning <- shiny::renderText("")
        })
        USER$logged <- FALSE
      }
    })

    ## export as 'public functions' :)
    USER$resetUSER <- resetUSER
    return(USER)
  })
}

## ================================================================================
## PasswordAuthenticationModule (ask login.name + password)
## ================================================================================

LoginCodeAuthenticationModule <- function(id,
                                          mail_creds,
                                          domain = NULL,
                                          credentials_file = NULL,
                                          allow_personal = TRUE,
                                          allow_new_users = TRUE) {
  shiny::moduleServer(id, function(input, output, session) {
    message("[AuthenticationModule] >>>> using secret authentication <<<<")

    ## user mail_creds="" for dry-run
    if (!file.exists(mail_creds)) {
      ## we continue but email is not working
      warning("[LoginCodeAuthenticationModule] ERROR : missing mail_creds file!!!")
    }
    if (!is.null(credentials_file) && credentials_file == FALSE) credentials_file <- NULL

    ns <- session$ns
    USER <- shiny::reactiveValues(
      method = "login-code",
      logged = FALSE,
      username = NA,
      email = NA,
      level = "",
      limit = "",
      options = opt, ## global
      user_dir = PGX.DIR ## global
    )

    email_sent <- FALSE
    login_code <- NULL

    login_modal <- splashLoginModal(
      ns = ns,
      with.email = TRUE,
      with.username = FALSE,
      with.password = FALSE,
      title = "Enter Email",
      subtitle = "To register or sign in, enter your email and we'll send you a login code.",
      button.text = "Send code!"
    )
    shiny::showModal(login_modal) ## need first time

    resetUSER <- function() {
      USER$logged <- FALSE
      USER$username <- NA
      USER$email <- NA
      USER$password <- NA
      USER$level <- ""
      USER$limit <- ""

      email_sent <<- FALSE
      login_code <<- NULL
      updateTextInput(session, "login_email", value = "")
      updateTextInput(session, "login_password", value = "")
      shiny::showModal(login_modal)
    }

    sendLoginCode <- function(user_email, login_code, mail_creds) {
      if (!file.exists(mail_creds)) {
        warning("[LoginCodeAuthenticationModule:sendLoginCode] WARNING: no mail_creds file")
        return(NULL)
      }
      blastula::smtp_send(
        blastula::compose_email(
          body = blastula::md(
            glue::glue(
              "Hello,",
              "<p>We received a request to sign in to Omics Playground using",
              "this email address. If you want to sign in with your",
              "{user_email} account, please use this login code:",
              "<p>{login_code}",
              "<p>If you did not request this code, you can safely ignore this email.",
              "<p>Thanks,",
              "<p>BigOmics Team",
              .sep = " "
            )
          ),
          footer = blastula::md(
            "BigOmics - Advanced omics analysis for everyone."
          )
        ),
        from = "support@bigomics.ch",
        to = user_email,
        subject = paste("Your login code to Omics Playground"),
        credentials = blastula::creds_file(mail_creds)
      )
    }

    output$showLogin <- shiny::renderUI({
      email_sent <<- FALSE
      login_code <<- NULL
      shiny::showModal(login_modal)
    })

    output$login_warning <- shiny::renderText("")

    email_waiter <- waiter::Waiter$new(
      id = ns("login_btn"),
      html = div(waiter::spin_3(),
        style = "transform: scale(0.6);"
      ),
      color = waiter::transparent(.8)
    )

    ## --------------------------------------
    ## Step 1: react on send email button
    ## --------------------------------------
    query_email <- shiny::reactive(shiny::getQueryString()$email)

    shiny::observeEvent(c(input$login_btn, query_email()),
      {
        if (is.null(query_email())) {
          shiny::req(input$login_email)
          login_email <- input$login_email
        } else {
          login_email <- query_email()
        }

        if (!email_sent) {
          login_email <- tolower(login_email)
          ## >>> We check here for email validaty and intercept the
          ## login process for not authorized people with wrong domain
          check <- checkEmail(
            email = login_email,
            domain = domain,
            credentials_file = credentials_file,
            check.personal = !allow_personal,
            check.existing = !allow_new_users
          )

          if (!check$valid) {
            output$login_warning <- shiny::renderText(check$msg)
            shinyjs::delay(4000, {
              output$login_warning <- shiny::renderText("")
            })
            return(NULL)
          }

          ## MAIL CODE TO USER
          ## login_code <- "hello123"
          ## login_code <<- paste0(sample(c(LETTERS), 6), collapse = "")
          login_code <<- paste(sapply(1:3, function(i) paste(sample(LETTERS, 4), collapse = "")), collapse = "-")

          info("[LoginCodeAuthenticationModule] sending login code", login_code, "to", login_email)
          email_waiter$show()
          sendLoginCode(login_email, login_code, mail_creds = mail_creds)
          email_waiter$hide()

          USER$email <- login_email
          USER$username <- login_email
          USER$logged <- FALSE
          email_sent <<- TRUE

          ## change buttons and field
          login_modal2 <- splashLoginModal(
            ns = ns,
            with.email = FALSE,
            with.username = FALSE,
            with.password = TRUE,
            hide.password = FALSE,
            title = "Enter Code",
            subtitle = "Enter the login code that we have just sent to you.",
            button.text = "Submit",
            add.cancel = TRUE,
            cancel.text = "Cancel"
          )
          shiny::showModal(login_modal2)
          updateTextInput(session, "login_email", value = "")
          updateTextInput(session, "login_password", value = "", placeholder = "enter code")

          shinyalert::shinyalert(
            title = "",
            text = "We have emailed you a login code. Please check your mailbox.",
            size = "xs"
          )
        }
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    )

    ## not sure why but using input$login_password directly does not
    ## work as the value does not reset for the next user (IK 8jul23)
    entered_code <- shiny::reactiveVal("")
    observeEvent(input$login_password, {
      entered_code(input$login_password)
    })

    ## --------------------------------------
    ## Step 2: react on submit CODE button
    ## --------------------------------------
    shiny::observeEvent(input$login_btn, {
      ## shiny::req(input$login_password)
      shiny::req(entered_code())

      if (email_sent) {
        # input_code <- input$login_password
        input_code <- entered_code()
        login.OK <- (input_code == login_code)

        if (!login.OK) {
          output$login_warning <- shiny::renderText("invalid code")
          shinyjs::delay(4000, {
            output$login_warning <- shiny::renderText("")
          })
          updateTextInput(session, "login_password", value = "")
          return(NULL)
        }

        if (login.OK) {
          output$login_warning <- shiny::renderText("")

          # create user_dir (always), set path, and set options
          USER$user_dir <- file.path(PGX.DIR, USER$email)
          create_user_dir_if_needed(USER$user_dir, PGX.DIR)
          if (!opt$ENABLE_USERDIR) {
            USER$user_dir <- file.path(PGX.DIR)
          }
          USER$options <- read_user_options(USER$user_dir)

          session$sendCustomMessage("set-user", list(user = USER$email))
          entered_code("") ## important for next user
          shiny::removeModal()

          USER$logged <- TRUE
        }
      }
    })

    shiny::observeEvent(input$cancel_btn, {
      resetUSER()
    })

    first_time <- TRUE
    observeEvent(USER$logged, {
      ## no need to show the modal if the user is logged this is due
      ## to persistence. But if it is the first time of the session
      ## we force reset/logout to delete sleeping (persistent?) logins.
      if (USER$logged && !first_time) {
        # create user_dir (always), set path, and set options
        USER$user_dir <- file.path(PGX.DIR, USER$email)
        create_user_dir_if_needed(USER$user_dir, PGX.DIR)
        if (!opt$ENABLE_USERDIR) {
          USER$user_dir <- file.path(PGX.DIR)
        }
        # set options
        USER$options <- read_user_options(USER$user_dir)
        return()
      }
      first_time <<- FALSE
      resetUSER()
    })

    ## export as 'public' functions
    USER$resetUSER <- resetUSER

    return(USER)
  })
}

## ================================================================================
## ================================= END OF FILE ==================================
## ================================================================================
