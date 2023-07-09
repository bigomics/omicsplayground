##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

create_or_read_user_options <- function(user_dir) {
  user_opt_file <- file.path(user_dir, "OPTIONS")
  new_opt <- opt ## opt from global
  if (!file.exists(user_opt_file)) {
    ## IK: no need
    ## file.copy(from = opt.file, to = user_opt_file)
  } else {
    user_opt <- playbase::pgx.readOptions(file = user_opt_file)
    for (opt_name in names(user_opt)) {
      new_opt[[opt_name]] <- user_opt[[opt_name]]
    }
  }
  new_opt
}

AuthenticationUI <- function(id) {
  ns <- shiny::NS(id) ## namespace
}

NoAuthenticationModule <- function(id,
                                   show_modal = TRUE,
                                   username = "",
                                   email = "") {
  shiny::moduleServer(
    id, function(input, output, session) {
      message("[NoAuthenticationModule] >>>> using no authentication <<<<")
      ns <- session$ns
      USER <- shiny::reactiveValues(
        method = "none",
        logged = FALSE,
        username = "",
        email = "",
        level = "",
        limit = "",
        options = opt
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
        USER$options <- create_or_read_user_options(PGX.DIR)
      })

      observeEvent(input$userLogout, {
        resetUSER()
      })

      return(USER)
    } ## end-of-server
  )
}


## ================================================================================
## FirebaseAuthenticationModule
## ================================================================================

upgrade.dialog <- function(ns, current.plan) {
  btn_basic <- "Go Basic!"
  btn_starter <- "Get Starter!"
  btn_premium <- "Get Premium!"
  if (current.plan == "free") btn_basic <- "Current Plan"
  if (current.plan == "starter") btn_starter <- "Current Plan"
  if (current.plan == "premium") btn_premium <- "Current Plan"

  modalDialog(
    title = h3("Find the right OmicsPlayground plan for you"),
    size = "m",
    div(
      class = "row",
      style = "padding-left:4rem;padding-right:4rem;text-align:center;",
      div(
        class = "col-md-4",
        style = "background:#F2FAFF;",
        HTML("<h4><b>Basic</b></h4>"),
        p("Try for free"),
        h3("Free!"),
        tags$ul(
          class = "list-unstyled",
          tags$li("Host up to 3 datasets"),
          tags$li("45 minutes time limit"),
          tags$li("Up to 25 samples / dataset"),
          tags$li("Up to 5 comparisons")
        ),
        shiny::actionButton(ns("get_basic"), btn_basic),
        br()
      ),
      div(
        class = "col-md-4",
        style = "background:#E8F8FF;",
        h4(HTML("<b>Starter</b>")),
        p("Great to start"),
        h3("Soon!"),
        tags$ul(
          class = "list-unstyled",
          tags$li("Host up to 10 datasets"),
          tags$li("3 hours time limit"),
          tags$li("Up to 100 samples / dataset"),
          tags$li("Up to 10 comparisons")
        ),
        shiny::actionButton(ns("get_starter"), btn_starter),
        br()
      ),
      div(
        class = "col-md-4",
        style = "background:#E2F4FF;",
        HTML("<h4><b>Premium</b></h4>"),
        p("For power users or small groups"),
        h3("Soon!"),
        tags$ul(
          class = "list-unstyled",
          tags$li("Host up to 100 datasets"),
          tags$li("8 hours time limit"),
          tags$li("Up to 2000 samples / dataset"),
          tags$li("Up to 100 comparisons")
        ),
        shiny::actionButton(ns("get_premium"), btn_premium),
        br()
      )
    ), ## content div
    div(
      style = "margin-top:3rem;text-align:center;",
      HTML("Looking for OmicsPlayground for <b>Enterprise</b>? <a href='mailto:info@bigomics.com'>Contact sales for info and pricing</a>.")
    ),
    footer = tagList(
      fillRow(
        flex = c(NA, 0.03, NA, 1, NA, NA),
        tags$label(
          class = "radio-inline",
          tags$input(
            id = "yearlyCheck",
            type = "radio",
            name = "yearly",
            onclick = "priceChange(name)",
            checked = TRUE
          ),
          "Billed yearly"
        ),
        br(),
        tags$label(
          class = "radio-inline",
          tags$input(
            id = "monthlyCheck",
            type = "radio",
            name = "monthly",
            onclick = "priceChange(name)"
          ),
          "Billed monthly"
        ),
        br(),
        shiny::actionButton(ns("manage"), "Manage Subscription"),
        modalButton("Dismiss")
      )
    )
  ) ## modalDialog
}

js.emailFeedbackMessage <- function(session, msg, type = "error") {
  session$sendCustomMessage(
    "email-feedback",
    list(
      type = type,
      msg = msg
    )
  )
}

checkAuthorizedDomain <- function(email, domain) {
  if (is.null(domain) || domain == "" || domain == "*") {
    return(TRUE)
  }
  domain1 <- strsplit(domain, split = "\\|")[[1]]
  domain1 <- paste0(paste0("@", domain1, "$"), collapse = "|")
  authorized <- grepl(domain1, email)
  authorized
}

checkAuthorizedUser <- function(email, credentials_file = NULL) {
  if (is.null(credentials_file) || credentials_file == FALSE) {
    return(TRUE)
  }
  if (!file.exists(credentials_file)) {
    return(TRUE)
  }
  CREDENTIALS <- read.csv(credentials_file, colClasses = "character")
  valid_user <- email %in% CREDENTIALS$email
  if (!valid_user) {
    return(FALSE)
  }
  sel <- match(email, CREDENTIALS$email)
  valid_date <- as.Date(CREDENTIALS$expiry[sel]) > as.Date(Sys.time())
  authorized <- valid_user && valid_date
  authorized
}

checkValidEmailFormat <- function(email) {
  grepl(".*@.*[.].*", email)
}

checkPersonalEmail <- function(email) {
  grepl("gmail|ymail|outlook|yahoo|hotmail|mail.com$|icloud|msn", email)
}

checkEmail <- function(email, domain = NULL, credentials_file = NULL, check.personal = TRUE) {
  dbg("[AuthenticationModule:checkEmail]")
  chk <- list()
  if (!checkValidEmailFormat(email)) {
    return(list(valid = FALSE, msg = "not a valid email"))
  }
  if (!checkAuthorizedDomain(email, domain)) {
    return(list(valid = FALSE, msg = "domain not authorized"))
  }
  if (!checkAuthorizedUser(email, credentials_file)) {
    return(list(valid = FALSE, msg = "not authorized user"))
  }
  if (check.personal) {
    if (heckPersonalEmail(email)) {
      return(list(valid = FALSE, msg = "No personal email allowed. Please provide your business, academic or institutional email."))
    }
  }
  list(valid = TRUE, "email ok")
}



FirebaseAuthenticationModule <- function(id,
                                         domain = NULL,
                                         credentials_file = NULL,
                                         firebase.rds = "firebase.rds") {
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
      options = opt
    )

    firebase <- firebase::FirebaseSocial$
      new(persistence = "local")

    firebase2 <- firebase::FirebaseEmailLink$
      new(persistence = "local")

    shinyjs::runjs("logout()")

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
      dbg("[FirebaseAuthenticationModule] *** signing out of firebase **** ")

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

    first_time <- TRUE
    observeEvent(USER$logged, {
      ## no need to show the modal if the user is logged this is due
      ## to persistence. But if it is the first time of the session
      ## we force reset/logout to delete sleeping logins.
      if (USER$logged && !first_time) {
        # set options
        USER$options <- create_or_read_user_options(
          file.path(PGX.DIR, USER$email)
        )
        return()
      }
      first_time <<- FALSE
      resetUSER()
    })

    observeEvent(input$userLogout, {
      message("[FirebaseAuthenticationModule] userLogout triggered!")
      resetUSER()
    })

    email_waiter <- waiter::Waiter$new(
      id = ns("emailSubmit"),
      html = div(waiter::spin_3(),
        style = "transform: scale(0.6);"
      ),
      color = waiter::transparent(.8)
    )

    sendEmailLink <- reactiveVal(NULL)

    observeEvent(input$emailSubmit, {
      shiny::req(input$emailInput)
      email_waiter$show()

      ## >>> We could check here for email validaty and intercept the
      ## login process for not authorized people with wrong domain
      ## or against a subscription list.
      check <- checkEmail(login_email, domain, credentials_file)
      if (!check$valid) {
        js.emailFeedbackMessage(session, check$msg, "error")
        sendEmailLink(NULL)
        email_waiter$hide()
        return(NULL)
      }

      ## >>> OK let's send auth request
      js.emailFeedbackMessage(session, "Email sent, check your inbox.", "success")
      sendEmailLink(input$emailInput)
    })

    observeEvent(sendEmailLink(),
      {
        firebase2$send_email(sendEmailLink())
        sendEmailLink(NULL)
        email_waiter$hide()
      },
      ignoreNULL = TRUE
    )

    observeEvent(firebase$get_signed_in(), {
      response <- firebase$get_signed_in()

      if (!response$success) {
        warning("[FirebaseAuthenticationModule] sign in NOT succesful")
        resetUSER()
        return()
      }

      ## even if the response is fine, we still need to check against
      ## the allowed domain or CREDENTIALS list again, especially if
      ## the user used the social buttons to login
      user_email <- response$response$email
      check2 <- checkEmail(user_email, domain, credentials_file)
      if (!check2$valid) {
        dbg("[Auth:observeEvent(firebase$get_signed_in()] check2$msg = ", check2$msg)
        shinyalert::shinyalert(
          title = "",
          text = paste("Sorry.", check2$msg),
          size = "xs"
        )
        ## js.emailFeedbackMessage(session, check2$msg, "error")
        return(NULL)
      }

      on.exit({
        dbg("[FirebaseAuthenticationModule] get_signed_in() on.exit")
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

      # set options
      USER$options <- create_or_read_user_options(
        file.path(PGX.DIR, USER$email)
      )

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
                                          firebase.rds = "firebase.rds") {
  shiny::moduleServer(id, function(input, output, session) {
    message("[AuthenticationModule] >>>> using email link (using Firebase) authentication <<<<")

    if (file.exists(firebase.rds)) {
      firebase_config <- firebase:::read_config(firebase.rds)
    } else {
      stop("[FATAL ERROR] no firebase.rds file found. please create.")
    }
    Sys.setenv(OMICS_GOOGLE_PROJECT = firebase_config$projectId)
    if (!is.null(credentials_file) && credentials_file == FALSE) credentials_file <- NULL

    ns <- session$ns
    USER <- shiny::reactiveValues(
      method = "email",
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
      options = opt
    )

    firebase <- firebase::FirebaseSocial$
      new(persistence = "local")

    firebase2 <- firebase::FirebaseEmailLink$
      new(persistence = "local")

    observeEvent(input$launchGoogle, {
      firebase$launch_google(flow = "popup")
    })

    shinyjs::runjs("logout()")

    resetUSER <- function() {
      message("[FirebaseAuthenticationModule] resetting USER... ")

      USER$logged <- FALSE
      USER$username <- ""
      USER$password <- ""
      USER$email <- ""
      USER$level <- ""
      USER$limit <- ""
      USER$token <- ""

      ## sign out (THIS LOOSES PERSISTENCE!)
      firebase$sign_out()
      dbg("[FirebaseAuthenticationModule] *** signing out of firebase **** ")

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
        dbg("[FirebaseAuthenticationModule] USER is already logged in! no modal")
        return()
      }

      first_time <<- FALSE
      message("[FirebaseAuthenticationModule] USER not logged in!")
      resetUSER()
    })

    observeEvent(input$userLogout, {
      message("[FirebaseAuthenticationModule] userLogout triggered!")
      resetUSER()
    })

    email_waiter <- waiter::Waiter$new(
      id = ns("emailSubmit"),
      html = div(waiter::spin_3(),
        style = "transform: scale(0.6);"
      ),
      color = waiter::transparent(.8)
    )

    sendEmailLink <- reactiveVal(NULL)

    observeEvent(input$emailSubmit, {
      email_waiter$show()

      on.exit({
        shiny::updateTextInput(session, "emailInput", value = "")
        email_waiter$hide()
      })

      if (input$emailInput == "") {
        js.emailFeedbackMessage(session, "Missing email", "error")
        return()
      }

      ## >>> We could check here for email validaty and intercept the
      ## login process for not authorized people with wrong domain
      authorized_domain <- checkAuthorizedDomain(input$emailInput, domain)
      if (!authorized_domain) {
        js.emailFeedbackMessage(session, "domain not authorized", "error")
        shiny::updateTextInput(session, "emailInput", value = "")
        return()
      }

      ## >>> We could cross-check here for valid email against a subscription list.
      authorized_user <- checkAuthorizedUser(email, credentials_file)
      if (!authorized_user) {
        js.emailFeedbackMessage(session, "user not authorized", "error")
        return()
      }

      ## if it is a new user we ask for business email, old users can go
      is_personal_email <- grepl("gmail|ymail|outlook|yahoo.com$|mail.com$|icloud.com$|msn.com$", input$emailInput)
      existing_user_dirs <- basename(list.dirs(pgx_dir))
      new_user <- !(input$emailInput %in% existing_user_dirs)
      if (is_personal_email && new_user) {
        js.emailFeedbackMessage(session, "No personal email allowed. Please use your business, academic or institutional email.", "error")
        return()
      }

      ## >>> OK let's send auth request
      js.emailFeedbackMessage(session, "Email sent, check your inbox.", "success")
      sendEmailLink(input$emailInput) ## can take a while...
    })

    observeEvent(sendEmailLink(),
      {
        firebase2$send_email(sendEmailLink())
        sendEmailLink(NULL)
        email_waiter$hide()
      },
      ignoreNULL = TRUE
    )

    observeEvent(firebase$get_signed_in(), {
      response <- firebase$get_signed_in()

      if (!response$success) {
        warning("[FirebaseAuthenticationModule] sign in NOT succesful")
        resetUSER()
        return()
      }

      on.exit({
        dbg("[FirebaseAuthenticationModule] get_signed_in() on.exit")
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

      # set options
      USER$options <- create_or_read_user_options(
        file.path(PGX.DIR, USER$email)
      )

      session$sendCustomMessage("get-permissions", list(ns = ns(NULL)))
    })

    return(USER)
  })
}

## ================================================================================
## PasswordAuthenticationModule (ask login.name + password)
## ================================================================================

PasswordAuthenticationModule <- function(id,
                                         credentials_file) {
  shiny::moduleServer(id, function(input, output, session) {
    message("[AuthenticationModule] >>>> using password authentication <<<<")

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
      options = opt
    )

    login_modal <- splashLoginModal(
      ns = ns,
      with.email = FALSE,
      with.username = TRUE,
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

      login_username <- input$login_username
      login_password <- input$login_password

      if (is.null(login_username) || login_username == "") {
        output$login_warning <- shiny::renderText("missing username")
        shinyjs::delay(2000, {
          output$login_warning <- shiny::renderText("")
        })
        return(NULL)
      }
      if (is.null(login_password) || login_password == "") {
        output$login_warning <- shiny::renderText("missing password")
        shinyjs::delay(2000, {
          output$login_warning <- shiny::renderText("")
        })
        return(NULL)
      }

      sel <- which(CREDENTIALS$username == login_username)[1]
      valid.user <- isTRUE(length(sel) > 0)
      valid.pw <- isTRUE(CREDENTIALS[sel, "password"] == input$login_password)
      valid.date <- isTRUE(Sys.Date() < as.Date(CREDENTIALS[sel, "expiry"]))
      login.OK <- (valid.user && valid.pw && valid.date)

      if (login.OK) {
        message("[PasswordAuthenticationModule::login] PASSED : login OK! ")
        output$login_warning <- shiny::renderText("")
        shiny::removeModal()
        sel <- which(CREDENTIALS$username == login_username)[1]
        cred <- CREDENTIALS[sel, ]
        USER$username <- cred$username
        USER$email <- cred$email
        USER$level <- cred$level
        USER$limit <- cred$limit
        USER$logged <- TRUE

        # set options
        USER$options <- create_or_read_user_options(
          file.path(PGX.DIR, USER$username)
        )

        session$sendCustomMessage("set-user", list(user = USER$username))
      } else {
        message("[PasswordAuthenticationModule::login] login invalid!")
        if (!valid.date) {
          output$login_warning <- shiny::renderText("Registration expired")
        }
        if (!valid.pw) {
          output$login_warning <- shiny::renderText("Invalid password")
        }
        if (!valid.user) {
          output$login_warning <- shiny::renderText("Invalid user")
        }
        shinyjs::delay(2000, {
          output$login_warning <- shiny::renderText("")
        })
        USER$logged <- FALSE
      }
    })

    observeEvent(input$userLogout, {
      resetUSER()
    })

    return(USER)
  })
}

## ================================================================================
## PasswordAuthenticationModule (ask login.name + password)
## ================================================================================

LoginCodeAuthenticationModule <- function(id,
                                          mail_creds,
                                          domain = NULL,
                                          credentials_file = NULL) {
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
      options = opt ## global
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
              "Hello,
<p>We received a request to sign in to Omics Playground using this email address. If you want to sign in with your {user_email} account, please use this login code:

<p>{login_code}

<p>If you did not request this code, you can safely ignore this email.

<p>Thanks,

<p>BigOmics Team"
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

    ## --------------------------------------
    ## Step 1: react on send email button
    ## --------------------------------------
    shiny::observeEvent(input$login_btn, {
      shiny::req(input$login_email)

      if (!email_sent) {
        login_email <- input$login_email

        ## >>> We check here for email validaty and intercept the
        ## login process for not authorized people with wrong domain
        check <- checkEmail(login_email, domain, credentials_file)
        if (!check$valid) {
          output$login_warning <- shiny::renderText(check$msg)
          shinyjs::delay(3000, {
            output$login_warning <- shiny::renderText("")
          })
          return(NULL)
        }


        ## MAIL CODE TO USER
        ## login_code <- "hello123"
        ## login_code <<- paste0(sample(c(letters,LETTERS,0:9),12),collapse='')
        login_code <<- paste0(sample(c(LETTERS), 6), collapse = "")
        dbg("[LoginCodeAuthenticationModule:observeEvent( input$login_btn] new_code = ", login_code)
        sendLoginCode(login_email, login_code, mail_creds = mail_creds)
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
    })

    ## not sure why but using input$login_password directly does not
    ## work as the value does not reset for the next user (IK 8jul23)
    entered_code <- reactiveVal("")
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
          shinyjs::delay(2000, {
            output$login_warning <- shiny::renderText("")
          })
          updateTextInput(session, "login_password", value = "")
          return(NULL)
        }

        if (login.OK) {
          message("[LoginCodeAuthenticationModule::login] 3 : login OK! ")
          output$login_warning <- shiny::renderText("")
          USER$logged <- TRUE
          USER$options <- create_or_read_user_options(
            file.path(PGX.DIR, USER$email)
          )
          session$sendCustomMessage("set-user", list(user = USER$username))
          entered_code("") ## important for next user

          shiny::removeModal()
        }
      }
    })

    shiny::observeEvent(input$cancel_btn, {
      resetUSER()
    })

    observeEvent(input$userLogout, {
      resetUSER()
    })

    return(USER)
  })
}

## ================================================================================
## UI FUNCTIONS
## ================================================================================

splashHelloModal <- function(name, msg = NULL, ns = NULL, duration = 3500) {
  if (is.null(ns)) {
    ns <- function(e) {
      return(e)
    }
  }
  message("[AuthenticationModule::splashHelloModel]")

  all.hello <- c(
    "Hello", "Salut", "Hola", "Pivet", "Ni hao", "Ciao", "Hi", "Hoi", "Hej",
    "Yassou", "Selam", "Hey", "Hei", "Grutzi", "Bonjour",
    "Namaste", "Salam", "Selamat", "Shalom", "Goeiedag", "Yaxshimusiz"
  )
  title <- paste(paste0(sample(all.hello, 3), "!"), collapse = " ")
  if (!is.null(name) && !is.na(name) && !name %in% c("NA", "")) {
    first.name <- strsplit(as.character(name), split = " ")[[1]][1]
    first.name <- paste0(
      toupper(substring(first.name, 1, 1)),
      substring(first.name, 2, 999)
    )
    title <- paste(paste0(sample(all.hello, 1), " ", first.name, "!"), collapse = " ")
  }
  subtitle <- "Have a good day!"
  subtitle <- "We wish you many great discoveries today!"
  if (!is.null(msg)) subtitle <- msg
  splash.title <- shiny::div(
    shiny::br(), br(), br(), br(),
    shiny::div(shiny::HTML(title), style = "font-size:70px;font-weight:700;line-height:1em;width:130%;"),
    shiny::br(),
    shiny::div(shiny::HTML(subtitle), style = "font-size:30px;line-height:1em;margin-top:0.6em;width:130%;"),
    shiny::br(), br(), br()
  )
  body <- shiny::tagList(
    shiny::div(id = "splash-title", splash.title)
  )
  m <- splashScreen("", body,
    ns = ns, easyClose = TRUE, fade = TRUE,
    buttons = FALSE, footer = FALSE
  )
  if (duration > 0) {
    cat("closing hello in", round(duration / 1000, 1), "seconds...\n")
    shinyjs::delay(duration, shiny::removeModal())
  }
  return(m)
}

splashLoginModal <- function(ns = NULL,
                             with.email = TRUE,
                             with.password = TRUE,
                             with.username = FALSE,
                             with.firebase = FALSE,
                             with.firebase_emailonly = FALSE,
                             button.text = "Login",
                             cancel.text = "cancel",
                             add.cancel = FALSE,
                             title = "Log in",
                             subtitle = "") {
  if (is.null(ns)) {
    ns <- function(e) {
      return(e)
    }
  }

  slogan <- list()
  slogan[[1]] <- c("Big Omics Data", "Isn't big anymore with BigOmics Playground")
  slogan[[2]] <- c("Great Discoveries", "Start on BigOmics Playground")
  slogan[[3]] <- c("Fasten Your Seat Belts!", "Hi-speed analytics")
  slogan[[4]] <- c("Do-it-yourself Omics Analytics", "Yes you can!")
  slogan[[5]] <- c("Twenty-Four Seven", "Your Playground doesn't go on coffee breaks")
  slogan[[6]] <- c("Analyze with confidence", "Be a data rockstar, a Freddie Mercury of omics!")
  slogan[[7]] <- c("Play-Explore-Discover", "Get deeper insights with BigOmics Playground")
  slogan[[8]] <- c("Skip the Queue", "Take the fast lane. Self-service analytics.")
  slogan[[9]] <- c("Look Ma! No help!", "I did it without a bioinformatician")
  slogan[[10]] <- c("Easy-peasy insight!", "Get insight from your data the easy way")
  slogan[[11]] <- c("Zoom-zoom-insight!", "Get faster insight from your data")
  slogan[[12]] <- c("Click-click-eureka!", "Owe yourself that <i>eureka!</i> moment")
  slogan[[13]] <- c("I Love Omics Data!", "Unleash your inner nerd with BigOmics Playground")
  slogan[[14]] <- c("More Omics Data", "Is all I want for my birthday")
  slogan[[15]] <- c("Keep Exploring", "Never stop discovering with BigOmics Playground")
  slogan[[16]] <- c("Real Bioinformaticians", "Do it with BigOmics Playground")
  slogan[[17]] <- c("Real Biologists", "Do it with BigOmics Playground")
  slogan[[18]] <- c("Ich bin doch nicht bl\u00F6d!", "Of course I use BigOmics Playground")
  slogan[[19]] <- c("Non sono mica scemo!", "Of course I use BigOmics Playground")
  slogan[[20]] <- c("The Unexplored Plan", "When you get into exploring, you realize that we live on a relatively unexplored plan. &ndash; E. O. Wilson")
  slogan[[21]] <- c("Explore More", "The more you explore, the more you learn and grow")
  slogan[[22]] <- c("Discover New Oceans", "Man cannot discover new oceans unless he has the courage to lose sight of the shore. &ndash; Andre Gide")
  slogan[[23]] <- c("Love Adventurous Life", "Be passionately curious about exploring new adventures. &ndash; Lailah Gifty Akita")
  slogan[[24]] <- c("Succes is Exploration", "Find the unknown. Learning is searching. Anything else is just waiting. &ndash; Dale Daute")
  slogan[[25]] <- c("Look Ma! No help!", "I did it without a bioinformagician")
  slogan[[26]] <- c("May the Force of Omics be with you", "Train hard youngling, one day a master you become")

  slogan <- sample(slogan, 1)[[1]]
  slogan.len <- nchar(paste(slogan, collapse = " "))
  if (slogan.len < 80) slogan[1] <- paste0("<br>", slogan[1])
  splash.title <- shiny::div(
    id = "splash-title",
    class = "text-white",
    shiny::div(shiny::HTML(slogan[1]), style = "font-size:3rem;font-weight:700;line-height:1em;width:130%;"),
    shiny::div(shiny::HTML(slogan[2]), style = "font-size:1.6rem;line-height:1.1em;margin-top:0.5em;width:130%;")
  )

  div.password <- div()
  div.email <- div()
  div.username <- div()
  div.firebase <- div()
  div.title <- div()
  div.subtitle <- div()

  if (with.email) {
    div.email <- div(
      id = "splash-email",
      textInput(ns("login_email"), NULL, placeholder = "your email")
    )
  }
  if (with.username) {
    div.username <- div(
      id = "splash-username",
      textInput(ns("login_username"), NULL, placeholder = "your username")
    )
  }
  if (with.password) {
    div.password <- div(
      id = "splash-password",
      passwordInput(ns("login_password"), NULL, placeholder = "your password")
    )
  }
  if (with.firebase) {
    div.firebase <- div(
      class = "card",
      div(
        class = "card-body",
        h1(
          title,
          class = "card-title pb-2"
        ),
        div(
          subtitle
        ),
        textInput(
          ns("emailInput"),
          "",
          placeholder = "Your email",
          width = "100%"
        ),
        actionButton(
          ns("emailSubmit"),
          "Send link",
          class = "btn-warning"
        ),
        p(
          id = "emailFeedbackShow"
        ),
        hr(style = "color:#888;opacity:1;margin-top:30px;"),
        h5(
          "or",
          class = "text-center pb-4 pt-1",
          style = "margin-top:-35px;background:white;width:50px;margin-left:auto;margin-right:auto;"
        ),
        div(
          class = "social-button google-button",
          actionLink(
            ns("launchGoogle"),
            HTML("&nbsp; Sign in with Google"),
            icon = shiny::icon("google", style = "font-size:18px;")
          )
        ),
        div(
          class = "social-button facebook-button",
          actionLink(
            ns("launchFacebook"),
            HTML("&nbsp; Sign in with Facebook"),
            icon = shiny::icon("facebook", style = "font-size:18px;")
          )
        ),
        ## div(
        ##   class = "social-button apple-button",
        ##   actionLink(
        ##     ns("launchApple"),
        ##     HTML("&nbsp; Sign in with Apple"),
        ##     icon = shiny::icon("apple", style="font-size:18px;")
        ##   )
        ## ),
        div(
          class = "social-button twitter-button",
          actionLink(
            ns("launchTwitter"),
            HTML("&nbsp; Sign in with Twitter"),
            icon = shiny::icon("twitter", style = "font-size:18px;")
          )
        )
      )
    )
  }
  if (with.firebase_emailonly) {
    div.firebase <- div(
      class = "card",
      div(
        class = "card-body",
        h1(
          title,
          class = "card-title pb-2"
        ),
        div(
          subtitle
        ),
        textInput(
          ns("emailInput"),
          "",
          placeholder = "Your email",
          width = "100%"
        ),
        actionButton(
          ns("emailSubmit"),
          "Send link",
          class = "btn-warning"
        ),
        p(
          id = "emailFeedbackShow"
        )
      )
    )
  }

  if (!is.null(title) && title != "") {
    div.title <- div(
      id = "splash-login-title",
      class = "pb-3",
      h1(title, style = "color:black;line-height:1em;")
    )
  }

  if (!is.null(subtitle) && subtitle != "") {
    div.subtitle <- div(
      id = "splash-login-subtitle",
      class = "pt-0 pb-2",
      h6(subtitle, style = "color:black;font-weight:400;")
    )
  }

  if (add.cancel) {
    div.button <- div(
      id = "splash-buttons",
      class = "pt-2",
      shiny::fillRow(
        flex = c(1, NA, NA, 1),
        br(),
        actionButton(ns("cancel_btn"), cancel.text,
          class = "btn-light btn-xl",
          style = "margin: 4px;"
        ),
        actionButton(ns("login_btn"), button.text,
          class = "btn-warning btn-xl",
          style = "margin: 4px;"
        ),
        br()
      )
    )
  } else {
    div.button <- div(
      id = "splash-buttons",
      class = "pt-2",
      actionButton(ns("login_btn"), button.text, class = "btn-warning btn-xl")
    )
  }

  ## splash.panel=div();ns=function(x)
  splash.content <- NULL
  if (with.firebase || with.firebase_emailonly) {
    splash.content <- div.firebase
  } else {
    splash.content <- shiny::wellPanel(
      ## style = "padding: 40px 20px; background-color: #ffffff22;",
      style = "padding: 35px 25px; background-color:white; color:black;",
      id = "splash-login",
      div.title,
      div.subtitle,
      div.username,
      div.email,
      div.password,
      div.button
    )
  }

  body <- div(
    id = "splash-content",
    splash.content
  )

  m <- splashScreen(title = splash.title, body = body, ns = ns)
  return(m)
}

splashscreen.buttons <- function() {
  tagList(
    shiny::tags$a(
      shiny::img(
        id = "splash-logo2",
        src = "static/bigomics-logo.png"
      ),
      href = "https://www.bigomics.ch",
      target = "_blank"
    ),
    div(
      class = "btn-group",
      role = "group",
      div(
        class = "btn-group",
        role = "group",
        tags$button(
          "Get support",
          id = "splash-toggle-support",
          type = "button",
          class = "btn btn-outline-primary dropdown-toggle",
          `data-bs-toggle` = "dropdown",
          `aria-expanded` = "false"
        ),
        tags$ul(
          class = "dropdown-menu",
          `aria-labelledby` = "splash-toggle-support",
          tags$li(
            shiny::tags$a(
              "Watch tutorials",
              class = "dropdown-item",
              href = "https://www.youtube.com/channel/UChGASaLbr63pxmDOeXTQu_A",
              target = "_blank"
            )
          ),
          tags$li(
            shiny::tags$a(
              "Read documentation",
              class = "dropdown-item",
              href = "https://omicsplayground.readthedocs.io",
              target = "_blank"
            ),
          ),
          tags$li(
            shiny::tags$a(
              "User forum",
              class = "dropdown-item",
              href = "https://groups.google.com/d/forum/omicsplayground",
              target = "_blank"
            )
          )
        )
      ),
      div(
        class = "btn-group",
        role = "group",
        tags$button(
          "I'm a developer",
          id = "splash-toggle-dev",
          type = "button",
          class = "btn btn-outline-primary dropdown-toggle",
          `data-bs-toggle` = "dropdown",
          `aria-expanded` = "false"
        ),
        tags$ul(
          class = "dropdown-menu",
          `aria-labelledby` = "splash-toggle-dev",
          tags$li(
            shiny::tags$a(
              "Get the source",
              class = "dropdown-item",
              href = "https://github.com/bigomics/omicsplayground",
              target = "_blank"
            )
          ),
          tags$li(
            shiny::tags$a(
              "Docker image",
              class = "dropdown-item",
              href = "https://hub.docker.com/r/bigomics/omicsplayground",
              target = "_blank"
            )
          ),
          tags$li(
            shiny::tags$a(
              "Buy us a coffee!",
              class = "dropdown-item",
              href = "https://www.buymeacoffee.com/bigomics",
              target = "_blank"
            )
          )
        )
      )
    )
  )
}

splashScreen <- function(title, body, ns = NULL, easyClose = FALSE, fade = FALSE,
                         buttons = TRUE, footer = TRUE) {
  if (is.null(ns)) {
    ns <- function(e) {
      return(e)
    }
  }

  div.buttons <- shiny::modalButton("Dismiss")
  if (buttons) {
    div.buttons <- splashscreen.buttons()
  }
  if (!footer) {
    div.buttons <- NULL
  }

  ## return modalDialog
  m <- modalDialog2(
    id = "splash-fullscreen",
    class = "bg-primary",
    header = div.buttons,
    shiny::div(
      class = "row",
      shiny::div(
        class = "col-md-4 offset-md-2",
        title,
        br(),
        br(),
        shiny::img(src = "static/mascotte-sc.png", class = "img-fluid", id = "splash-image"),
      ),
      shiny::div(
        class = "col-md-3 offset-md-2",
        shiny::div(
          id = "splash-panel",
          body,
          div(textOutput(ns("login_warning")), style = "color:white;font-size:1.2em;"),
        ),
      )
    ),
    footer = NULL,
    size = "fullscreen",
    easyClose = easyClose,
    fade = fade
  ) ## end of modalDialog

  return(m)
}
