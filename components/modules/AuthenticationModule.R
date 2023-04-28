##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

AuthenticationUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
}

NoAuthenticationModule <- function( id,
                                   show_modal = TRUE,
                                   username = "",
                                   email = "") {
  shiny::moduleServer( id, function(input, output, session) {
    
    message("[NoAuthenticationModule] >>>> using no authentication <<<<")
    ns <- session$ns
    USER <- shiny::reactiveValues(
        logged=FALSE,
        name="",
        email="",
        level="",
        limit=""
    )

    resetUSER <- function() {
        USER$logged <- FALSE
        USER$name <- ""
        USER$email <- ""
        USER$level <- ""
        USER$limit <- ""
        if(show_modal) {
            m <- splashLoginModal(
              ns=ns,
              with.email=FALSE,
              with.password=FALSE,
              alt=h5("Ready to explore your data?",style="color:white;"),
              login.text="Sure I am!"
            )
            shiny::showModal(m)
        } else {
            USER$logged <- TRUE
        }
        USER$name   <- username
        USER$email  <- email
    }

    output$showLogin <- shiny::renderUI({
        resetUSER()
    })

    output$login_warning = shiny::renderText("")

    shiny::observeEvent( input$login_btn, {
        shiny::removeModal()
        ##shiny::showModal(splashHelloModal(name=USER$name,ns=ns))
        USER$logged <- TRUE
        ##Sys.sleep(3);cat("wait 3 seconds to close...\n");removeModal()
        session$sendCustomMessage("set-user", list(user = USER$email))
    })

    observeEvent( input$userLogout, {
        resetUSER()
    })

    rt <- list(
        name   = shiny::reactive(USER$name),
        email  = shiny::reactive(USER$email),
        level  = shiny::reactive(USER$level),
        logged = shiny::reactive(USER$logged),
        limit  = shiny::reactive(USER$limit),
        stripe_id  = shiny::reactive(''),
        href  = shiny::reactive('')
    )
    return(rt)
  } ## end-of-server
 )
}


##================================================================================
## FirebaseAuthenticationModule
##================================================================================

upgrade.dialog <- function(ns, current.plan) {

    btn_basic   <- "Go Basic!"
    btn_starter <- "Get Starter!"
    btn_premium <- "Get Premium!"
    if(current.plan=='free') btn_basic <- "Current Plan"
    if(current.plan=='starter') btn_starter <- "Current Plan"
    if(current.plan=='premium') btn_premium <- "Current Plan"

    modalDialog(
        title = h3("Find the right OmicsPlayground plan for you"),
        size = "m",
        ##p("Do you want to remove the 60 minutes time limit? Do you want to be able to save more datasets?"),
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
            shiny::actionButton(ns("get_basic"),btn_basic),
            br()
          ),
          div(
            class = "col-md-4",
            style = "background:#E8F8FF;",
            h4(HTML("<b>Starter</b>")),
            p("Great to start"),
            ##            h3("CHF49 / month", id = "starter-pricing"),
            h3("Soon!"),
            tags$ul(
              class = "list-unstyled",
              tags$li("Host up to 10 datasets"),
              tags$li("3 hours time limit"),
              tags$li("Up to 100 samples / dataset"),
              tags$li("Up to 10 comparisons")
            ),
            shiny::actionButton(ns("get_starter"),btn_starter),
            ##shiny::actionButton(ns("get_starter"),"Get Starter!", onClick='upgrade_plan()'),
            br()
          ),
          div(
              class = "col-md-4",
              style = "background:#E2F4FF;",
              HTML("<h4><b>Premium</b></h4>"),
              p("For power users or small groups"),
              ##h3("CHF490 / month", id = "premium-pricing"),
              h3("Soon!"),
              tags$ul(
                       class = "list-unstyled",
                       tags$li("Host up to 100 datasets"),
                       tags$li("8 hours time limit"),
                       tags$li("Up to 2000 samples / dataset"),
                       tags$li("Up to 100 comparisons")
                   ),
              ##shiny::actionButton(ns("get_premium"),"Get Premium!", onClick='get_premium()'),
              shiny::actionButton(ns("get_premium"),btn_premium),
              br()
          )
        ),  ## content div
        div(
            style = "margin-top:3rem;text-align:center;",
            HTML("Looking for OmicsPlayground for <b>Enterprise</b>? <a href='mailto:info@bigomics.com'>Contact sales for info and pricing</a>.")
        ),
        footer = tagList(
            fillRow(
                flex=c(NA,0.03,NA,1,NA,NA),
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
                shiny::actionButton(ns("manage"),"Manage Subscription"),
                modalButton("Dismiss")
            )
        )
    )  ## modalDialog
}

FirebaseAuthenticationModule <- function(id, firebase.rds="firebase.rds") {
  shiny::moduleServer( id, function(input, output, session) {
    message("[AuthenticationModule] >>>> using FireBase authentication <<<<")

    dbg("[AuthenticationModule] getwd = ",getwd())
    dbg("[AuthenticationModule] file.exists('firebase.rds') = ",file.exists(firebase.rds))

    if(file.exists(firebase.rds)) {
      firebase_config <- firebase:::read_config(firebase.rds)
    } else {
      stop("[FATAL ERROR] no firebase.rds file found. please create.")
    }
    Sys.setenv(OMICS_GOOGLE_PROJECT = firebase_config$projectId)

    ns <- session$ns
    USER <- shiny::reactiveValues(
        logged = FALSE,
        name = "",
        password = "",
        email = "",
        level = "",
        limit = "",
        token = NULL,
        uid = NULL,
        stripe_id = NULL,
        href = NULL
    )

    firebase <- firebase::FirebaseSocial$
        new(persistence = "local")

    firebase2 <- firebase::FirebaseEmailLink$
        new(persistence = "local")

    observeEvent(input$launchGoogle, {
        firebase$launch_google(flow = "popup")
    })

    resetUSER <- function() {

        message("[FirebaseAuthenticationModule] resetting USER... ")

        USER$logged <- FALSE
        USER$name <- ""
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
            with.firebase = FALSE,
            with.firebase_emailonly = TRUE,
            alt = h5("Ready to explore your data?",style="color:white;"),
            login.text = "Sure I am!"
        )

        ## login modal
        shiny::showModal(m)
    }

    first_time = TRUE
    observeEvent(USER$logged, {

        ## no need to show the modal if the user is logged this is due
        ## to persistence. But if it is the first time of the session
        ## we force reset/logout to delete sleeping logins.
        if(USER$logged && !first_time) {
            dbg("[FirebaseAuthenticationModule] USER is already logged in! no modal")
            return()
        }

        first_time <<- FALSE
        message("[FirebaseAuthenticationModule] USER not logged in!")
        resetUSER()

    })

    observeEvent( input$userLogout, {
        message("[FirebaseAuthenticationModule] userLogout triggered!")
        resetUSER()
    })

    observeEvent( input$emailSubmit, {
        if(input$emailInput == ""){
            session$sendCustomMessage(
                "email-feedback",
                list(
                    type = "error",
                    msg = "Missing email"
                )
            )
            return()
        }

        ## >>> We could check here for email validaty and intercept the
        ## login process for not authorized people with wrong domain
        ## or against a subscription list.
        authorized <- grepl("@bigomics.ch$",input$emailInput)
        if(FALSE && !authorized) {
          shinyalert::shinyalert("We're sorry...","You are not authorized to log in. Please contact your systems administrator.")
          resetUSER()
          return()
        }
        
        is_personal_email <- grepl("gmail|ymail|outlook|yahoo|mail.com$|icloud",input$emailInput)
        if(TRUE && is_personal_email) {
          shinyalert::shinyalert("We're sorry...","No personal email allowed. Please log in with your business, academic or institutional email.")
          shiny::updateTextInput(session, "emailInput", value="")
          return()
        }

        ## >>> OK let's send auth request
        session$sendCustomMessage(
            "email-feedback",
            list(
                type = "success",
                msg = "Email sent, check your inbox."
            )
        )
        firebase2$send_email(input$emailInput)
    })

    observeEvent( firebase$get_signed_in(), {

        response <- firebase$get_signed_in()

        if(!response$success) {
            warning("[FirebaseAuthenticationModule] sign in NOT succesful")
            resetUSER()
            return()
        }

        on.exit({
            dbg("[FirebaseAuthenticationModule] get_signed_in() on.exit")
            if(USER$logged) removeModal()
        })

##        for(i in 1:length(response$response)) {
##            dbg("[FirebaseAuthentication] ",names(response$response)[i],"=",response$response[[i]])
##        }

        ## >>> We could check here for email validaty and intercept the
        ## login process for not authorized people with wrong domain
        ## or against a subscription list.
        authorized <- grepl("@bigomics.ch$",response$response$email)
        if(FALSE && !authorized) {
          shinyalert::shinyalert("We're sorry","You are not authorized to log in. Please contact your systems administrator.")
          USER$logged <- FALSE
          return(NULL)
        }
        
        USER$logged <- TRUE
        USER$uid <- as.character(response$response$uid)
        USER$name  <- response$response$displayName
        USER$email <- response$response$email

        if(!is.null(USER$name))  USER$name  <- as.character(USER$name)
        if(!is.null(USER$email)) USER$email <- as.character(USER$email)
        if(is.null(USER$name))   USER$name  <- ""
        if(is.null(USER$email))  USER$email <- ""

        dbg("[FirebaseAuthenticationModule@firebase$get_signed_in] user.name=='' = ",USER$name=='' )
        dbg("[FirebaseAuthenticationModule@firebase$get_signed_in] user.email=='' = ",USER$email=='' )

        session$sendCustomMessage("get-permissions", list(ns = ns(NULL)))
    })

    observeEvent(input$stripeId, {
        USER$stripe_id <- input$stripeId$id
        USER$href      <- input$stripeId$href
    })

    observeEvent(input$permissions, {
        perm <- input$permissions

        USER$level <- "free"
        if(perm$success)
            USER$level <- "premium"

        session$sendCustomMessage(
            "set-user",
            list(
                user = USER$email,
                level = USER$level,
                pricing = Sys.getenv("OMICS_STRIPE_PREMIUM_PRICE")
            )
        )
    })

    observeEvent( input$firebaseUpgrade, {
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
        session$sendCustomMessage('manage-sub', content$url)
    })

    rt <- list(
        name   = shiny::reactive(USER$name),
        email  = shiny::reactive(USER$email),
        level  = shiny::reactive(USER$level),
        logged = shiny::reactive(USER$logged),
        limit  = shiny::reactive(USER$limit),
        stripe_id  = shiny::reactive(USER$stripe_id),
        href  = shiny::reactive(USER$href)
    )
    return(rt)
  }
)}

##================================================================================
## EmailLinkeAuthenticationModule (using Firebase, no stripe)
##================================================================================

EmailLinkAuthenticationModule <- function(id, pgx_dir, firebase.rds="firebase.rds") {
  shiny::moduleServer( id, function(input, output, session) {
    message("[AuthenticationModule] >>>> using email link (using Firebase) authentication <<<<")

    dbg("[AuthenticationModule] getwd = ",getwd())
    dbg("[AuthenticationModule] file.exists('firebase.rds') = ",file.exists(firebase.rds))

    if(file.exists(firebase.rds)) {
      firebase_config <- firebase:::read_config(firebase.rds)
    } else {
      stop("[FATAL ERROR] no firebase.rds file found. please create.")
    }
    Sys.setenv(OMICS_GOOGLE_PROJECT = firebase_config$projectId)

    ns <- session$ns
    USER <- shiny::reactiveValues(
        logged = FALSE,
        name = "",
        password = "",
        email = "",
        level = "",
        limit = "",
        token = NULL,
        uid = NULL,
        stripe_id = NULL,
        href = NULL
    )

    firebase <- firebase::FirebaseSocial$
        new(persistence = "local")

    firebase2 <- firebase::FirebaseEmailLink$
        new(persistence = "local")

    observeEvent(input$launchGoogle, {
        firebase$launch_google(flow = "popup")
    })

    resetUSER <- function() {

        message("[FirebaseAuthenticationModule] resetting USER... ")

        USER$logged <- FALSE
        USER$name <- ""
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
            with.firebase = FALSE,
            with.firebase_emailonly = TRUE,
            alt = h5("Ready to explore your data?",style="color:white;"),
            login.text = "Sure I am!"
        )

        ## login modal
        shiny::showModal(m)
    }

    first_time = TRUE
    observeEvent(USER$logged, {

        ## no need to show the modal if the user is logged this is due
        ## to persistence. But if it is the first time of the session
        ## we force reset/logout to delete sleeping logins.
        if(USER$logged && !first_time) {
            dbg("[FirebaseAuthenticationModule] USER is already logged in! no modal")
            return()
        }

        first_time <<- FALSE
        message("[FirebaseAuthenticationModule] USER not logged in!")
        resetUSER()

    })

    observeEvent( input$userLogout, {
        message("[FirebaseAuthenticationModule] userLogout triggered!")
        resetUSER()
    })

    observeEvent( input$emailSubmit, {
        if(input$emailInput == ""){
            session$sendCustomMessage(
                "email-feedback",
                list(
                    type = "error",
                    msg = "Missing email"
                )
            )
            return()
        }

        ## >>> We could check here for email validaty and intercept the
        ## login process for not authorized people with wrong domain
        ## or against a subscription list.
        authorized <- grepl("@bigomics.ch$",input$emailInput)
        if(FALSE && !authorized) {
          shinyalert::shinyalert("We're sorry...","You are not authorized to log in. Please contact your systems administrator.")
          resetUSER()          
          return()
        }

        ## if it is a new user we ask for business email, old users can go
        is_personal_email <- grepl("gmail|ymail|outlook|yahoo|mail.com$|icloud",input$emailInput)
        existing_user_dirs <- basename(list.dirs(pgx_dir))
        new_user <- !(input$emailInput %in% existing_user_dirs)
        if(TRUE && is_personal_email && new_user) {
          shinyalert::shinyalert("We're sorry...","No personal email allowed. Please provide your business, academic or institutional email.")
          shiny::updateTextInput(session, "emailInput", value="")
          return()
        }
        
        ## >>> OK let's send auth request
        session$sendCustomMessage(
            "email-feedback",
            list(
                type = "success",
                msg = "Email sent, check your inbox."
            )
        )
        firebase2$send_email(input$emailInput)
    })

    observeEvent( firebase$get_signed_in(), {

        response <- firebase$get_signed_in()

        if(!response$success) {
            warning("[FirebaseAuthenticationModule] sign in NOT succesful")
            resetUSER()
            return()
        }

        on.exit({
            dbg("[FirebaseAuthenticationModule] get_signed_in() on.exit")
            if(USER$logged) removeModal()
        })

##      for(i in 1:length(response$response)) {
##          dbg("[FirebaseAuthentication] ",names(response$response)[i],"=",response$response[[i]])
##      }

        ## >>> We could check here for email validaty and intercept the
        ## login process for not authorized people with wrong domain
        ## or against a subscription list.
        authorized <- grepl("@bigomics.ch$",response$response$email)
        if(FALSE && !authorized) {
          shinyalert::shinyalert("We're sorry","You are not authorized to log in. Please contact your systems administrator.")
          USER$logged <- FALSE
          resetUSER()          
          return(NULL)
        }
        
        USER$logged <- TRUE
        USER$uid <- as.character(response$response$uid)
        USER$name  <- response$response$displayName
        USER$email <- response$response$email

        if(!is.null(USER$name))  USER$name  <- as.character(USER$name)
        if(!is.null(USER$email)) USER$email <- as.character(USER$email)
        if(is.null(USER$name))   USER$name  <- ""
        if(is.null(USER$email))  USER$email <- ""

        dbg("[FirebaseAuthenticationModule@firebase$get_signed_in] user.name=='' = ",USER$name=='' )
        dbg("[FirebaseAuthenticationModule@firebase$get_signed_in] user.email=='' = ",USER$email=='' )

        session$sendCustomMessage("get-permissions", list(ns = ns(NULL)))
    })

    rt <- list(
        name   = shiny::reactive(USER$name),
        email  = shiny::reactive(USER$email),
        level  = shiny::reactive(USER$level),
        logged = shiny::reactive(USER$logged),
        limit  = shiny::reactive(USER$limit),
        href   = shiny::reactive(USER$href)
    )
    return(rt)
  }
)}

##================================================================================
## PasswordAuthenticationModule (ask login.name + password)
##================================================================================

credentials.file='CREDENTIALS'
PasswordAuthenticationModule <- function(id, credentials.file) {
  shiny::moduleServer( id, function(input, output, session) {  
    message("[AuthenticationModule] >>>> using local Email+Password authentication <<<<")

    ns <- session$ns
    USER <- shiny::reactiveValues(
                       logged=FALSE,
                       ## username=NA,
                       email=NA,
                       password=NA,
                       level=NA,
                       limit=NA)

    resetUSER <- function() {

        USER$logged <- FALSE
        ## USER$username <- NA
        USER$email <- NA
        USER$password <- NA
        USER$level <- ""
        USER$limit <- ""

        m <- splashLoginModal(
            ns=ns,
            with.email=TRUE,
            with.username=FALSE,
            with.password=TRUE,
            login.text="Let me in!"
            )
        shiny::showModal(m)

    }

    CREDENTIALS <- read.csv(credentials.file,colClasses="character")
    head(CREDENTIALS)

    output$showLogin <- shiny::renderUI({
        m <- splashLoginModal(
            ns=ns,
            with.email=TRUE,
            with.username=FALSE,
            with.password=TRUE,
            login.text="Let me in!"
        )
        shiny::showModal(m)
    })

    output$login_warning = shiny::renderText("")

    shiny::observeEvent( input$login_btn, {

        login.OK   = FALSE
        valid.date = FALSE
        valid.user = FALSE

        login_email    <- input$login_email
        login_password <- input$login_password

        if( is.null(login_email) || is.null(login_password)) return(NULL)
        if( login_email=="" || login_password=="") return(NULL)
        sel <- tail(which( CREDENTIALS$email == login_email),1)
        dbg("[AuthenticationModule:input$login_btn] CREDENTIALS$email = ",CREDENTIALS$email)
        dbg("[AuthenticationModule:input$login_btn] sel = ",sel)

        valid.user <- isTRUE(CREDENTIALS$email[sel] == login_email) && length(sel)>0
        valid.pw   <- isTRUE(CREDENTIALS[sel,"password"] == input$login_password)
        valid.date <- isTRUE(Sys.Date() < as.Date(CREDENTIALS[sel,"expiry"]) )
        login.OK = (valid.user && valid.pw && valid.date)

        message("--------- password login ---------")
        ##message("input.username = ",input$login_username)
        message("input.email    = ",input$login_email)
        message("input.password = ",input$login_password)
        message("user.password  = ",CREDENTIALS[sel,"password"])
        message("user.expiry    = ",CREDENTIALS[sel,"expiry"])
        message("user.email     = ",CREDENTIALS[sel,"email"])
        message("user.limit     = ",CREDENTIALS[sel,"limit"])
        message("valid.user     = ",valid.user)
        message("valid.date     = ",valid.date)
        message("valid.pw       = ",valid.pw)
        message("----------------------------------")

        if (login.OK) {

            message("[PasswordAuthenticationModule::login] PASSED : login OK! ")

            output$login_warning = shiny::renderText("")
            shiny::removeModal()
            ##USER$name   <- input$login_username
            ##USER$email <- CREDENTIALS[sel,"email"]
            cred <- CREDENTIALS[sel,]
            ##USER$username  <- cred$username
            USER$email     <- cred$email
            USER$level     <- cred$level
            USER$limit     <- cred$limit

            ## Here you can perform some user-specific functions, site
            ## news, or 2nd hello modal...
            ##shiny::showModal(splashHelloModal(USER$name,ns=ns))
            ##removeModal()
            USER$logged <- TRUE
            session$sendCustomMessage("set-user", list(user = USER$email))

        } else {
            message("[PasswordAuthenticationModule::login] REFUSED : invalid login! ")
            if(!valid.date) {
                output$login_warning = shiny::renderText("Registration expired")
            }
            if(!valid.pw) {
                output$login_warning = shiny::renderText("Invalid password")
            }
            if(!valid.user) {
                output$login_warning = shiny::renderText("Invalid user")
            }
            shinyjs::delay(2000, {output$login_warning <- shiny::renderText("")})
            USER$logged <- FALSE
        }
    })

    observeEvent( input$userLogout, {
        resetUSER()
    })

    ## module reactive return value
    rt <- list(
        name   = shiny::reactive(""),
        email  = shiny::reactive(USER$email),
        level  = shiny::reactive(USER$level),
        logged = shiny::reactive(USER$logged),
        limit  = shiny::reactive(USER$limit)
    )
    return(rt)
  }
)}

##================================================================================
## HELPER FUNCTIONS
##================================================================================

splashHelloModal <- function(name, msg=NULL, ns=NULL, duration=3500)
{
    if(is.null(ns)) ns <- function(e) return(e)
    message("[AuthenticationModule::splashHelloModel]")

    all.hello = c("Hello","Salut","Hola","Pivet","Ni hao","Ciao","Hi","Hoi","Hej",
                  "Yassou","Selam","Hey","Hei","Grutzi","Bonjour",
                  "Namaste","Salam","Selamat","Shalom","Goeiedag","Yaxshimusiz")
    title = paste(paste0(sample(all.hello,3),"!"), collapse=" ")
    if(!is.null(name) && !is.na(name) && !name %in% c("NA","")) {
        first.name <- strsplit(as.character(name),split=" ")[[1]][1]
        first.name <- paste0(toupper(substring(first.name,1,1)),
                             substring(first.name,2,999))
        ##title = paste0(sample(all.hello,1)," ",first.name,"!")
        title = paste(paste0(sample(all.hello,1)," ",first.name,"!"),collapse=" ")
    }
    subtitle = "Have a good day!"
    subtitle = "We wish you many great discoveries today!"
    if(!is.null(msg)) subtitle <- msg
    splash.title <- shiny::div(
        shiny::br(),br(),br(),br(),
        shiny::div(shiny::HTML(title),style="font-size:70px;font-weight:700;line-height:1em;width:130%;"),
        shiny::br(),
        shiny::div(shiny::HTML(subtitle),style="font-size:30px;line-height:1em;margin-top:0.6em;width:130%;"),
        shiny::br(),br(),br()
    )
    body <- shiny::tagList(
        shiny::div(id="splash-title",splash.title)
    )
    m <- splashScreen("", body, ns=ns, easyClose=TRUE, fade=TRUE,
        buttons=FALSE, footer=FALSE)
    if(duration>0) {
        cat("closing hello in",round(duration/1000,1),"seconds...\n")
        shinyjs::delay(duration, shiny::removeModal())
    }
    return(m)
}

splashLoginModal <- function(ns=NULL, with.email=TRUE, with.password=TRUE,
                             with.username=FALSE, 
                             with.firebase=FALSE,
                             with.firebase_emailonly=FALSE,
                             login.text="Login", alt=NULL)
{
    if(is.null(ns)) ns <- function(e) return(e)
    message("[AuthenticationModule::splashLoginModal]")

    titles <- list()
    titles[[1]] = c("Big Omics Data","Isn't big anymore with BigOmics Playground")
    titles[[2]] = c("Great Discoveries","Start on BigOmics Playground")
    titles[[3]] = c("Fasten Your Seat Belts!","Hi-speed analytics")
    titles[[4]] = c("Do-it-yourself Omics Analytics","Yes you can!")
    titles[[5]] = c("Twenty-Four Seven","Your Playground doesn't go on coffee breaks")
    titles[[6]] = c("Analyze with confidence","Be a data rockstar, a Freddie Mercury of omics!")
    titles[[7]] = c("Play-Explore-Discover","Get deeper insights with BigOmics Playground")
    titles[[8]] = c("Skip the Queue","Take the fast lane. Self-service analytics.")
    titles[[9]] = c("Look Ma! No help!","I did it without a bioinformatician")
    titles[[10]] = c("Easy-peasy insight!","Get insight from your data the easy way")
    titles[[11]] = c("Zoom-zoom-insight!","Get faster insight from your data")
    titles[[12]] = c("Click-click-eureka!","Owe yourself that <i>eureka!</i> moment")
    titles[[13]] = c("I Love Omics Data!","Unleash your inner nerd with BigOmics Playground")
    titles[[14]] = c("More Omics Data","Is all I want for my birthday")
    titles[[15]] = c("Keep Exploring","Never stop discovering with BigOmics Playground")
    titles[[16]] = c("Real Bioinformaticians","Do it with BigOmics Playground")
    titles[[17]] = c("Real Biologists","Do it with BigOmics Playground")
    titles[[18]] = c("Ich bin doch nicht bl\u00F6d!","Of course I use BigOmics Playground")
    titles[[19]] = c("Non sono mica scemo!","Of course I use BigOmics Playground")
    titles[[20]] = c("The Unexplored Plan","When you get into exploring, you realize that we live on a relatively unexplored plan. &ndash; E. O. Wilson")
    titles[[21]] = c("Explore More","The more you explore, the more you learn and grow")
    titles[[22]] = c("Discover New Oceans","Man cannot discover new oceans unless he has the courage to lose sight of the shore. &ndash; Andre Gide")
    titles[[23]] = c("Love Adventurous Life","Be passionately curious about exploring new adventures. &ndash; Lailah Gifty Akita")
    titles[[24]] = c("Succes is Exploration","Find the unknown. Learning is searching. Anything else is just waiting. &ndash; Dale Daute")
    titles[[25]] = c("Look Ma! No help!","I did it without a bioinformagician")
    titles[[26]] = c("May the Force of Omics be with you","Train hard youngling, one day a master you become")
    title <- titles[[length(titles)]]
    title <- sample(titles,1)[[1]]
    title.len <- nchar(paste(title,collapse=' '))
    if(title.len < 80) title[1] <- paste0("<br>",title[1])
    splash.title <- shiny::div(
        class = "text-white",
        shiny::div(shiny::HTML(title[1]),style="font-size:3.2rem;font-weight:700;line-height:1em;width:130%;"),
        shiny::div(shiny::HTML(title[2]),style="font-size:1.8rem;line-height:1.1em;margin-top:0.5em;width:130%;")
    )

    div.password <- div()
    div.email <- div()
    div.username <- div()
    div.firebase <- div()

    if(with.email) {
        div.email <- div(
            id="splash-email",
            textInput(ns("login_email"),NULL,placeholder="login with your email")
        )
    }
    if(with.username) {
        div.email <- div(
            id="splash-username",
            textInput(ns("login_username"),NULL,placeholder="login with your username")
        )
    }
    if(with.password) {
        div.password <- div(
            id="splash-password",
            passwordInput(ns("login_password"),NULL,placeholder="enter your password")
        )
    }
    if(with.firebase) {
        div.firebase <- div(
            class = "card",
            div(
                class = "card-body",
                h1(
                    "Sign in",
                    class = "card-title pb-2"
                ),
                div(
                    "Enter your email and we'll send you a link."
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
                hr( style="color:#888;opacity:1;margin-top:30px;"),
                h5(
                    "or",
                    class = "text-center pb-4 pt-1",
                    style = "margin-top:-35px;background:white;width:50px;margin-left:auto;margin-right:auto;"
                ),
                div(
                    class = "social-button google-button",
                    actionLink(
                        ns("launchGoogle"),
                        "Sign in with Google",
                        icon = icon("google")
                    )
                ),
                ## div(
                ##     class = "social-button apple-button",
                ##     actionLink(
                ##         ns("launchApple"),
                ##         "Sign in with Apple",
                ##         icon = icon("apple")
                ##     )
                ## ),
                div(
                    class = "social-button facebook-button",
                    actionLink(
                        ns("launchFacebook"),
                        "Sign in with Facebook",
                        icon = icon("facebook")
                    )
                )
            )
        )
    }
    if(with.firebase_emailonly) {
        div.firebase <- div(
            class = "card",
            div(
                class = "card-body",
                h1(
                    "Sign in",
                    class = "card-title pb-2"
                ),
                div(
                    "Enter your email and we'll send you a link."
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

    div.alt <- div()
    if(!is.null(alt)) div.alt <- alt
    top <- shiny::HTML(rep("<br>",sum(c(!with.email,!with.username,!with.password))))

    div.button <- div(
        id="splash-buttons",
        class = "pb-4",
        actionButton(ns("login_btn"),login.text,class="btn-warning btn-xl shadow blink")
    )

    ##splash.panel=div();ns=function(x)x
    if(with.firebase || with.firebase_emailonly) {
        splash.content <- div.firebase
    } else {
        splash.content <- div(
            id="splash-login",
            ## top,
            div.username,
            div.email,
            div.password,
            div.alt,
            div.button
        )
    }

    splash.subtitle <- NULL

    body <- div(
      id="splash-content",
      splash.content
    )

    title <- div(
        div(id="splash-title",splash.title),
        div(id="splash-subtitle",splash.subtitle)
    )

    m <- splashScreen(title, body, ns=ns)

    return(m)
}

splashscreen.buttons <- function() {
    tagList(
        shiny::tags$a(
            shiny::img(
                id="splash-logo2",
                src="static/bigomics-logo.png"
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

splashScreen <- function(title, body, ns=NULL, easyClose=FALSE, fade=FALSE,
                         buttons=TRUE, footer=TRUE)
{

  if(is.null(ns)) ns <- function(e) return(e)
  message("[AuthenticationModule::monsterFullScreen]")

  div.buttons = shiny::modalButton("Dismiss")
  if(buttons) {
    div.buttons <- splashscreen.buttons()
  }
  if(!footer) {
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
            shiny::img(src="static/mascotte-sc.png", class = "img-fluid", id = "splash-image"),
        ),
        shiny::div(
            class = "col-md-3 offset-md-2",
            shiny::div(
                id="splash-panel",
                body,
                textOutput(ns("login_warning")),
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

