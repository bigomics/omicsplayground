##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

AuthenticationUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
}


NoAuthenticationModule <- function(input, output, session, show_modal=TRUE, username="", email="")
{
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
              login.text="Start!"
            )
            ## shinyjs::delay(2000, {output$login_warning <- shiny::renderText("")})
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

    observeEvent( input$firebaseLogout, {
        dbg("[NoAuthenticationModule] observe::input$firebaseLogout() reacted")
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

FirebaseAuthenticationModule <- function(input, output, session)
{
    message("[AuthenticationModule] >>>> using FireBase (email+password) authentication <<<<")

    dbg("[AuthenticationModule] getwd = ",getwd())
    dbg("[AuthenticationModule] file.exists('firebase.rds') = ",file.exists("firebase.rds"))

    if(file.exists("firebase.rds")) {
      firebase_config <- firebase:::read_config("firebase.rds")
    } else {
      stop("[ERROR] no firebase.rds file found. please create.")
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
        firebase$launch_google()
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

        ## sign out
        firebase$sign_out()
        
        m <- splashLoginModal(
            ns = ns,
            with.email = FALSE,
            with.username = FALSE,
            with.password = FALSE,
            with.register = FALSE,
            with.firebase = TRUE,            
            login.text = "Start!"
        )

        ## login modal
        shiny::showModal(m)
    }
    
    first_time = TRUE
    observeEvent(USER$logged, {
        
        ## no need to show the modalif the user is logged this is due
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

    observeEvent( input$firebaseLogout, {    
        message("[FirebaseAuthenticationModule] firebaseLogout triggered!")
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
        
        dbg("[FirebaseAuthenticationModule] observe::get_signed_in() reacted")

        response <- firebase$get_signed_in()
        
        if(!response$success) {
            dbg("[FirebaseAuthenticationModule] sign in NOT succesful")                        
            resetUSER() 
            return()
        } 

        on.exit({
            dbg("[FirebaseAuthenticationModule] get_signed_in() on.exit")            
            removeModal()      
        })
        
        dbg("[FirebaseAuthenticationModule] names(response) = ",names(response))        
        dbg("[FirebaseAuthenticationModule] names(response$response) = ",names(response$response))
        for(i in 1:length(response$response)) {
            dbg("[FirebaseAuthenticationModule] ",names(response$response)[i],"=",response$response[[i]])
        }

        t1=1637262031144
        t0 <- response$response[['createdAt']]
        t1 <- response$response[['lastLoginAt']]
        if(is.null(t0)) t0 <- 0
        if(is.null(t1)) t1 <- 0        
        t0 <- as.POSIXct(as.numeric(t0)/1000, origin="1970-01-01")
        t1 <- as.POSIXct(as.numeric(t1)/1000, origin="1970-01-01")           
        dbg("[FirebaseAuthenticationModule] createdAt=",t0)
        dbg("[FirebaseAuthenticationModule] lastLoginAt=",t1)
        dbg("[FirebaseAuthenticationModule] TIMEOUT=",TIMEOUT)        

        delta.secs <- as.numeric(Sys.time() - t1 , units='secs')
        delta.secs
        WAIT_TIME = 3600
        WAIT_TIME = TIMEOUT + 60*5
        WAIT_TIME = 60*5        

        ## NEED RETHINK!!! Not working very well because lastLoginAt
        ## is often the current login time. We would actually need the
        ## last logout time instead.
        ##
        if( FALSE && TIMEOUT>0 && delta.secs < WAIT_TIME ) {
            wait.mins <- format((WAIT_TIME - delta.secs)/60, digits=0)
            msg <- paste("You need to wait",wait.mins,"minutes before you can login again.")
            msg <- paste(msg,"\nLast login:",t1)
            shinyalert::shinyalert("Bummer...", msg, callbackR = resetUSER)            
            return()
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
        dbg("[FirebaseAuthenticationModule] observe::firebaseUpgrade reacted")        
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


##================================================================================
## PasswordAuthenticationModule (ask login.name + password)
##================================================================================

credentials.file='CREDENTIALS'
PasswordAuthenticationModule <- function(input, output, session,
                                         credentials.file)
{
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
            with.password=TRUE)
        ## shinyjs::delay(2000, {output$login_warning <- shiny::renderText("")})
        shiny::showModal(m)
        
    }

    CREDENTIALS <- read.csv(credentials.file,colClasses="character")
    head(CREDENTIALS)

    output$showLogin <- shiny::renderUI({
        m <- splashLoginModal(
            ns=ns,
            with.email=TRUE,
            with.username=FALSE,
            with.password=TRUE)
        ## shinyjs::delay(2000, {output$login_warning <- shiny::renderText("")})
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
            ##shinyjs::delay(2000, shinyjs::hide("login_warning", anim = TRUE, animType = "fade"))
            shinyjs::delay(2000, {output$login_warning <- shiny::renderText("")})
            USER$logged <- FALSE
        }
        ##hide("login_warning")
    })

    observeEvent( input$firebaseLogout, {
        dbg("[NoAuthenticationModule] observe::input$firebaseLogout() reacted")
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


##================================================================================
## HELPER FUNCTIONS
##================================================================================

splashHelloModal <- function(name, msg=NULL, ns=NULL, duration=3500)
{
    if(is.null(ns)) ns <- function(e) return(e)
    message("[AuthenticationModule::splashHelloModel]")
    
    all.hello = c("Hello","Salut","Hola","Privet","Ni hao","Ciao","Hi","Hoi","Hej",
                  "Yassou","Selam","Hey","Hei","Grutzi","Салам","Bonjour",
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
        shiny::div(shiny::HTML(title),style="font-size:70px;font-weight:700;line-height:1em;"),
        shiny::br(),
        shiny::div(shiny::HTML(subtitle),style="font-size:30px;"),
        shiny::br(),br(),br()
    )
    body <- shiny::tagList(
        shiny::div(id="splash-title",splash.title)
    )
    m <- splashScreen(body, ns=ns, easyClose=TRUE, fade=TRUE,
                              buttons=FALSE, footer=FALSE)    
    if(duration>0) {
        cat("closing hello in",round(duration/1000,1),"seconds...\n")
        shinyjs::delay(duration, shiny::removeModal())
    }
    return(m)
}

splashLoginModal <- function(ns=NULL, with.email=TRUE, with.password=TRUE,
                             with.username=FALSE, with.register=FALSE,
                             with.firebase=FALSE,
                             login.text="Login", alt=NULL)
{
    if(is.null(ns)) ns <- function(e) return(e)
    message("[AuthenticationModule::splashLoginModal]")

    titles <- list()
    titles[[1]] = c("Big Omics Data","Isn't big anymore with Omics Playground")
    titles[[2]] = c("Great Discoveries","Start on the Omics Playground")
    titles[[3]] = c("Fasten Your Seat Belts!","Hi-speed analytics")
    titles[[4]] = c("Do-it-yourself Omics Analytics","Yes you can!")
    titles[[5]] = c("Twenty-Four Seven","Your Playground doesn't go on coffee breaks")
    titles[[6]] = c("Analyze with confidence","Be a data rockstar, a Freddie Mercury of omics!")
    titles[[7]] = c("Play-Explore-Discover","Get deeper insights with Omics Playground")
    titles[[8]] = c("Skip the Queue","Take the fast lane. Self-service analytics.")
    titles[[9]] = c("Look Ma! No help!","I did it without a bioinformatician")
    titles[[10]] = c("Easy-peasy insight!","Get insight from your data the easy way")
    titles[[11]] = c("Zoom-zoom-insight!","Get faster insight from your data")
    titles[[12]] = c("Click-click-eureka!","Owe yourself that <i>eureka!</i> moment")
    titles[[13]] = c("I Love Omics Data!","Unleash your inner nerd with Omics Playground")
    titles[[14]] = c("More Omics Data","Is all I want for Christmas")
    titles[[15]] = c("Keep Exploring","Never stop discovering with Omics Playground")
    titles[[16]] = c("Real Bioinformaticians","Do it with Omics Playground")
    titles[[17]] = c("Real Biologists","Do it with Omics Playground")
    titles[[18]] = c("Ich bin doch nicht bl\u214d!","Of course I use Omics Playground")
    titles[[19]] = c("Non sono mica scemo!","Of course I use Omics Playground")
    ## below from https://www.quotesweekly.com/keep-exploring-quotes/
    titles[[20]] = c("The Unexplored Plan","When you get into exploring, you realize that we live on a relatively unexplored plan. &ndash; E. O. Wilson")
    titles[[21]] = c("Explore More","The more you explore, the more you learn and grow.<br>&ndash; Nitesh Nishad")
    titles[[22]] = c("Discover New Oceans","Man cannot discover new oceans unless he has the courage to lose sight of the shore. &ndash; Andre Gide")
    titles[[23]] = c("Adventurous Life","Love adventurous life. Be passionately curious about exploring new adventures. &ndash; Lailah Gifty Akita")
    titles[[24]] = c("Succes is Exploration","The first thing you have to find is the unknown. Learning is searching. Anything else is just waiting. &ndash; Dale Daute")
    titles[[25]] = c("Look Ma! No help!","I did it without a bioinformagician")
    titles[[26]] = c("May the Force of Omics be with you","Train hard youngling, one day a master you become")    
    title <- titles[[length(titles)]]
    title <- sample(titles,1)[[1]]
    title.len <- nchar(paste(title,collapse=' '))
    if(title.len < 80) title[1] <- paste0("<br>",title[1])
    ##if(title.len < 40) title[1] <- paste0("<br>",title[1])
    splash.title <- shiny::div(
        shiny::br(),br(),
        shiny::div(shiny::HTML(title[1]),style="font-size:70px;font-weight:700;line-height:1em;"),
        shiny::br(),
        shiny::div(shiny::HTML(title[2]),style="font-size:28px;"),
        shiny::br(),br(),br()
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
            id = "firebaseAuth",
            div(
                id = "firebaseBtns",
                actionButton(ns("launchGoogle"), "Google", icon = icon("google"), class = "btn-warning"),
                "  ",
                tags$a(class = "btn btn-default", onclick = "toggleEmail()", icon("envelope"), "Email"),
            ),
            div(
                id = "emailLinkWrapper",
                div(
                    class = "row",
                    div(
                        class = "col-md-9",
                        textInput(
                            ns("emailInput"),
                            "",
                            placeholder = "Your email",
                            width = "100%"
                        )
                    ),
                    div(
                        class = "col-md-3",
                        br(),
                        actionButton(
                            ns("emailSubmit"),
                            "Send"
                        )
                    )
                ),
                p(
                    id = "emailFeedbackShow",
                    class = "white"
                )
            )
        )
    }

    div.alt <- div()
    if(!is.null(alt)) div.alt <- alt
    top <- shiny::HTML(rep("<br>",sum(c(!with.email,!with.username,!with.password))))

    div.button <- div(
        id="splash-buttons",
        actionButton(ns("login_btn"),login.text,class="btn-warning")
    )
    if(with.register) {
      div.button <- div(
        id="splash-buttons",        
        actionButton(ns("login_btn"),login.text, class = "btn-outline-primary"),
        actionButton(ns("register_btn"),"Register",class="btn-primary")
      )
    }
    
    ##splash.panel=div();ns=function(x)x
    if(with.firebase) {
      splash.content <- div(
        "please login",
        div.firebase
      )
    } else {
      splash.content <- div(
        id="splash-login",
        br(),br(),top,
        div.username,
        div.email,
        div.password,
        div.alt,
        br(),
        div.button
        )
    }

    splash.subtitle <- NULL
    
    body <- div(
      id="splash-panel",
      div(id="splash-title",splash.title),
      div(id="splash-content",splash.content),      
      div(id="splash-subtitle",splash.subtitle)
    )

    m <- splashScreen(body, ns=ns)

    return(m)
}

splashscreen.buttons <- function() {
  shiny::div(
    shiny::tags$a(
      shiny::img(id="splash-logo2", src=base64enc::dataURI(file="www/bigomics-logo.png")),
      href = "https://www.bigomics.ch",
      target = "_blank"
    ),    
    shiny::tags$a(
      icon("book"),
      "Read-the-docs", 
      class = "btn btn-outline-primary",
      href = "https://omicsplayground.readthedocs.io",
      target = "_blank"
    ),
    shiny::tags$a(
      icon("youtube"),
      "Watch tutorials",
      class = "btn btn-outline-primary",
      href = "https://www.youtube.com/channel/UChGASaLbr63pxmDOeXTQu_A",
      target = "_blank"
    ),
    shiny::tags$a(
      icon("github"),
      "Get the source", 
      class = "btn btn-outline-primary",
      href = "https://github.com/bigomics/omicsplayground",
      target = "_blank"
    ),
    shiny::tags$a(
      icon("docker"),
      "Docker image", 
      class = "btn btn-outline-primary",
      href = "https://hub.docker.com/r/bigomics/omicsplayground",
      target = "_blank"
    ),
    shiny::tags$a(
      icon("users"),
      "User forum", 
      class = "btn btn-outline-primary",
      href = "https://groups.google.com/d/forum/omicsplayground",
      target = "_blank"
    ),
    shiny::tags$a(
      icon("coffee"),
      "Buy us a coffee!",
      class = "btn btn-outline-primary",
      href = "https://www.buymeacoffee.com/bigomics",
      target = "_blank"
    )
  )

}

splashScreen <- function(body, ns=NULL, easyClose=FALSE, fade=FALSE,
                         buttons=TRUE, footer=TRUE)
{
    
  if(is.null(ns)) ns <- function(e) return(e)
  message("[AuthenticationModule::monsterFullScreen]")
  
  div.footer = shiny::modalButton("Dismiss")
  if(buttons) {
    div.footer <- splashscreen.buttons()      
  }
  if(!footer) {
    div.footer <- NULL
  }
  
  ## return modalDialog
  m <- modalDialog2(
    id = "splash-fullscreen",
    class = "bg-primary",
    shiny::div(
      shiny::img(src=base64enc::dataURI(file="www/mascotte-sc.png"),id="splash-image"),
      body,
      shiny::br(),
      shiny::div(id="splash-warning",textOutput(ns("login_warning")),style="color:red;"),
      style="height: 32rem; width: 100%;"                
    ),
    footer = div(div.footer,id="splash-footer"),
    size = "fullscreen",
    easyClose = easyClose,
    fade = fade
  ) ## end of modalDialog
  
  return(m)
}


splashScreen.SAVE <- function(body, ns=NULL, easyClose=FALSE, fade=FALSE,
                              buttons=TRUE, footer=TRUE)
{
    
    if(is.null(ns)) ns <- function(e) return(e)
    message("[AuthenticationModule::splashModal]")

    div.footer = shiny::modalButton("Dismiss")
    if(buttons) {
        div.footer = bigomics.buttons()
    }
    if(!footer) {
        div.footer <- NULL
    }

    ## return modalDialog
    particlesjs.conf <- rjson::fromJSON(file="resources/particlesjs-config.json")
    m <- shiny::modalDialog(
        id = "modal-splash",
        class = "bg-primary",
        shiny::div(
            id="particles-target",
            shiny::div(id="splash-logo", shiny::img(src=base64enc::dataURI(file="www/logo.png"),
                                      width=32,height=32)),
            body,
            ##firebase::useFirebaseUI(),
            shiny::br(),
            shiny::div(id="splash-warning",textOutput(ns("login_warning")),style="color:red;"),
            style="height: 32rem; width: 100%;"                
        ),
        footer = div.footer,
        particlesjs::particles(config=particlesjs.conf, target_id = "particles-target", timeout = 1000),
        size="xl", easyClose=easyClose, fade=fade
    ) ## end of modalDialog
    
    return(m)
}
