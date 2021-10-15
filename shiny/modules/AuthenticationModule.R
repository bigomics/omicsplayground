##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

if(0) {
    type="password"
    CREDENTIALS = data.frame(username=c("ivo","stefan"),
                             password=c("iii","sss"),
                             expiry=c("2025-01-01","2020-01-01"),
                             row.names=1)
    credentials.file = "CREDENTIALS"
    .setSmtpServer("~/bigomics/server-conf/smtp_server/infomaniak.env")
    smtp.env = "~/bigomics/server-conf/smtp_server/infomaniak.env"
    smtp.env = "~/bigomics/server-conf/smtp_server/mailjet.env"
}

.setSmtpServer <- function(smtp.env="./smtp-server.env") {
    if(!file.exists(smtp.env)) return(NULL)
    settings <- read.table(smtp.env)[,2]
    smtp <- sapply( as.character(settings), strsplit, split="=")
    names(smtp) <- sapply(smtp,function(s) s[1])
    smtp
    do.call(Sys.setenv, lapply(smtp, function(s) s[2]))
    Sys.getenv()[grep("SMTP",names(Sys.getenv()))]
}

AuthenticationUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::showModal(shiny::uiOutput(ns("showLogin")))
}


NoAuthenticationModule <- function(input, output, session, username=NULL)
{
    message("[NoAuthenticationModule] >>>> using no authentication <<<<")
    ns <- session$ns    
    USER <- shiny::reactiveValues(logged=FALSE, name="", email="",
                                  password=NA, registered=NA,
                                  level="")    
    USER$name <- username
    
    output$showLogin <- shiny::renderUI({
        m <- splashLoginModal(
            ns=ns, with.email=FALSE, with.password=FALSE, login.text="Start")
        ## shinyjs::delay(2000, {output$login_warning <- shiny::renderText("")})
        shiny::showModal(m)
    })

    output$login_warning = shiny::renderText("")

    shiny::observeEvent( input$login_btn, {           
        shiny::removeModal()
        shiny::showModal(splashHelloModal(name=USER$name,ns=ns))
        USER$logged <- TRUE
        ##Sys.sleep(3);cat("wait 3 seconds to close...\n");removeModal()        
    })
    
    rt <- list(
        name   = shiny::reactive(USER$name),
        level  = shiny::reactive(USER$level),
        logged = shiny::reactive(USER$logged),
        limit  = shiny::reactive(USER$limit)
    )
    return(rt)
}
##================================================================================
## FirebaseAuthenticationModule
##================================================================================

FirebaseAuthenticationModule <- function(input, output, session)
{
    message("[FirebaseAuthenticationModule] >>>> using FireBase (email+password) authentication <<<<")

    ns <- session$ns    
    USER <- shiny::reactiveValues(
        logged = FALSE, 
        name = NA, 
        password = NA, 
        email = NA, 
        level = NA
    )    

    firebase <- firebase::FirebaseUI$
        new(persistence = "session")$ # instantiate
        set_providers( # define providers
            email_link = TRUE, 
            ##email = TRUE,
            google = TRUE
        )
    firebase$set_tos_url("https://bigomics.ch/terms")
    firebase$set_privacy_policy_url("https://bigomics.ch/privacy")    
    
    output$showLogin <- shiny::renderUI({

        message("[FirebaseAuthenticationModule] showLogin... ")
        
        m <- splashLoginModal(
            ns = ns,
            with.email = FALSE,
            with.username = FALSE,
            with.password = FALSE,
            with.register = FALSE,
            with.firebase = TRUE,            
            login.text = "Start!"
        )

        on.exit({
            dbg("[FirebaseAuthenticationModule] on.exit")            
            firebase$launch()
        })

        # no need to show the modal
        # if the user is logged
        # this is due to persistence
        if(USER$logged) {
            dbg("[FirebaseAuthenticationModule] USER is already logged in! no modal")                        
            return()
        }

        dbg("[FirebaseAuthenticationModule] showing Firebase login modal")                                
        shiny::tagList(
                   shiny::showModal(m)
               )
    })
    
    observeEvent( firebase$get_signed_in(), {

        response <- firebase$get_signed_in()

        dbg("[FirebaseAuthenticationModule] get_signed_in() reacted")
        dbg("[FirebaseAuthenticationModule] response$success = ",response$success)
        
        if(!response$success) {
            dbg("[FirebaseAuthenticationModule] sign in NOT succesful")                        
            return()
        } else {
            dbg("[FirebaseAuthenticationModule] sign in SUCCESSFUL!")            
        }

        on.exit({
            dbg("[FirebaseAuthenticationModule] get_signed_in() on.exit")            
            removeModal()            
        })
        
        USER$logged <- TRUE
        USER$name <- response$response$displayName
        USER$email <- response$response$email
    })
    
    rt <- list(
        name   = shiny::reactive(USER$name),
        level  = shiny::reactive(USER$level),
        logged = shiny::reactive(USER$logged),
        limit = shiny::reactive(USER$limit)
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
    message("[NoAuthenticationModule] >>>> using local Email+Password authentication <<<<")

    ns <- session$ns    
    USER <- shiny::reactiveValues(logged=FALSE, username=NA, password=NA, level=NA, limit=NA)    
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
        
        login_username <- input$login_username
        if(is.null(login_username) || login_username =="") {
            login_username <- input$login_email
        }
        login_password <- input$login_password

        if( is.null(login_username) || is.null(login_password)) return(NULL)
        if( login_username=="" || login_password=="") return(NULL)            
        sel <- tail(which( CREDENTIALS$username == login_username),1)       
        valid.user <- isTRUE(login_username %in% CREDENTIALS$username)
        valid.pw   <- isTRUE(CREDENTIALS[sel,"password"]==input$login_password)
        valid.date <- isTRUE(Sys.Date() < as.Date(CREDENTIALS[sel,"expiry"]) )
        login.OK = (valid.user && valid.pw && valid.date)
        
        message("--------- password login ---------")
        message("input.username = ",input$login_username)
        message("input.email    = ",input$login_email)
        message("input.password = ",input$login_password)
        message("user.password  = ",CREDENTIALS[sel,"password"])
        message("user.expiry    = ",CREDENTIALS[sel,"expiry"])
        message("user.name      = ",CREDENTIALS[sel,"username"])
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
            USER$username  <- cred$username
            USER$level     <- cred$level
            USER$limit     <- cred$limit
            ##USER$expiry    <- cred$expiry
            ##USER$password  <- cred$password
            
            ## Here you can perform some user-specific functions, or site news
            shiny::showModal(splashHelloModal(USER$name,ns=ns))
            ##removeModal()
            USER$logged <- TRUE            

        } else {

            message("[PasswordAuthenticationModule::login] REFUSED : invalid login! ")

            if(!valid.date) {
                output$login_warning = shiny::renderText("Registration expired")
            } else {
                output$login_warning = shiny::renderText("Invalid username or password")
            }
            ##shinyjs::delay(2000, shinyjs::hide("login_warning", anim = TRUE, animType = "fade"))
            shinyjs::delay(2000, {output$login_warning <- shiny::renderText("")})
            USER$logged <- FALSE
        }
        ##hide("login_warning")
    })


    ## module reactive return value
    rt <- list(
        name   = shiny::reactive(USER$username),
        level  = shiny::reactive(USER$level),
        logged = shiny::reactive(USER$logged),
        limit  = shiny::reactive(USER$limit)        
    )
    return(rt)
}

##================================================================================
## RegisterAuthenticationModule (just ask email, no password)
##================================================================================

register.file="../logs/register.log"
RegisterAuthenticationModule <- function(input, output, session, register.file)
{
    message("[NoAuthenticationModule] >>>> using Register authentication <<<<")

    ns <- session$ns
    dir.create("../logs",showWarnings=FALSE)        
    register.file    
    if(!file.exists(register.file)) {
        ee <- data.frame(email="", name="", password="",
                         date=NA, expiry=NA, registered=FALSE, level="basic")
        write.table(ee[0,], file=register.file, append=FALSE, sep=",",
                    row.names=FALSE, col.names=TRUE)
    }
    REGISTERED <- read.csv(register.file, colClasses="character",
                            stringsAsFactors=FALSE)
    head(REGISTERED)    
    SURVEY.LOG = "../logs/survey.log"
    DIALOG_SLEEP = 0.3
    USER <- shiny::reactiveValues(logged=FALSE, name="", email="",
                           password=NA, registered=FALSE, level="")
        
    output$showLogin <- shiny::renderUI({
        message("[AuthenticationModule::UI] USER$logged = ",USER$logged)
        ## If login is OK then do nothing
        if(USER$logged) {
            return(NULL)
        } else {
            showLoginDialog()             
        }
    })

    showLoginDialog <- function() {
        ##removeModal()
        alt <- shiny::actionLink(ns("create_account"),"or create an account",
                          style="color:white;", class="white-link")
        shiny::showModal(splashLoginModal(
            ns=ns,
            with.password = FALSE,
            with.username = FALSE,
            with.email = TRUE,
            alt = alt))
    }
    
    shiny::observeEvent( input$login_btn, {           
        
        message("[AuthenticationModule::login] login_btn pressed")                
        message("[AuthenticationModule::login] email = ",input$login_email)
        email <- input$login_email
        ## check login
        email.exists <- email %in% REGISTERED$email
        if(!email.exists) {
            output$login_warning = shiny::renderText("invalid email")
            shinyjs::delay(2000, {output$login_warning <- shiny::renderText("")}) 
            return(NULL)
        }

        sel <- tail(which(REGISTERED$email == email),1)
        user.registered <- REGISTERED[sel,"registered"]
        USER$name <- REGISTERED[sel,"name"]
        USER$email <- email
        USER$registered <- user.registered
        USER$logged <- TRUE
        
        message("[AuthenticationModule::login] email.exists = ",email.exists)
        message("[AuthenticationModule::login] user.registered = ",user.registered)
        shiny::removeModal()
        Sys.sleep(DIALOG_SLEEP)
        shiny::showModal(splashHelloModal(
            ns = ns,
            name = USER$name,
            msg = "Welcome back. We wish you many great discoveries today!"
        ))
        
    })
        
    shiny::observeEvent( input$create_account, {
        shiny::removeModal()
        Sys.sleep(DIALOG_SLEEP)
        m <- createAccountModal(USER)
        shiny::showModal(m)
    })
    
    validEmailFormat <- function(email) {
        has.atsign <- grep("@",email)
        domain <- strsplit(email,split='@')[[1]][2]
        has.dot <- length(strsplit(domain,split="[.]")[[1]]) > 1
        has.atsign && has.dot
    }

    createAccountModal <- function(user) {
        message("[AuthenticationModule::createAccount] user$name = ",user$name)            

        m <- shiny::modalDialog(
            id = "auth_dialog",
            ##title = "Create an Account",
            shiny::tagList(
                shiny::fillRow(
                    flex = c(0.12,0.5,0.2,0.5,0.15),
                    height = 450,
                    shiny::br(),
                    shiny::tagList(
                        shiny::h3("Try premium?"),
                        shiny::HTML("Good news! For limited time, all registered users will receive free access to our <b>Premium plan</b> which includes:<br><br>"),
                        shiny::HTML("<ul><li>Drug enrichment<li>Biomarker analysis<li>Connectivity mapping<li>WGCNA analysis<li>And more!</ul>"),
                        shiny::br(),
                        ##h4("Recommend us to your friends"),
                        ##textAreaInput(ns("register_friends"),"Email(s) of friends:",rows=3)
                        ),
                    shiny::br(),
                    shiny::tagList(
                        ##h3("Tell us a bit about yourself"),
                        shiny::h3("Create a free account"),                            
                        shiny::textInput(ns("register_name"),"First name"),
                        shiny::textInput(ns("register_lastname"),"Last name"),
                        shiny::textInput(ns("register_email"),"E-mail",value=USER$email),
                        shiny::textInput(ns("register_organization"),"Organization"),
                        shiny::textInput(ns("register_jobtitle"),"Job title"),                        
                        shiny::selectInput(ns("register_institutiontype"),"Institution type",
                                    ## select=FALSE, selectize=FALSE,
                                    c("",sort(c("Pharma/Biotech",
                                                "Hospital/Medical Center",
                                                "Government",
                                                "Service provider",
                                                "Start-up")),"other")),
                        shiny::selectInput(ns("register_hear"),"How did you hear about us?",
                                    c("",sort(c("From coworker/friend",
                                                "LinkedIn",
                                                "Twitter",
                                                "Journal article",
                                                "Web article/blog",
                                                "Web search",
                                                "BigOmics website",
                                                "Event",
                                                "Contacted by Sales")), "other")),
                        shiny::br(),br()
                    ),
                    shiny::br()
                )
            ),
            footer = shiny::tagList(
                shiny::div(shiny::textOutput(ns("register_warning")),style="color:red; text-align:center;"),
                shiny::actionLink(ns("register_cancel"), "Cancel"),shiny::HTML("&nbsp;&nbsp;"),
                shiny::actionButton(ns("register_submit"), "Register")
            ),                    
            ##footer = NULL,
            easyClose = FALSE,
            size = "m"
        )

        shiny::observeEvent( input$register_cancel, {
            shiny::removeModal()	
            Sys.sleep(DIALOG_SLEEP)
            showLoginDialog() 
        }, once=TRUE, ignoreInit=TRUE)

        shiny::observeEvent( input$register_submit, {
            
            message("[AuthenticationModule::createAccount] REGISTER pressed")
            message("[AuthenticationModule::createAccount] REGISTER pressed")
            qq <- c(
                name = input$register_name,
                lastname = input$register_lastname,
                email = input$register_email,
                organization = input$register_organization,                
                jobtitle = input$register_jobtitle,
                organization = input$register_organization,                    
                institutiontype = input$register_institutiontype,                    
                hear = input$register_hear
            )
            
            if(any(is.null(qq) | qq=="")) {
                message("[AuthenticationModule::createAccount] missing fields")
                output$register_warning = shiny::renderText("Please fill in all required fields")
                shinyjs::delay(2000, {output$register_warning <- shiny::renderText("")})
                return(NULL)
            }

            email.ok = validEmailFormat( input$register_email )
            if(!email.ok) {
                output$register_warning = shiny::renderText("Error: invalid email format")
                shinyjs::delay(2000, {output$register_warning <- shiny::renderText("")}) 
                return(NULL)
            }

            email.exists <- (input$register_email %in% REGISTERED$email)
            if(email.exists) {
                output$register_warning = shiny::renderText("Error: email exists")
                shinyjs::delay(2000, {output$register_warning <- shiny::renderText("")}) 
                return(NULL)
            }

            message("[AuthenticationModule::createAccount] all answered!")
            lapply(1:length(qq), function(i) message(names(qq)[i],"=",qq[i]))
            ##
            ## save somewhere!
            ##
            USER$name  <- input$register_name
            USER$email <- input$register_email
            USER$password <- ""
            USER$registered <- TRUE
            USER$level <- "free"
            updateRegister(USER, register.file)
            updateSurveyLog(USER, qq)                
            USER$logged = TRUE

            ## success
            shiny::removeModal()
            shiny::showModal(splashHelloModal(
                ns = ns,                
                name = USER$name,
                msg = "Thank you for registering. We wish you many great discoveries today!"
            ))
            
        })
        m
    }   

    updateRegister <- function(user, register.file) {
        ##today = as.character(Sys.Date())
        today = as.character(Sys.time())
        ##expiry <- as.character(Sys.Date()+365)
        expiry <- "2099-01-01"
        ee <- data.frame(
            email = user$email,
            name = user$name,
            password = user$password,
            date = today,
            expiry = expiry,
            registered = user$registered,
            level = user$level)

        if(nrow(REGISTERED)>0) {
            ee <- ee[,match(colnames(REGISTERED),colnames(ee))]            
            ee
            REGISTERED <- rbind(REGISTERED, ee)
        } else {
            REGISTERED <- ee
        }
        suppressWarnings(
            write.table(ee, file=register.file, append=TRUE, sep=",",
                        row.names=FALSE, col.names=!file.exists(register.file))
        )
    }

    updateSurveyLog <- function(user, answers) {
        ##qa.list <- list(...)
        today = as.character(Sys.time())        
        QA <- cbind(question=names(answers), answer=answers)
        ee <- data.frame( email = user$email, date = today, QA)
        rownames(ee) <- NULL        
        suppressWarnings(
            write.table(ee, file=SURVEY.LOG, append=TRUE, sep=",",
                        row.names=FALSE,
                        col.names=!file.exists(SURVEY.LOG))
        )        
    }    
    ##hide("login_warning")
    output$login_warning = shiny::renderText("")
    output$register_warning = shiny::renderText("")
    res <- list(
        name   = shiny::reactive(USER$name),
        level  = shiny::reactive(USER$level),
        logged = shiny::reactive(USER$logged),
        limit = shiny::reactive(USER$limit)
    )
    return(res)
}


##================================================================================
## EmailAuthenticationModule (ask login.name + password + CAPTCHA verification)
##================================================================================

credentials.file="./CREDENTIALS"
EmailAuthenticationModule.SAVE <- function(input, output, session,
                                           credentials.file,
                                           type="splash")
{
    message("[NoAuthenticationModule] >>>> using Email authentication <<<<")

    ns <- session$ns
    dir.create("../logs",showWarnings=FALSE)            
    credentials.file    
    if(!file.exists(credentials.file)) {
        ee <- data.frame(email="", name="", password="",
                         date=NA, expiry=NA, registered=FALSE, level="basic")
        write.table(ee[0,], file=credentials.file, append=FALSE, sep=",",
                    row.names=FALSE, col.names=TRUE)
    }
    CREDENTIALS <- read.csv(credentials.file, colClasses="character",
                            stringsAsFactors=FALSE)
    head(CREDENTIALS)    
    SURVEY.LOG = "../logs/survey.log"
    USER <- shiny::reactiveValues(logged=FALSE, name="", email="",
                           password=NA, verified=FALSE, level="")
    
    SMTP_SERVER  <- Sys.getenv("SMTP_SERVER")
    SMTP_USER    <- Sys.getenv("SMTP_USER")
    SMTP_SECRET  <- Sys.getenv("SMTP_SECRET")
    has.smtp <- SMTP_SERVER!="" && SMTP_USER!="" && SMTP_SECRET!=""
    has.smtp
    message("[AuthenticationModule] has.smtp = ",has.smtp)
    
    DIALOG_SLEEP = 0.3
    ##waiter_show(spin_fading_circles())
    waiter <- Waiter$new(id="animated_captcha", color=transparent(0.5))

    validEmailFormat <- function(email) {
        has.atsign <- grep("@",email)
        domain <- strsplit(email,split='@')[[1]][2]
        has.dot <- length(strsplit(domain,split="[.]")[[1]]) > 1
        has.atsign && has.dot
    }
    
    output$showLogin <- shiny::renderUI({
        message("[AuthenticationModule::UI] USER$logged = ",USER$logged)
        ## If login is OK then do nothing
        if(USER$logged) {
            return(NULL)
        } else {
            showDialog("login")
        }
    })
        
    shiny::observeEvent( input$login_btn, {           
        
        message("[AuthenticationDialog::login] login_btn pressed")                
        message("[AuthenticationModule::login] email = ",input$login_email)
        email <- input$login_email

        ## check login
        valid.email <- validEmailFormat(email)
        if(!valid.email) {
            message("[AuthenticationModule::login] email not valid")
            output$login_warning = shiny::renderText("invalid email format")
            shinyjs::delay(2000, {output$login_warning <- shiny::renderText("")}) 
            return(NULL)
        }
        
        email.exists <- email %in% CREDENTIALS$email
        user.registered <- FALSE
        nlogin <- sum(CREDENTIALS$email == email)
        ##USER <- shiny::reactiveValues(logged=FALSE, name=NA, email=NA, verified=FALSE)
        if(email.exists) {
            sel <- tail(which(CREDENTIALS$email == email),1)
            user.registered <- CREDENTIALS[sel,"registered"]
            USER$name <- CREDENTIALS[sel,"name"]
        }
        USER$email <- email
        USER$registered <- user.registered
        USER$logged <- FALSE
        
        message("[AuthenticationModule::login] email.exists = ",email.exists)
        message("[AuthenticationModule::login] user.registered = ",user.registered)
        message("[AuthenticationModule::login] nlogin = ",nlogin)
        
        if(email.exists && user.registered) {
            ## existing user, already registered
            for(i in 1:10) shiny::removeModal()
            Sys.sleep(DIALOG_SLEEP)
            updateCredentials(USER, credentials.file)
            if(FALSE && nlogin == 2) {
                message("[AuthenticationModule::login] >>> asking user experience")
                showDialog("ask.experience")
            } else {
                message("[AuthenticationModule::login] >>> welcome  user")                
                showDialog("welcome")
            }
            USER$logged <- TRUE
        } else  {
            ## fist time or unregistered user
            message("[AuthenticationModule::login] >>> creating account")
            showDialog("create.account")
        }
    })

    verifyUser <- function() {
        auth.code <- paste0(sample(c("A","T","G","C"), 7, replace=TRUE),collapse="")
        auth.code        
        USER$password <- auth.code
        if(has.smtp) {
            smtpAuthCode(email, auth.code)
            showDialog("verify")
        } else {            
            showDialog("captcha", auth.code)
        }
        auth.code
    }
    
    ##observeEvent( input$verify_recaptcha, {})    
    shiny::observeEvent( input$verify_btn, {
        message("[AuthenticationModule::register] verify_btn pressed")
        message("[AuthenticationModule::register] USER$email = ",USER$email)
        message("[AuthenticationModule::register] USER$password = ",USER$password)        
        message("[AuthenticationModule::register] verify_code = ",input$verify_code)

        if(input$verify_code == "") return(NULL)

        ok = (USER$password == input$verify_code)
        message("[AuthenticationModule::register] ok = ",ok)        
        
        if(ok) {
            output$verify_warning = shiny::renderText("")
            output$captcha_warning = shiny::renderText("")
            ##removeModal()	
            for(i in 1:10) shiny::removeModal()
            Sys.sleep(DIALOG_SLEEP)
            updateCredentials(USER, credentials.file)
            showDialog("welcome")
            USER$logged <- TRUE
        } else {
            output$verify_warning = shiny::renderText("Invalid verification")
            output$captcha_warning = shiny::renderText("Invalid verification")
            shinyjs::delay(2000, {output$verify_warning <- shiny::renderText("")})
            shinyjs::delay(2000, {output$captcha_warning <- shiny::renderText("")})            
        }        
    })

    shiny::observeEvent( input$cancel_btn, {
        shiny::removeModal()	
        Sys.sleep(DIALOG_SLEEP)
        showDialog("login")
    })
    
    ##---------------------------------------------------
    ## functions
    ##---------------------------------------------------
    
    showDialog <- function(stage, data=NULL)
    {
        message("[AuthenticationModule::showDialog] stage = ",stage)
        ##removeModal()

        output$verify_warning = shiny::renderText("")
        output$login_warning = shiny::renderText("")
        output$register_warning = shiny::renderText("")

        m <- NULL
        if(stage=="login" && type=="simple") {
            m <- shiny::modalDialog(
                id = "auth_dialog",
                title = "Log in to your Playground",
                shiny::tagList( 
                    shiny::p( "Enter your email to log in"),
                    shiny::textInput(ns("login_email"), "E-mail:"),
                    shiny::div(shiny::actionButton(ns("login_btn"), "Login"),style="text-align: center;")
                ),
                footer = NULL,
                size = "s"
            )
        } else if(stage=="login" && type=="splash") {
            m <- splashLoginModal(ns=ns, with.password=FALSE)            
        } else if(stage=="verify") {
            m <- shiny::modalDialog(
                id = "auth_dialog",
                title = "Verify your email",
                shiny::tagList( 
                    shiny::p("Please enter the verification code that we sent to your email:"),
                    shiny::textInput(ns("verify_code"), "Verification code:"),
                    shiny::div( shiny::actionButton(ns("cancel_btn"), "Cancel"),
                        shiny::actionButton(ns("verify_btn"), "Submit"),
                        style="text-align: center;")
                ),
                footer = shiny::div(shiny::textOutput(ns("verify_warning")),style="color: red;"),
                size = "s"
            )
        } else if(stage=="captcha") {            
            m <- shiny::modalDialog(
                id = "auth_dialog",
                title = "Please verify you're human",
                shiny::tagList( 
                    shiny::p("Copy the amino acid letters into the captcha box below:"),
                    ##div(shiny::img(src=base64enc::dataURI(file=tmpfile), width="240px"), style="text-align: center;"),
                    shiny::imageOutput(ns("animated_captcha"),height=50,width=240),
                    shiny::br(),
                    shiny::textInput(ns("verify_code"), "Enter captcha:"),                    
                    shiny::div(
                        shiny::actionButton(ns("verify_recaptcha"), "Refresh captcha"), 
                        shiny::actionButton(ns("verify_btn"), "Submit"), style="text-align: center;")
                ),
                footer = shiny::div(shiny::textOutput(ns("captcha_warning")),style="color: red;"),                
                size = "s"
            )
            ##unlink(tmpfile)
        } else if(stage=="welcome") {
            m <- splashHelloModal(USER$name, ns=ns)            
        } else if(stage=="create.account") {
            m <- createAccountModal(USER)
        } else if(stage=="ask.experience") {
            m <- askExperienceModal(USER)             
        } else {
            m <- shiny::modalDialog(
                id = "auth_dialog",
                title = "Error",
                shiny::tagList( shiny::p("error.") ),
                footer = shiny::div(shiny::HTML("error"),style="color: red;"),
                size = "s"
            )
        }
        shiny::showModal(m)
    }

    animateCaptcha <- function() {
        ##graphics::text(0.5,0.5,data,cex=4)
        waiter$show()
        n = 7
        dd <- rep(NA,n)
        xx <- seq(0,1,1/length(dd))*1.03
        aa <- rep(0,length(dd))
        par(mar=c(0,0,0,0))
        empty <- function(i) {
            base::plot(0,0,pch="",xlim=c(-0.1,0.95),ylim=c(0,1),bty="n",axes=FALSE)
            if(i>1) for(i in 1:(i-1)) graphics::text(xx[i],0.5,dd[i],srt=aa[i],cex=4)
        }
        for(i in 1:length(dd)) {
            for(j in 1:4) {
                empty(i)
                a <- runif(1,-35,35)
                d <- sample(c("A","T","C","G"), 1, replace=TRUE)                
                graphics::text(xx[i],0.5,d,srt=a,cex=4)
                ##Sys.sleep(0.01)
            }
            aa[i] <- a
            dd[i] <- d
        }
        USER$password <- paste(dd,collapse="")        
        return(USER$password)
    }
    
    output$animated_captcha <- shiny::renderImage({
        input$verify_recaptcha

        tmpfile="./captcha.png"
        tmpfile <- tempfile()        
        message("[AuthenticationModule:animated_captcha] render GIF")        
        saveGIF( animateCaptcha(),
                movie.name = tmpfile, interval = 0.01, fps=70,
                nmax = 100, ani.width = 240, ani.height = 50, loop=1)
        ## Return a list containing information about the image
        message("[AuthenticationModule:animated_captcha] saved GIF")
        list(
            src = tmpfile,
            contentType = "image/gif",
            width = 240,
            height = 50,
            alt = "animated captcha"
        )
    })
    
    createAccountModal.NOTUSED <- function(user) {
        message("[AuthenticationModule::createAccount] user$name = ",user$name)            
        m <- shiny::modalDialog(
            id = "auth_dialog",
            ##title = "Tell us a bit about yourself",
            shiny::tagList(
                shiny::fillRow(
                    flex = c(0.02,0.4,0.03,0.6),
                    height = 380,
                    shiny::br(),
                    shiny::tagList(
                        shiny::h3("Try premium?"),
                        shiny::HTML("Good news! For limited time, all registered users will receive free trial of our <b>Premium plan</b> which includes:<br><br>"),
                        shiny::HTML("<ul><li>Drug enrichment<li>Biomarker analysis<li>Connectivity mapping<li>WGCNA analysis<li>And more!</ul>"),
                        shiny::br(),
                        ##h4("Recommend us to your friends"),
                        ##textAreaInput(ns("ask_friends"),"Email(s) of friends:",rows=3)
                        ),
                    shiny::br(),
                    shiny::tagList(
                        shiny::h3("Tell us a bit about yourself"),                            
                        shiny::textInput(ns("ask_name"),"First name"),
                        shiny::textInput(ns("ask_lastname"),"Last name"),
                        shiny::textInput(ns("ask_email"),"E-mail",value=USER$email),
                        shiny::textInput(ns("ask_organization"),"Organization"),
                        shiny::textInput(ns("ask_jobtitle"),"Job title"),                        
                        shiny::selectInput(ns("ask_institutiontype"),"Institution type?",
                                    ## select=FALSE, selectize=FALSE,
                                    c("",sort(c("Pharma/Biotech",
                                                "Hospital/Medical Center",
                                                "Government",
                                                "Service provider",
                                                "Start-up")),"other")),
                        shiny::selectInput(ns("ask_hear"),"How did you hear about us?",
                                    c("",sort(c("From coworker/friend",
                                                "LinkedIn",
                                                "Twitter",
                                                "Journal article",
                                                "Web article/blog",
                                                "Web search",
                                                "BigOmics website",
                                                "Event",
                                                "Contacted by Sales")), "other"))
                    ),
                    shiny::br(),br()
                )
            ),
            footer = shiny::tagList(
                shiny::div(shiny::textOutput(ns("ask_warning")),style="color: red;"),
                shiny::actionLink(ns("ask_btn_skip"), "Skip"),shiny::HTML("&nbsp;&nbsp;"),
                shiny::actionButton(ns("ask_register"), "Register")
            ),                    
            ##footer = NULL,
            easyClose = FALSE,
            size = "m"
        )
        shiny::observeEvent( input$ask_btn_skip, {
            shiny::removeModal()
        })
        shiny::observeEvent( input$ask_register, {
            message("[AuthenticationModule::createAccount] REGISTER pressed")
            questions <- c(
                name = input$ask_name,
                lastname = input$ask_lastname,
                email = input$ask_email,
                organization = input$ask_organization,                
                jobtitle = input$ask_jobtitle,
                organization = input$ask_organization,                    
                institutiontype = input$ask_institutiontype,                    
                hear = input$ask_hear
            )

            if(any(is.null(questions) | questions=="")) {
                message("[AuthenticationModule::createAccount] missing fields")
                output$ask_warning = shiny::renderText("Please fill in all required fields")
                shinyjs::delay(2000, {output$ask_warning <- shiny::renderText("")})                
                return(NULL)
            } else {
                message("[AuthenticationModule::createAccount] all answered!")
                lapply(1:length(questions),
                       function(i) message(names(questions)[i],"=",questions[i]))
                ##
                ## save somewhere!
                ##
                USER$name <- input$ask_name
                updateCredentials(USER, credentials.file)                
                updateSurveyLog(USER, questions)                
            }
            USER$logged = TRUE
        })
        m
    }   
    
    askExperienceModal <- function(user) {
            message("[AuthenticationModule::askExperience]")            
            m <- shiny::modalDialog(
                id = "auth_dialog",
                ##title = "Tell us a bit about yourself",
                shiny::tagList(
                    shiny::fillRow(
                        flex = c(0.02,0.6,0.04,0.33),
                        height = 380,
                        shiny::br(),
                        shiny::tagList(
                            shiny::h3("Tell us how you like us"),br(),
                            shiny::radioButtons(
                                ns("askxp_satisfied"),
                                "Are you so far satisfied with your Playground experience?",
                                choiceValue=c(1,2,3,4,5), inline=TRUE,
                                choiceNames=c("no","so-so","ok","satisfied","love it!")),
                            shiny::radioButtons(
                                ns("askxp_recommend"),"Would you recommend us to a friend",
                                choiceValue=c(1,2,3),
                                choiceNames=c("not really","not sure","definitely!"),
                                inline=TRUE),
                            shiny::radioButtons(
                                ns("askxp_useagain"),"Would you use us again?",
                                c("no","depends","yes"), inline=TRUE),
                            shiny::selectInput(
                                ns("askxp_problems"),"What was your biggest problem?",
                                c("Data preparation","Slow website","Server disconnects",
                                  "Too complex", "No problems","Other"))
                            ##selectInput(ns("ask_cosize"),"How big is the company?",
                            ##            c("micro (<10)","small (<50)","medium (<250)",
                            ##              "large (> 250")),
                        ), shiny::br(),
                        shiny::tagList(
                            shiny::h3("Try premium?"),br(),
                            shiny::HTML("This is your chance! Complete this survey and unlock a free trial of our <b>Premium plan</b> which includes:<br><br>"),
                            shiny::HTML("<ul><li>Drug enrichment<li>Biomarker analysis<li>Connectivity mapping<li>And more!</ul>"),
                            shiny::br(),
                            ##h4("Recommend us to your friends"),
                            ##textAreaInput(ns("ask_friends"),"Email(s) of friends:",rows=3)
                        )
                    ),
                    shiny::div(
                        shiny::actionLink(ns("askxp_btn_skip"), "Skip"),shiny::HTML("&nbsp;&nbsp;"),
                        shiny::actionButton(ns("askxp_submit"), "Submit"),
                        style="text-align: right;"
                    )                    
                ),
                footer = NULL,
                easyClose = FALSE,
                size = "m"
            )

            shiny::observeEvent( input$askxp_btn_skip, {
                shiny::removeModal()
            })
            shiny::observeEvent( input$askxp_submit, {
                message("[AuthenticationModule::askExperience] Submit pressed")
                answers <- c(
                    askxp_satisfied = input$askxp_satisfied,
                    askxp_recommend = input$askxp_recommend,
                    askxp_useagain = input$askxp_useagain,
                    askxp_problems = input$askxp_problems
                )                    
                updateSurveyLog(USER, answers)
                shiny::removeModal()
            })
            m
    }   
        
    updateCredentials <- function(user, credentials.file) {
        ##today = as.character(Sys.Date())
        today = as.character(Sys.time())
        ##expiry <- as.character(Sys.Date()+365)
        expiry <- "2099-01-01"
        ee <- data.frame(
            email = user$email,
            name = user$name,
            password = user$password,
            date = today,
            expiry = expiry,
            registered = user$registered,
            level = user$level)

        if(nrow(CREDENTIALS)>0) {
            ee <- ee[,match(colnames(CREDENTIALS),colnames(ee))]            
            ee
            CREDENTIALS <- rbind(CREDENTIALS, ee)
        } else {
            CREDENTIALS <- ee
        }
        suppressWarnings(
            write.table(ee, file=credentials.file, append=TRUE, sep=",",
                        row.names=FALSE, col.names=!file.exists(credentials.file))
        )
    }

    updateSurveyLog <- function(user, answers) {
        ##qa.list <- list(...)
        today = as.character(Sys.time())
        
        QA <- cbind(question=names(answers), answer=answers)
        ee <- data.frame(
            email = user$email,
            date = today,
            QA
        )
        rownames(ee) <- NULL        
        suppressWarnings(
            write.table(ee, file=SURVEY.LOG, append=TRUE, sep=",",
                        row.names=FALSE,
                        col.names=!file.exists(SURVEY.LOG))
        )        
    }
    
    ##auth.code="GATTACA";email="ivo.kwee@gmail.com";name="Ivo Kwee"
    smtpAuthCode <- function(email, auth.code) {


        
        message("[AuthenticationModule::smtpAuthCode] *************************")
        message("[AuthenticationModule::smtpAuthCode] ****** SENDING EMAIL ****")
        message("[AuthenticationModule::smtpAuthCode] *************************")
        message("[AuthenticationModule::smtpAuthCode] email = ",email)
        message("[AuthenticationModule::smtpAuthCode] code = ",auth.code)

        SMTP_SERVER  <- Sys.getenv("SMTP_SERVER")
        SMTP_USER    <- Sys.getenv("SMTP_USER")
        SMTP_SECRET  <- Sys.getenv("SMTP_SECRET")        
        message("[AuthenticationModule::smtpAuthCode] SMTP_SERVER = ",SMTP_SERVER)
        message("[AuthenticationModule::smtpAuthCode] SMTP_USER = ",SMTP_USER)
        message("[AuthenticationModule::smtpAuthCode] SMTP_SECRET = ",SMTP_SECRET)        

        msg <- paste0(
            "Welcome to the Omics Playground!\n\n",            
            "Your login email is:\t",email,"\n",
            "Your authentication code is:\t",auth.code,"\n\n"
        )
        
        email.msg <- emayili::envelope() %>%
            from("no-reply@bigomics.ch") %>%
            to(email) %>%
            subject("Omics Playground authentication code") %>%
            text(msg)
        
        smtp <- emayili::server(
                             host = SMTP_SERVER,
                             port = 25,
                             username = SMTP_USER,
                             password = SMTP_SECRET)
        
        ##smtp(email.msg, verbose=FALSE)
        smtp(email.msg, verbose=TRUE)

    }
    
    ##hide("login_warning")
    output$login_warning = shiny::renderText("")

    ## shiny::observeEvent( input$logout, {
    ##     ##updateTextInput(session, ".username", value=NULL)
    ##     shinyjs::reset(ns("login_username"))
    ##     shinyjs::reset(ns("login_password"))
    ##     USER$logged <- FALSE
    ## })

    res <- list(
        name   = shiny::reactive(USER$name),
        logged = shiny::reactive(USER$logged)
    )
    return(res)
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
    m <- particlesSplashModal(body, ns=ns, easyClose=TRUE, fade=TRUE,
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
    titles[[4]] = c("Do-it-yourself Analytics","Yes you can!")
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
    titles[[18]] = c("Ich bin doch nicht blöd!","Of course I use Omics Playground")
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
        div.firebase <- firebase::useFirebaseUI()
    }

    div.alt <- div()
    if(!is.null(alt)) div.alt <- alt
    top <- HTML(rep("<br>",sum(c(!with.email,!with.username,!with.password))))

    div.button <- div(
        id="splash-buttons",
        actionButton(ns("login_btn"),login.text,class="red-button")
    )
    if(with.register) {
        div.button <- div(
            id="splash-buttons",
            actionButton(ns("login_btn"),login.text,class="red-button"),
            actionButton(ns("register_btn"),"Register",class="red-button")
        )
    }
    
    ##splash.panel=div();ns=function(x)x
    if(with.firebase) {
        splash.panel <- div(div.firebase)
    } else {
        splash.panel <- div(
            id="splash-panel",            
            br(),br(),top,
            div.username,
            div.email,
            div.password,
            div.alt,
            br(),
            div.button
        )
    }

    body <- tagList(
        div(id="splash-title",splash.title),
        splash.panel
    )

    m <- particlesSplashModal(body, ns=ns)

    return(m)
}

particlesSplashModal <- function(body, ns=NULL, easyClose=FALSE, fade=FALSE,
                                 buttons=TRUE, footer=TRUE)
{
    
    if(is.null(ns)) ns <- function(e) return(e)
    message("[AuthenticationModule::splashModal]")

    div.footer = shiny::modalButton("Dismiss")
    if(buttons) {
        div.footer = shiny::tagList(
            shiny::actionButton(ns("action1"),"Read-the-docs", icon=icon("book"),
                         onclick="window.open('https://omicsplayground.readthedocs.io','_blank')"),
            shiny::actionButton(ns("action2"),"Watch tutorials", icon=icon("youtube"),
                         onclick="window.open('https://www.youtube.com/channel/UChGASaLbr63pxmDOeXTQu_A','_blank')"),
            shiny::actionButton(ns("action3"),"Get the source", icon=icon("github"),
                         onclick="window.open('https://github.com/bigomics/omicsplayground','_blank')"),
            shiny::actionButton(ns("action4"),"Docker image", icon=icon("docker"),
                         onclick="window.open('https://hub.docker.com/r/bigomics/omicsplayground','_blank')"),
            shiny::actionButton(ns("action5"),"User forum", icon=icon("users"),
                         onclick="window.open('https://groups.google.com/d/forum/omicsplayground','_blank')"),
            ##actionButton(ns("action_beer"),"Buy us a beer!", icon=icon("beer"),
            ##             onclick="window.open('https://www.buymeacoffee.com/bigomics','_blank
            shiny::actionButton(ns("action_beer"),"Buy us a coffee!", icon=icon("coffee"),
                         onclick="window.open('https://www.buymeacoffee.com/bigomics','_blank')")
            ## shiny::modalButton("Let's start!")
            ## shiny::actionButton(ns("action_play"), "Let's play!", class="red-button")
        )
    }
    if(!footer) {
        div.footer <- NULL
    }

    ## return modalDialog
    particlesjs.conf <- rjson::fromJSON(file="resources/particlesjs-config.json")
    m <- shiny::modalDialog(
        id = "modal-splash",
        shiny::div(
            id="particles-target",
            ##img(src = base64enc::dataURI(file="www/splash.png"),
            ##    width="100%", height="auto%", style="position:absolute;"),
            shiny::div(id="splash-logo", shiny::img(src=base64enc::dataURI(file="www/logo.png"),
                                      width=32,height=32)),
            body,
            ##firebase::useFirebaseUI(),
            shiny::br(),
            shiny::div(id="splash-warning",textOutput(ns("login_warning")),style="color:red;"),
            style="height: 500px; width: 100%;"                
        ),
        footer = div.footer,
        particlesjs::particles(config=particlesjs.conf, target_id = "particles-target", timeout = 1000),
        size="m", easyClose=easyClose, fade=fade
    ) ## end of modalDialog
    
    return(m)
}
