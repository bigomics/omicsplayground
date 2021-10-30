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


NoAuthenticationModule <- function(input, output, session, username="", email="")
{
    message("[AuthenticationModule] >>>> using no authentication <<<<")
    ns <- session$ns    
    USER <- shiny::reactiveValues(
                       logged=FALSE,
                       name="",
                       email="",
                       level="",
                       limit="")    

    resetUSER <- function() {
        USER$logged <- FALSE
        USER$name <- ""
        USER$email <- ""
        USER$level <- ""
        USER$limit <- ""
    }
    
    output$showLogin <- shiny::renderUI({
        resetUSER()
        m <- splashLoginModal(
            ns=ns, with.email=FALSE, with.password=FALSE, login.text="Start")
        ## shinyjs::delay(2000, {output$login_warning <- shiny::renderText("")})
        shiny::showModal(m)

        USER$name   <- username
        USER$email  <- email
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
        m <- splashLoginModal(
            ns=ns, with.email=FALSE, with.password=FALSE, login.text="Start")
        ## shinyjs::delay(2000, {output$login_warning <- shiny::renderText("")})
        shiny::showModal(m)
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


upgrade.dialog = shiny::HTML("
Do you want to remove the 60 minutes time limit? Do you want to be able to save more datasets?
<br><br><center><table width=90% style='background-color:#F4FAFF;'><tr>
<th>BASIC<br></th>
<th>STARTER<br></th>
<th>PREMIUM</th>
<th>ENTERPRISE</th></tr>
<tr><td>Try out for free</td>
<td>Great to start</td>
<td>For professionals</td>
<td>Enterprise-Ready</td>
<tr><td><h3><b>FREE</b></h3></td>
<td><h3><b>Soon!</b></h3></td>
<td><h3><b>Soon!</b></h3></td>
<td><h3><b>Contact us!</b></h3></td>
</tr><tr><td>&nbsp;</tr>
<tr>
<td><ul><li>Host up to 3 datasets</li><li>60 minutes time limit</li>
<li>Up to 25 samples/dataset</li><li>Up to 5 comparisons</li></ul>
<td><ul><li>Host up to 10 datasets</li><li>3 hours time limit</li>
<li>Up to 100 samples/dataset</li><li>Up to 10 comparisons</li></ul>
<td><ul><li>Host up to 100 datasets</li><li>8 hours time limit</li>
<li>Up to 2000 samples/dataset</li><li>Up to 100 comparisons</li></ul>
<td><ul><li>Host unlimited datasets</li><li>No time limit</li>
<li>Up to 2000 samples/dataset</li><li>Up to 100 comparisons</li></ul>

<tr>
<td>
<td><a onClick='upgrade_plan()' style='font-weight:bold;color:#2a9d8f;cursor:pointer;' id='authentication-upgrade'>Get Starter!</a>
<td><a onClick='upgrade_plan()' style='font-weight:bold;color:#2a9d8f;cursor:pointer;' id='authentication-upgrade'>Get Premium!</a>
<td><a style='font-weight:bold;color:#2a9d8f' href='mailto:info@bigomics.ch'>Send Email</a>
</table></center><br><br>
")

FirebaseAuthenticationModule <- function(input, output, session)
{
    message("[AuthenticationModule] >>>> using FireBase (email+password) authentication <<<<")

    firebase_config <- firebase:::read_config("firebase.rds")
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

    firebase <- firebase::FirebaseUI$
        new(persistence = "local")$ # instantiate
        set_providers( # define providers
            email_link = TRUE, 
            google = TRUE
        )$
        set_privacy_policy_url(
            "https://bigomics.ch/privacy/"
        )$
        set_tos_url(
            "https://bigomics.ch/terms/"
        )
    firebase$set_tos_url("https://bigomics.ch/terms")
    firebase$set_privacy_policy_url("https://bigomics.ch/privacy")    
    
    resetUSER <- function() {
        USER$logged <- FALSE
        USER$name <- ""
        USER$password <- ""
        USER$email <- ""
        USER$level <- ""
        USER$limit <- ""
        USER$token <- ""
    }

    first_time = TRUE
    
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

        # no need to show the modal
        # if the user is logged
        # this is due to persistence
        if(USER$logged && first_time) {
            dbg("[FirebaseAuthenticationModule] USER is already logged in & first time = TRUE")
            dbg("[FirebaseAuthenticationModule] signing out any previous user")
            firebase$sign_out()
            resetUSER()
            first_time <<- FALSE            
            ##return()
        } else if(USER$logged) {
            dbg("[FirebaseAuthenticationModule] USER is already logged in! no modal")                        
            return()
        }
        
        on.exit({
            dbg("[FirebaseAuthenticationModule] on.exit")            
            firebase$launch()
        })

        dbg("[FirebaseAuthenticationModule] showing Firebase login modal")                                
        shiny::tagList(
            shiny::showModal(m)
        )
    })

    observeEvent( input$firebaseLogout, {    

        dbg("[FirebaseAuthenticationModule] observe::input$firebaseLogout reacted")        
        
        on.exit({
            firebase$launch()
        })

        dbg("[FirebaseAuthenticationModule] signing out from Firebase")
        firebase$sign_out()

        dbg("[FirebaseAuthenticationModule] reset user")        
        resetUSER()
        
        m <- splashLoginModal(
            ns = ns,
            with.email = FALSE,
            with.username = FALSE,
            with.password = FALSE,
            with.register = FALSE,
            with.firebase = TRUE,            
            login.text = "Start!"
        )
        
        shiny::showModal(m)
    })


    observeEvent(firebase$get_signed_in(), {

        dbg("[FirebaseAuthenticationModule] observe::get_signed_in() reacted")

        response <- firebase$get_signed_in()

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
        USER$uid <- as.character(response$response$uid)
        USER$name  <- response$response$displayName
        USER$email <- response$response$email

        dbg("[FirebaseAuthenticationModule@firebase$get_signed_in] is.null(user.name) = ",is.null(USER$name) )
        dbg("[FirebaseAuthenticationModule@firebase$get_signed_in] is.null(user.email) = ",is.null(USER$email) )
        if(!is.null(USER$name))  USER$name  <- as.character(USER$name)
        if(!is.null(USER$email)) USER$email <- as.character(USER$email)

        if(is.null(USER$name))  USER$name  <- ""
        if(is.null(USER$email)) USER$email <- ""

        dbg("[FirebaseAuthenticationModule@firebase$get_signed_in] user.name==''  = ",USER$name=='' )
        dbg("[FirebaseAuthenticationModule@firebase$get_signed_in] user.email=='' = ",USER$email=='' )
        
        session$sendCustomMessage(
            "get-permissions",
            list(
                ns = ns(NULL)
            )
        )
    })

    observeEvent(input$stripeId, {
        USER$stripe_id <- input$stripeId$id
        USER$href <- input$stripeId$href
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
        shinyalert::shinyalert(
                        title = "Coming Soon!",
                        text = upgrade.dialog,
                        html=TRUE,
                        animation = FALSE,
                        size = 'l',
                        immediate = TRUE
                    )
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
        m <- splashLoginModal(
            ns=ns,
            with.email=TRUE,
            with.username=FALSE,
            with.password=TRUE)
        ## shinyjs::delay(2000, {output$login_warning <- shiny::renderText("")})
        shiny::showModal(m)
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
## RegisterAuthenticationModule (just ask email, no password)
##================================================================================

register.file="../logs/register.log"
RegisterAuthenticationModule <- function(input, output, session, register.file)
{
    message("[AuthenticationModule] >>>> using Register authentication <<<<")

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
    USER <- shiny::reactiveValues(
                       logged=FALSE,
                       name="",
                       email="",
                       password=NA,
                       registered=FALSE,
                       level="")
        

    resetUSER <- function() {
        USER$logged <- FALSE
        USER$name <- ""
        USER$email <- ""
        USER$password <- NA
        USER$registered <- FALSE
        USER$level <- ""
    }

    output$showLogin <- shiny::renderUI({
        message("[AuthenticationModule::UI] USER$logged = ",USER$logged)
        ## If login is OK then do nothing
        if(USER$logged) {
            return(NULL)
        } else {
            showLoginDialog()             
        }
    })

    observeEvent( input$firebaseLogout, {
        dbg("[RegisterAuthenticationModule] observe::input$firebaseLogout() reacted")
        resetUSER()
        REGISTERED <<- read.csv(register.file, colClasses="character",
                               stringsAsFactors=FALSE)
        showLoginDialog()                     
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
        email.exists <- (email %in% REGISTERED$email)
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
        session$sendCustomMessage("set-user", list(user = USER$email))
        
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
            USER$level <- "FREE"
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
        email  = shiny::reactive(USER$email),
        level  = shiny::reactive(USER$level),
        logged = shiny::reactive(USER$logged),
        limit = shiny::reactive(USER$limit)
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
    top <- shiny::HTML(rep("<br>",sum(c(!with.email,!with.username,!with.password))))

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
