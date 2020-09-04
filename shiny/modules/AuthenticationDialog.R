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
    ns <- NS(id)  ## namespace
    showModal(uiOutput(ns("showLogin")))
}


NoAuthenticationDialog <- function(input, output, session, username=NULL)
{
    ns <- session$ns    
    USER <- reactiveValues(logged=FALSE, name="", email="",
                           password=NA, registered=NA, level="")    
    USER$name <- username
    
    output$showLogin <- renderUI({
        m <- splashLoginModal(ns=ns, with.email=FALSE, with.password=FALSE)
        ## shinyjs::delay(2000, {output$login_warning <- renderText("")})
        showModal(m)
    })

    output$login_warning = renderText("")

    observeEvent( input$login_btn, {           
        removeModal()
        showModal(splashHelloModal(name=USER$name,ns=ns))
        USER$logged <- TRUE
    })

    observeEvent(input$hello_letsplay,{
        removeModal()
    })
    
    rt <- list(
        level  = reactive(USER$level),
        name = reactive(USER$name),
        logged = reactive(USER$logged)
    )
    return(rt)
}

PasswordAuthenticationDialog <- function(input, output, session,
                                         credentials.file)
{
    ns <- session$ns    
    USER <- reactiveValues(logged=FALSE, username=NA, password=NA, level=NA)    
    credentials <- read.csv(credentials.file,colClasses="character")
    head(credentials)

    output$showLogin <- renderUI({
        m <- splashLoginModal(
            ns=ns,
            with.email=FALSE,
            with.username=TRUE,
            with.password=TRUE)
        ## shinyjs::delay(2000, {output$login_warning <- renderText("")})
        showModal(m)
    })

    output$login_warning = renderText("")

    observeEvent( input$login_btn, {           

        cat("name=",input$login_username,"\n")
        cat("email=",input$login_email,"\n")
        cat("password=",input$login_password,"\n")
        cat("Logged=",USER$logged,"\n")

        login.OK   = FALSE
        valid.date = FALSE
        valid.user = FALSE
        
        if( is.null(input$login_username) || is.null(input$login_password)) return(NULL)
        if( input$login_username=="" || input$login_password=="") return(NULL)    
        username <- input$login_username
        sel <- tail(which( credentials$username == username),1)
        valid.user <- isTRUE(username %in% credentials$username)
        valid.pw   <- isTRUE(credentials[sel,"password"]==input$login_password)
        valid.date <- isTRUE( Sys.Date() < as.Date(credentials[sel,"expiry"]) )
        login.OK = (valid.user && valid.pw && valid.date)
        message("--------- password login ---------")
        message("input.username=",input$login_username)
        message("input.password=",input$login_password)
        message("user.password=",credentials[sel,"password"])
        message("user.expiry=",credentials[sel,"expiry"])
        message("user.name=",credentials[sel,"name"])
        message("valid.user=",valid.user)
        message("valid.date=",valid.date)
        message("valid.pw=",valid.pw)
        message("----------------------------------")
        
        if (login.OK) {
            output$login_warning = renderText("")
            removeModal()
            ##USER$name   <- input$login_username
            ##USER$email <- credentials[sel,"email"]
            USER$name  <- username
            ## Here you can perform some user-specific functions, or site news
            showModal(splashHelloModal(USER$name,ns=ns))
            ##removeModal()
            USER$logged <- TRUE            
        } else {
            if(!valid.date) {
                output$login_warning = renderText("Registration expired")
            } else {
                output$login_warning = renderText("Invalid username or password")
            }
            ##shinyjs::delay(2000, hide("login_warning", anim = TRUE, animType = "fade"))
            shinyjs::delay(2000, {output$login_warning <- renderText("")})
            USER$logged <- FALSE
        }
        ##hide("login_warning")
    })

    observeEvent(input$hello_letsplay,{
        removeModal()
    })

    rt <- list(
        level  = reactive(USER$level),
        name   = reactive(USER$username),
        logged = reactive(USER$logged)
    )
    return(rt)
}

register.file="../logs/register.log"
RegisterAuthenticationDialog <- function(input, output, session, register.file)
{
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
    USER <- reactiveValues(logged=FALSE, name="", email="",
                           password=NA, registered=FALSE, level="")
        
    output$showLogin <- renderUI({
        message("[AuthenticationDialog::UI] USER$logged = ",USER$logged)
        ## If login is OK then do nothing
        if(USER$logged) {
            return(NULL)
        } else {
            showLoginDialog()             
        }
    })

    showLoginDialog <- function() {
        ##removeModal()
        alt <- actionLink(ns("create_account"),"or create an account",
                          style="color:white;", class="white-link")
        showModal(splashLoginModal(
            ns=ns, with.password=FALSE, with.username=FALSE,
            alt = alt))
    }
    
    observeEvent( input$login_btn, {           
        
        message("[AuthenticationDialog::login] login_btn pressed")                
        message("[AuthenticationDialog::login] email = ",input$login_email)
        email <- input$login_email
        ## check login
        email.exists <- email %in% REGISTERED$email
        if(!email.exists) {
            output$login_warning = renderText("invalid email")
            shinyjs::delay(2000, {output$login_warning <- renderText("")}) 
            return(NULL)
        }

        sel <- tail(which(REGISTERED$email == email),1)
        user.registered <- REGISTERED[sel,"registered"]
        USER$name <- REGISTERED[sel,"name"]
        USER$email <- email
        USER$registered <- user.registered
        USER$logged <- TRUE
        
        message("[AuthenticationDialog::login] email.exists = ",email.exists)
        message("[AuthenticationDialog::login] user.registered = ",user.registered)
        ##showDialog("welcome")
        removeModal()
        Sys.sleep(DIALOG_SLEEP)
        showModal(splashHelloModal(
            ns = ns,
            name = USER$name,
            msg = "Welcome back. We wish you many great discoveries today!"
        ))
        observeEvent(input$hello_letsplay,{
            removeModal()
        })
        
    })
        
    observeEvent( input$create_account, {
        removeModal()
        Sys.sleep(DIALOG_SLEEP)
        m <- createAccountModal(USER)
        showModal(m)
    })
    
    validEmailFormat <- function(email) {
        has.atsign <- grep("@",email)
        domain <- strsplit(email,split='@')[[1]][2]
        has.dot <- length(strsplit(domain,split="[.]")[[1]]) > 1
        has.atsign && has.dot
    }

    createAccountModal <- function(user) {
        message("[AuthenticationDialog::createAccount] user$name = ",user$name)            
        m <- modalDialog(
            id = "auth_dialog",
            ##title = "Create an Account",
            tagList(
                fillRow(
                    flex = c(0.12,0.5,0.2,0.5,0.15),
                    height = 400,
                    br(),
                    tagList(
                        h3("Try premium?"),
                        HTML("Good news! For limited time, all registered users will receive free access to our <b>Premium plan</b> which includes:<br><br>"),
                        HTML("<ul><li>Drug enrichment<li>Biomarker analysis<li>Connectivity mapping<li>And more!</ul>"),
                        br(),
                        ##h4("Recommend us to your friends"),
                        ##textAreaInput(ns("register_friends"),"Email(s) of friends:",rows=3)
                        ),
                    br(),
                    tagList(
                        ##h3("Tell us a bit about yourself"),
                        h3("Create a free account"),                            
                        textInput(ns("register_name"),"First name"),
                        textInput(ns("register_email"),"E-mail",value=USER$email),
                        selectInput(ns("register_company"),"What kind of company is it?",
                                    ## select=FALSE, selectize=FALSE,
                                    c("",sort(c("Academia","Hospital",
                                                "Non-profit organization",
                                                "Other industry (SME)",
                                                "Other industry (large)",
                                                "Pharmaceutical industry","Start-up",
                                                "Self-employed")),"other")),
                        selectInput(ns("register_role"),"What is your role there?",
                                    c("","Biologist","Bioinformatician","Manager","Other")),
                        selectInput(ns("register_hear"),"How did you hear about us?",
                                    c("",sort(c("Coworker/friend","Event",
                                                "Contacted by Sales",
                                                "Article","Blog","Web search",
                                                "LinkedIn","Twitter",
                                                "Newsleter")),"other")),
                        br(),br()
                    ),
                    br()
                )
            ),
            footer = tagList(
                div(textOutput(ns("register_warning")),style="color: red;"),
                actionLink(ns("register_cancel"), "Cancel"),HTML("&nbsp;&nbsp;"),
                actionButton(ns("register_submit"), "Register")
            ),                    
            ##footer = NULL,
            easyClose = FALSE,
            size = "l"
        )

        observeEvent( input$register_cancel, {
            removeModal()	
            Sys.sleep(DIALOG_SLEEP)
            showLoginDialog() 
        }, once=TRUE, ignoreInit=TRUE)

        observeEvent( input$register_submit, {
            
            message("[AuthenticationDialog::createAccount] REGISTER pressed")
            qq <- c(name = input$register_name,
                    company = input$register_company,
                    role = input$register_role,
                    hear = input$register_hear)

            if(any(is.null(qq) | qq=="")) {
                message("[AuthenticationDialog::createAccount] missing fields")
                output$register_warning = renderText("Please fill in all required fields")
                shinyjs::delay(2000, {output$register_warning <- renderText("")})
                return(NULL)
            }

            email.ok = validEmailFormat( input$register_email )
            if(!email.ok) {
                output$register_warning = renderText("Error: invalid email format")
                shinyjs::delay(2000, {output$register_warning <- renderText("")}) 
                return(NULL)
            }

            email.exists <- (input$register_email %in% REGISTERED$email)
            if(email.exists) {
                output$register_warning = renderText("Error: email exists")
                shinyjs::delay(2000, {output$register_warning <- renderText("")}) 
                return(NULL)
            }

            message("[AuthenticationDialog::createAccount] all answered!")
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
            
            answers <- c(
                register_company = input$register_company,
                register_role = input$register_role,
                register_hear = input$register_hear
            )
            updateSurveyLog(USER, answers)                
            USER$logged = TRUE

            ## success
            removeModal()
            showModal(splashHelloModal(
                ns = ns,                
                name = USER$name,
                msg = "Thank you for registering. We wish you many great discoveries today!"
            ))

            observeEvent(input$hello_letsplay,{
                removeModal()
            })
            
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
    output$login_warning = renderText("")
    output$register_warning = renderText("")
    res <- list(
        level  = reactive(USER$level),
        name   = reactive(USER$name),
        logged = reactive(USER$logged)
    )
    return(res)
}

credentials.file="./CREDENTIALS"
FullAuthenticationDialog.SAVE <- function(input, output, session,
                                          credentials.file,
                                          type="splash")
{
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

    USER <- reactiveValues(logged=FALSE, name="", email="",
                           password=NA, verified=FALSE, level="")
    
    SMTP_SERVER  <- Sys.getenv("SMTP_SERVER")
    SMTP_USER    <- Sys.getenv("SMTP_USER")
    SMTP_SECRET  <- Sys.getenv("SMTP_SECRET")
    has.smtp <- SMTP_SERVER!="" && SMTP_USER!="" && SMTP_SECRET!=""
    has.smtp
    message("[AuthenticationDialog] has.smtp = ",has.smtp)
    
    DIALOG_SLEEP = 0.3
    ##waiter_show(spin_fading_circles())
    waiter <- Waiter$new(id="animated_captcha", color=transparent(0.5))

    validEmailFormat <- function(email) {
        has.atsign <- grep("@",email)
        domain <- strsplit(email,split='@')[[1]][2]
        has.dot <- length(strsplit(domain,split="[.]")[[1]]) > 1
        has.atsign && has.dot
    }
    
    output$showLogin <- renderUI({
        message("[AuthenticationDialog::UI] USER$logged = ",USER$logged)
        ## If login is OK then do nothing
        if(USER$logged) {
            return(NULL)
        } else {
            showDialog("login")
        }
    })
        
    observeEvent( input$login_btn, {           
        
        message("[AuthenticationDialog::login] login_btn pressed")                
        message("[AuthenticationDialog::login] email = ",input$login_email)
        email <- input$login_email

        ## check login
        valid.email <- validEmailFormat(email)
        if(!valid.email) {
            message("[AuthenticationDialog::login] email not valid")
            output$login_warning = renderText("invalid email format")
            shinyjs::delay(2000, {output$login_warning <- renderText("")}) 
            return(NULL)
        }
        
        email.exists <- email %in% CREDENTIALS$email
        user.registered <- FALSE
        nlogin <- sum(CREDENTIALS$email == email)
        ##USER <- reactiveValues(logged=FALSE, name=NA, email=NA, verified=FALSE)
        if(email.exists) {
            sel <- tail(which(CREDENTIALS$email == email),1)
            user.registered <- CREDENTIALS[sel,"registered"]
            USER$name <- CREDENTIALS[sel,"name"]
        }
        USER$email <- email
        USER$registered <- user.registered
        USER$logged <- FALSE
        
        message("[AuthenticationDialog::login] email.exists = ",email.exists)
        message("[AuthenticationDialog::login] user.registered = ",user.registered)
        message("[AuthenticationDialog::login] nlogin = ",nlogin)
        
        if(email.exists && user.registered) {
            ## existing user, already registered
            for(i in 1:10) removeModal()
            Sys.sleep(DIALOG_SLEEP)
            updateCredentials(USER, credentials.file)
            if(FALSE && nlogin == 2) {
                message("[AuthenticationDialog::login] >>> asking user experience")
                showDialog("ask.experience")
            } else {
                message("[AuthenticationDialog::login] >>> welcome  user")                
                showDialog("welcome")
            }
            USER$logged <- TRUE
        } else  {
            ## fist time or unregistered user
            message("[AuthenticationDialog::login] >>> creating account")
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
    observeEvent( input$verify_btn, {
        message("[AuthenticationDialog::register] verify_btn pressed")
        message("[AuthenticationDialog::register] USER$email = ",USER$email)
        message("[AuthenticationDialog::register] USER$password = ",USER$password)        
        message("[AuthenticationDialog::register] verify_code = ",input$verify_code)

        if(input$verify_code == "") return(NULL)

        ok = (USER$password == input$verify_code)
        message("[AuthenticationDialog::register] ok = ",ok)        
        
        if(ok) {
            output$verify_warning = renderText("")
            output$captcha_warning = renderText("")
            ##removeModal()	
            for(i in 1:10) removeModal()
            Sys.sleep(DIALOG_SLEEP)
            updateCredentials(USER, credentials.file)
            showDialog("welcome")
            USER$logged <- TRUE
        } else {
            output$verify_warning = renderText("Invalid verification")
            output$captcha_warning = renderText("Invalid verification")
            shinyjs::delay(2000, {output$verify_warning <- renderText("")})
            shinyjs::delay(2000, {output$captcha_warning <- renderText("")})            
        }        
    })

    observeEvent( input$cancel_btn, {
        removeModal()	
        Sys.sleep(DIALOG_SLEEP)
        showDialog("login")
    })
    
    ##---------------------------------------------------
    ## functions
    ##---------------------------------------------------
    
    showDialog <- function(stage, data=NULL)
    {
        message("[AuthenticationDialog::showModalDialog] stage = ",stage)
        ##removeModal()

        output$verify_warning = renderText("")
        output$login_warning = renderText("")
        output$register_warning = renderText("")

        m <- NULL
        if(stage=="login" && type=="simple") {
            m <- modalDialog(
                id = "auth_dialog",
                title = "Log in to your Playground",
                tagList( 
                    p( "Enter your email to log in"),
                    textInput(ns("login_email"), "E-mail:"),
                    div(actionButton(ns("login_btn"), "Login"),style="text-align: center;")
                ),
                footer = NULL,
                size = "s"
            )
        } else if(stage=="login" && type=="splash") {
            m <- splashLoginModal(ns=ns, with.password=FALSE)            
        } else if(stage=="verify") {
            m <- modalDialog(
                id = "auth_dialog",
                title = "Verify your email",
                tagList( 
                    p("Please enter the verification code that we sent to your email:"),
                    textInput(ns("verify_code"), "Verification code:"),
                    div( actionButton(ns("cancel_btn"), "Cancel"),
                        actionButton(ns("verify_btn"), "Submit"),
                        style="text-align: center;")
                ),
                footer = div(textOutput(ns("verify_warning")),style="color: red;"),
                size = "s"
            )
        } else if(stage=="captcha") {            
            m <- modalDialog(
                id = "auth_dialog",
                title = "Please verify you're human",
                tagList( 
                    p("Copy the amino acid letters into the captcha box below:"),
                    ##div(img(src=base64enc::dataURI(file=tmpfile), width="240px"), style="text-align: center;"),
                    imageOutput(ns("animated_captcha"),height=50,width=240),
                    br(),
                    textInput(ns("verify_code"), "Enter captcha:"),                    
                    div(
                        actionButton(ns("verify_recaptcha"), "Refresh captcha"), 
                        actionButton(ns("verify_btn"), "Submit"), style="text-align: center;")
                ),
                footer = div(textOutput(ns("captcha_warning")),style="color: red;"),                
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
            m <- modalDialog(
                id = "auth_dialog",
                title = "Error",
                tagList( p("error.") ),
                footer = div(HTML("error"),style="color: red;"),
                size = "s"
            )
        }
        showModal(m)
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
            plot(0,0,pch="",xlim=c(-0.1,0.95),ylim=c(0,1),bty="n",axes=FALSE)
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
    
    output$animated_captcha <- renderImage({
        input$verify_recaptcha
        require(animation)
        tmpfile="./captcha.png"
        tmpfile <- tempfile()        
        message("[AuthenticationDialog:animated_captcha] render GIF")        
        saveGIF( animateCaptcha(),
                movie.name = tmpfile, interval = 0.01, fps=70,
                nmax = 100, ani.width = 240, ani.height = 50, loop=1)
        ## Return a list containing information about the image
        message("[AuthenticationDialog:animated_captcha] saved GIF")
        list(
            src = tmpfile,
            contentType = "image/gif",
            width = 240,
            height = 50,
            alt = "animated captcha"
        )
    })
    
    createAccountModal <- function(user) {
        message("[AuthenticationDialog::createAccount] user$name = ",user$name)            
        m <- modalDialog(
            id = "auth_dialog",
            ##title = "Tell us a bit about yourself",
            tagList(
                fillRow(
                    flex = c(0.02,0.4,0.03,0.6),
                    height = 380,
                    br(),
                    tagList(
                        h3("Try premium?"),
                        HTML("Good news! For limited time, all registered users will receive free trial of our <b>Premium plan</b> which includes:<br><br>"),
                        HTML("<ul><li>Drug enrichment<li>Biomarker analysis<li>Connectivity mapping<li>And more!</ul>"),
                        br(),
                        ##h4("Recommend us to your friends"),
                        ##textAreaInput(ns("ask_friends"),"Email(s) of friends:",rows=3)
                        ),
                    br(),
                    tagList(
                        h3("Tell us a bit about yourself"),                            
                        textInput(ns("ask_name"),"First name"),
                        textInput(ns("ask_email"),"E-mail",value=USER$email),
                        selectInput(ns("ask_company"),"What kind of company is it?",
                                    ## select=FALSE, selectize=FALSE,
                                    c("",sort(c("Academia","Hospital","Non-profit organization",
                                      "Other industry (SME)","Other industry (large)",
                                      "Pharmaceutical industry","Start-up",
                                      "Self-employed")),"other")),
                        selectInput(ns("ask_role"),"What is your role there?",
                                    c("","Biologist","Bioinformatician","Manager","Other")),
                        selectInput(ns("ask_hear"),"How did you hear about us?",
                                    c("",sort(c("From coworker/friend","LinkedIn","Twitter",
                                      "Journal article","Web article/blog",
                                      "Web search","BigOmics website", "Event",
                                      "Contacted by Sales")), "other"))
                    ),
                    br(),br()
                )
            ),
            footer = tagList(
                div(textOutput(ns("ask_warning")),style="color: red;"),
                actionLink(ns("ask_btn_skip"), "Skip"),HTML("&nbsp;&nbsp;"),
                actionButton(ns("ask_register"), "Register")
            ),                    
            ##footer = NULL,
            easyClose = FALSE,
            size = "m"
        )
        observeEvent( input$ask_btn_skip, {
            removeModal()
        })
        observeEvent( input$ask_register, {
            message("[AuthenticationDialog::createAccount] REGISTER pressed")
            qq <- c(name = input$ask_name,
                    company = input$ask_company,
                    role = input$ask_role,
                    hear = input$ask_hear)

            if(any(is.null(qq) | qq=="")) {
                message("[AuthenticationDialog::createAccount] missing fields")
                output$ask_warning = renderText("Please fill in all required fields")
                shinyjs::delay(2000, {output$ask_warning <- renderText("")})                
                return(NULL)
            } else {
                message("[AuthenticationDialog::createAccount] all answered!")
                lapply(1:length(qq), function(i) message(names(qq)[i],"=",qq[i]))
                ##
                ## save somewhere!
                ##
                USER$name <- input$ask_name
                updateCredentials(USER, credentials.file)                

                answers <- c(
                    ask_company = input$ask_company,
                    ask_role = input$ask_role,
                    ask_hear = input$ask_head
                )
                updateSurveyLog(USER, answers)                
            }
            USER$logged = TRUE
        })
        m
    }   
    
    askExperienceModal <- function(user) {
            message("[AuthenticationDialog::askExperience]")            
            m <- modalDialog(
                id = "auth_dialog",
                ##title = "Tell us a bit about yourself",
                tagList(
                    fillRow(
                        flex = c(0.02,0.6,0.04,0.33),
                        height = 380,
                        br(),
                        tagList(
                            h3("Tell us how you like us"),br(),
                            radioButtons(
                                ns("askxp_satisfied"),
                                "Are you so far satisfied with your Playground experience?",
                                choiceValue=c(1,2,3,4,5), inline=TRUE,
                                choiceNames=c("no","so-so","ok","satisfied","love it!")),
                            radioButtons(
                                ns("askxp_recommend"),"Would you recommend us to a friend",
                                choiceValue=c(1,2,3),
                                choiceNames=c("not really","not sure","definitely!"),
                                inline=TRUE),
                            radioButtons(
                                ns("askxp_useagain"),"Would you use us again?",
                                c("no","depends","yes"), inline=TRUE),
                            selectInput(
                                ns("askxp_problems"),"What was your biggest problem?",
                                c("Data preparation","Slow website","Server disconnects",
                                  "Too complex", "No problems","Other"))
                            ##selectInput(ns("ask_cosize"),"How big is the company?",
                            ##            c("micro (<10)","small (<50)","medium (<250)",
                            ##              "large (> 250")),
                        ), br(),
                        tagList(
                            h3("Try premium?"),br(),
                            HTML("This is your chance! Complete this survey and unlock a free trial of our <b>Premium plan</b> which includes:<br><br>"),
                            HTML("<ul><li>Drug enrichment<li>Biomarker analysis<li>Connectivity mapping<li>And more!</ul>"),
                            br(),
                            ##h4("Recommend us to your friends"),
                            ##textAreaInput(ns("ask_friends"),"Email(s) of friends:",rows=3)
                        )
                    ),
                    div(
                        actionLink(ns("askxp_btn_skip"), "Skip"),HTML("&nbsp;&nbsp;"),
                        actionButton(ns("askxp_submit"), "Submit"),
                        style="text-align: right;"
                    )                    
                ),
                footer = NULL,
                easyClose = FALSE,
                size = "m"
            )

            observeEvent( input$askxp_btn_skip, {
                removeModal()
            })
            observeEvent( input$askxp_submit, {
                message("[AuthenticationDialog::askExperience] Submit pressed")
                answers <- c(
                    askxp_satisfied = input$askxp_satisfied,
                    askxp_recommend = input$askxp_recommend,
                    askxp_useagain = input$askxp_useagain,
                    askxp_problems = input$askxp_problems
                )                    
                updateSurveyLog(USER, answers)
                removeModal()
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
        require(emayili)
        require(magrittr) 
        
        message("[AuthenticationDialog::smtpAuthCode] *************************")
        message("[AuthenticationDialog::smtpAuthCode] ****** SENDING EMAIL ****")
        message("[AuthenticationDialog::smtpAuthCode] *************************")
        message("[AuthenticationDialog::smtpAuthCode] email = ",email)
        message("[AuthenticationDialog::smtpAuthCode] code = ",auth.code)

        SMTP_SERVER  <- Sys.getenv("SMTP_SERVER")
        SMTP_USER    <- Sys.getenv("SMTP_USER")
        SMTP_SECRET  <- Sys.getenv("SMTP_SECRET")        
        message("[AuthenticationDialog::smtpAuthCode] SMTP_SERVER = ",SMTP_SERVER)
        message("[AuthenticationDialog::smtpAuthCode] SMTP_USER = ",SMTP_USER)
        message("[AuthenticationDialog::smtpAuthCode] SMTP_SECRET = ",SMTP_SECRET)        

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
    output$login_warning = renderText("")

    ## observeEvent( input$logout, {
    ##     ##updateTextInput(session, ".username", value=NULL)
    ##     reset(ns("login_username"))
    ##     reset(ns("login_password"))
    ##     USER$logged <- FALSE
    ## })

    res <- list(
        name   = reactive(USER$name),
        logged = reactive(USER$logged)
    )
    return(res)
}

splashHelloModal <- function(name, msg=NULL, ns=NULL, duration=-4000)
{
    if(is.null(ns)) ns <- function(e) return(e)
    message("[AuthenticationDialog::splashHelloModel]")
    
    all.hello = c("Hello","Salut","Hola","Privet","Nǐ hǎo","Ciao","Hi","Hoi","Hej",
                  "Yassou","Selam","Hey","Hei","Grützi","Салам")
    title = paste( paste0(sample(all.hello,3),"!"), collapse=" ")
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
    splash.title <- div(
        br(),br(),br(),br(),
        div(HTML(title),style="font-size:70px;font-weight:700;line-height:1em;"),
        br(),
        div(HTML(subtitle),style="font-size:30px;"),
        br(),br(),br(),
        actionButton(ns("hello_letsplay"),"Let's play!",class="red-button")
        ## actionButton("hello_letsplay","Let's play!",class="red-button")        
    )
    body <- tagList(
        div(id="splash-title",splash.title)
    )
    m <- particlesSplashModal(body, ns=ns, easyClose=TRUE, fade=TRUE,
                              buttons=FALSE, footer=FALSE)    
    if(duration>0) shinyjs::delay(duration, removeModal())
    return(m)
}

splashLoginModal <- function(ns=NULL, with.email=TRUE, with.password=TRUE,
                             with.username=FALSE, alt=NULL)
{
    if(is.null(ns)) ns <- function(e) return(e)
    message("[AuthenticationDialog::splashLoginModal]")

    titles <- list()
    titles[[1]] = c("Big Omics Data","Isn't big anymore with Omics Playground")
    titles[[2]] = c("Great Discoveries","Start on the Omics Playground")
    titles[[3]] = c("Fasten Your Seat Belts!","Hi-speed analytics")
    titles[[4]] = c("Do-it-yourself Analytics","Yes you can!")
    titles[[5]] = c("Twenty-Four Seven","Your Playground doesn't go on coffee breaks")
    titles[[6]] = c("Analyze with confidence","Be a data rockstar, a Freddie Mercury of omics!")
    titles[[7]] = c("Play, Explore, Discover","Get deeper insights with Omics Playground")
    titles[[8]] = c("Skip the Queue","Take the fast lane. Self-service analytics.")
    titles[[9]] = c("Look Ma! No help!","I did it without a bioinformatician")
    titles[[10]] = c("Easy-peasy insight!","Get insight from your data the easy way")
    titles[[11]] = c("Zoom-zoom-insight!","Get faster insight from your data")
    titles[[12]] = c("Click-click-eureka!","Owe yourself that <i>eureka!</i> moment")
    titles[[13]] = c("I Love Omics Data!","Unleash your inner nerd with Omics Playground")
    titles[[14]] = c("More Omics Data","Is all I want for Christmas")
    titles[[15]] = c("Keep Exploring","Never stop discovering with Omics Playground")
    titles[[16]] = c("Real Bioinformaticians","Do it with Omics Playground (or in assembler...)")
    titles[[17]] = c("Real Biologists","Do it with Omics Playground")
    titles[[18]] = c("Ich bin doch nicht blöd!","Of course I use Omics Playground")
    titles[[19]] = c("Non sono mica scemo!","Of course I use Omics Playground")
    ## below from https://www.quotesweekly.com/keep-exploring-quotes/
    titles[[20]] = c("The Unexplored Plan","When you get into the whole field of exploring, you realize that we live on a relatively unexplored plan. – E. O. Wilson")
    titles[[21]] = c("Explore More","The more you explore, the more you learn and grow<br>– Nitesh Nishad")
    titles[[22]] = c("Discover New Oceans","Man cannot discover new oceans unless he has the courage to lose sight of the shore – Andre Gide")
    titles[[23]] = c("Adventurous Life","Love adventurous life. Be passionately curious about exploring new adventures. – Lailah Gifty Akita")
    titles[[24]] = c("Succes is Exploration","The first thing you have to find is the unknown. Learning is searching; anything else is just waiting. – Dale Daute")
    
    title <- titles[[length(titles)]]
    title <- sample(titles,1)[[1]]
    splash.title <- div(
        br(),br(),
        div(HTML(title[1]),style="font-size:70px;font-weight:700;line-height:1em;"),
        br(),
        div(HTML(title[2]),style="font-size:28px;"),
        br(),br(),br()
    )

    div.password <- div()
    div.email <- div()
    div.username <- div()
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
    div.alt <- div()
    if(!is.null(alt)) div.alt <- alt
    top <- HTML(rep("<br>",sum(c(!with.email,!with.username,!with.password))))
    
    ##splash.panel=div();ns=function(x)x
    splash.panel <- div(
        br(),br(),top,
        div.username,
        div.email,
        div.password,
        div.alt,
        br(),br(),
        ##actionButton(ns("login_btn"),NULL,icon=icon("sign-in-alt")),
        actionButton(ns("login_btn"),"Login",class="red-button")
    )
    body <- tagList(
        div(id="splash-title",splash.title),
        div(id="splash-panel",splash.panel)
    )
    m <- particlesSplashModal(body, ns=ns)
    return(m)
}

particlesSplashModal <- function(body, ns=NULL, easyClose=FALSE, fade=FALSE,
                                 buttons=TRUE, footer=TRUE)
{
    require(particlesjs)
    type <- type[1]    
    if(is.null(ns)) ns <- function(e) return(e)
    message("[AuthenticationDialog::splashModal]")

    div.footer = modalButton("Dismiss")
    if(buttons) {
        div.footer = tagList(
            actionButton(ns("action1"),"Read-the-docs", icon=icon("book"),
                         onclick="window.open('https://omicsplayground.readthedocs.io','_blank')"),
            actionButton(ns("action2"),"Watch tutorials", icon=icon("youtube"),
                         onclick="window.open('https://www.youtube.com/watch?v=_Q2LJmb2ihU&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-','_blank')"),
            actionButton(ns("action3"),"Get the source", icon=icon("github"),
                         onclick="window.open('https://github.com/bigomics/omicsplayground','_blank')"),
            actionButton(ns("action4"),"Docker image", icon=icon("docker"),
                         onclick="window.open('https://hub.docker.com/r/bigomics/omicsplayground','_blank')"),
            actionButton(ns("action5"),"User forum", icon=icon("users"),
                         onclick="window.open('https://groups.google.com/d/forum/omicsplayground','_blank')"),
            ##actionButton(ns("action_beer"),"Buy us a beer!", icon=icon("beer"),
            ##             onclick="window.open('https://www.buymeacoffee.com/bigomics','_blank
            actionButton(ns("action_beer"),"Buy us a coffee!", icon=icon("coffee"),
                         onclick="window.open('https://www.buymeacoffee.com/bigomics','_blank')")
            ## modalButton("Let's start!")
            ## actionButton(ns("action_play"), "Let's play!", class="red-button")
        )
    }
    if(!footer) {
        div.footer <- NULL
    }

    ## return modalDialog
    particlesjs.conf <- rjson::fromJSON(file="resources/particlesjs-config.json")
    m <- modalDialog(
        id = "modal-splash",
        div(
            id="particles-target",
            ##img(src = base64enc::dataURI(file="www/splash.png"),
            ##    width="100%", height="auto%", style="position:absolute;"),
            div(id="splash-logo", img(src=base64enc::dataURI(file="www/logo.png"),
                                      width=32,height=32)),
            body,
            div(id="splash-warning", textOutput(ns("login_warning")),style="color: red;"),
            style="height: 500px; width: 100%;"                
        ),
        footer = div.footer,
        particlesjs::particles(config=particlesjs.conf, target_id = "particles-target", timeout = 1000),
        size="l", easyClose=easyClose, fade=fade
    ) ## end of modalDialog
    
    return(m)
}


if(0) {
    
    ui <- fluidPage(
        plotOutput("results")
    )
    
    server <- function(input, output, session) {
        ##auth <- AuthenticationDialog(input, output, session, "password") 
        auth <- callModule(AuthenticationDialog, "auth", type="password")
        ## Render the plot
        output$results <- renderPlot({
            ##auth$showLogin()
            if(auth$logged()==FALSE) {
                AuthenticationUI("auth")
                message("[authenticateDialog] done! email=",email)
            } else {
                plot(sin)
                cat("class(auth)=",class(auth),"\n")
                cat("names(auth)=",names(auth),"\n")
                cat("class(auth$name)=",class(auth$name),"\n")
                title(paste0("hello ",auth$name(),"!"),cex.main=1.5)
            }
        })
    }
    
    shinyApp(ui, server)

}
