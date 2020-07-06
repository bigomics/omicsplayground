
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

AuthenticationDialog <- function(input, output, session,
                                 credentials.file)
{
    ns <- session$ns
    
    credentials.file    
    if(!file.exists(credentials.file)) {
        ee <- data.frame(email="", name="", password="",
                         date=NA, expiry=NA, verified=FALSE)
        write.table(ee[0,], file=credentials.file, append=FALSE, sep=",",
                    row.names=FALSE, col.names=TRUE)
    }
    CREDENTIALS <- read.csv(credentials.file, stringsAsFactors=FALSE)
    head(CREDENTIALS)
    
    USER <- reactiveValues(logged=FALSE, name="", email="",
                           password=NA, verified=FALSE)
    
    SMTP_SERVER  <- Sys.getenv("SMTP_SERVER")
    SMTP_USER    <- Sys.getenv("SMTP_USER")
    SMTP_SECRET  <- Sys.getenv("SMTP_SECRET")
    has.smtp <- SMTP_SERVER!="" && SMTP_USER!="" && SMTP_SECRET!=""
    has.smtp
    message("[AuthenticationDialog] has.smtp = ",has.smtp)
    
    DIALOG_SLEEP = 0.3

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
        cat("[AuthenticationDialog::login] email = ",input$login_email,"\n")
        
        ## check login
        email <- input$login_email
        email.exists <- email %in% CREDENTIALS$email
        email.verified <- FALSE
        ##USER <- reactiveValues(logged=FALSE, name=NA, email=NA, verified=FALSE)
        if(email.exists) {
            sel <- tail(which(CREDENTIALS$email == email),1)
            email.verified <- CREDENTIALS[sel,"verified"]
            USER$name <- CREDENTIALS[sel,"name"]
        }
        USER$email <- email
        USER$verified <- email.verified
        USER$logged <- TRUE
        
        cat("[AuthenticationDialog::login] email.exists = ",email.exists,"\n")
        cat("[AuthenticationDialog::login] email.verified = ",email.verified,"\n")
        
        if(email.exists && email.verified) {

            for(i in 1:10) removeModal()
            Sys.sleep(DIALOG_SLEEP)
            updateCredentials(credentials.file, USER)
            showDialog("welcome")
            
        } else  {

            auth.code <- paste0(sample(c("A","T","G","C"),12,replace=TRUE),collapse="")
            auth.code        
            USER$password <- auth.code
            
            if(has.smtp) {
                smtpAuthCode(email, auth.code)
                showDialog("verify")
            } else {            
                showDialog("captcha", auth.code)
            }
        }
    })

    observeEvent( input$verify_btn, {
        message("[AuthenticationDialog::register] verify_btn pressed")
        message("[AuthenticationDialog::register] USER$email = ",USER$email)
        message("[AuthenticationDialog::register] USER$password = ",USER$password)        
        message("[AuthenticationDialog::register] verify_code = ",input$verify_code)

        USER$verified = (USER$password == input$verify_code)
        message("[AuthenticationDialog::register] verified = ",USER$verified)        
        
        if(USER$verified) {
            output$verify_warning = renderText("")
            output$captcha_warning = renderText("")
            ##removeModal()	
            for(i in 1:10) removeModal()
            Sys.sleep(DIALOG_SLEEP)
            updateCredentials(credentials.file, USER)
            showDialog("welcome")
        } else {
            output$verify_warning = renderText("Invalid verification")
            output$captcha_warning = renderText("Invalid verification")
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
        if(stage=="login") {
            m <- modalDialog(
                id = "auth_dialog",
                title = "Log in to your Playground",
                tagList( 
                    p( "Enter your email to log in to Omics Playground."),
                    textInput(ns("login_email"), "E-mail:"),
                    div(actionButton(ns("login_btn"), "Login"),style="text-align: center;")
                ),
                footer = NULL,
                size = "s"
            )
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

            tmpfile <- tempfile()
            png(tmpfile, w=400, h=50)
            par(mar=c(0,0,0,0))
            frame()
            graphics::text(0.5,0.5,data,cex=4)
            dev.off()
            
            m <- modalDialog(
                id = "auth_dialog",
                title = "Verify captcha",
                tagList( 
                    p("Please copy the letters into the captcha below:"),
                    div(img(src=base64enc::dataURI(file=tmpfile), width="250px"), style="text-align: center;"),
                    br(),
                    textInput(ns("verify_code"), "Enter captcha:"),                    
                    div( actionButton(ns("verify_btn"), "Next"), style="text-align: center;")
                ),
                footer = div(textOutput(ns("captcha_warning")),style="color: red;"),                
                size = "s"
            )
            unlink(tmpfile)
        } else if(stage=="welcome") {

            message("[AuthenticationDialog::showDialog] welcome ",USER$name)            
            all.hello = c("Hello","Salut","Hola","Privet","Nǐ hǎo","Ciao","Hi","Hoi","Hej",
                          "Yassou","Selam","Hey","Hei")
            title = paste( paste0(sample(all.hello,3),"!"), collapse=" ")
            if(!is.na(USER$name) && !USER$name %in% c("NA","")) {
                first.name <- strsplit(as.character(USER$name),split=" ")[[1]][1]
                first.name <- paste0(toupper(substring(first.name,1,1)),substring(first.name,2,999))
                title = paste0(sample(all.hello,1)," ",first.name,"!")
            }
            toon <- pgx.randomCartoon()
            m <- modalDialog(
                id = "auth_dialog",
                title = title,
                tagList( 
                    HTML("Welcome to the Playground! ",
                         "We are continuously improving the Playground and have a bunch of new features. ",
                         "If you have any questions, please ask our ",
                         "<a href='https://groups.google.com/forum/#!forum/omicsplayground'>user forum</a>. ",
                         "We hope you will discover lots of interesting things today! "),
                    br(),br(),
                    div(img(src=base64enc::dataURI(file=toon$img), width="475px"), style="text-align: center;")
                ),
                ## footer = NULL,
                easyClose = TRUE,
                size = "m"
            )
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

    updateCredentials <- function(credentials.file, user) {

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
            verified = user$verified)

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

PasswordAuthenticationDialog <- function(input, output, session,
                                         credentials.file)
{
    ns <- session$ns    
    USER <- reactiveValues(logged=FALSE, name=NA, email=NA, valid=NA)

    credentials <- read.csv(credentials.file)
    dim(credentials)

    output$showLogin <- renderUI({
        modalDialog(
            id = "login",
            title = "Login to Omics Playground",
            tagList( 
                textInput(ns("login_email"), "E-mail:"),
                passwordInput(ns("login_password"), "Password:"),
                div( actionButton(ns("login_btn"), "Login"),
                    style="text-align: center;")
            ),
            footer = textOutput(ns("login_warning")),
            size = "s"
        )
    })

    output$login_warning = renderText("")

    observeEvent( input$logout, {
        ##updateTextInput(session, ".username", value=NULL)
        reset(ns("login_email"))
        reset(ns("login_password"))
        USER$logged <- FALSE
    })

    observeEvent( input$login_btn, {           

        cat("email=",input$login_email,"\n")
        cat("password=",input$login_password,"\n")
        cat("Logged=",USER$logged,"\n")

        login.OK   = FALSE
        valid.date = FALSE
        valid.user = FALSE
        
        
        if( is.null(input$login_email) || is.null(input$login_password)) return(NULL)
        if( input$login_email=="" || input$login_password=="") return(NULL)    
        email <- input$login_email
        sel <- tail(which( credentials$email == email),1)
        valid.user <- isTRUE(email %in% credentials$email)
        valid.pw   <- isTRUE(credentials[sel,"password"]==input$login_password)
        valid.date <- isTRUE( Sys.Date() < as.Date(credentials[sel,"expiry"]) )
        login.OK = (valid.user && valid.pw && valid.date)
        message("--------- password login ---------")
        message("input.email=",input$login_email)
        message("input.password=",input$login_password)
        message("user.password=",credentials[sel,"password"])
        message("user.expiry=",credentials[sel,"expiry"])
        message("valid.user=",valid.user)
        message("valid.date=",valid.date)
        message("valid.pw=",valid.pw)
        message("----------------------------------")

        
        if (login.OK) {
            output$login_warning = renderText("")
            removeModal()
            ##USER$name   <- input$login_username
            USER$email <- input$login_email
            USER$logged <- TRUE
            
            ## Here you can perform some user-specific functions, or site news
            if(0) {
                removeModal()
                showModal(modalDialog(
                    paste("Welcome",USER$name,"!"),
                    footer=NULL, size="m"))
                Sys.sleep(3)
            }
            removeModal()
            
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

    rt <- list(
        name = reactive(USER$name),
        logged = reactive(USER$logged)
    )
    return(rt)
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
