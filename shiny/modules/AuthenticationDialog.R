
if(0) {
    type="password"
    CREDENTIALS = data.frame(username=c("ivo","stefan"),
                             password=c("iii","sss"),
                             expiry=c("2025-01-01","2020-01-01"),
                             row.names=1)
}

AuthenticationUI <- function(id) {
    ns <- NS(id)  ## namespace
    showModal(uiOutput(ns("showLogin")))
}

AuthenticationDialog <- function(input, output, session,
                                 type="password", credentials=NULL) {

    ns <- session$ns
    
    USER <- reactiveValues(logged=FALSE, name=NA, email=NA, valid=NA)
    SMTP_SERVER  <- 
    
    if(type=="password") {

        output$showLogin <- renderUI({
            modalDialog(
                id = "login",
                title = "Login to Omics Playground",
                tagList( 
                    textInput(ns("login_username"), "Username:"),
                    passwordInput(ns("login_password"), "Password:"),
                    div( actionButton(ns("login_btn"), "Login"),
                        style="text-align: center;")
                ),
                footer = textOutput(ns("login_warning")),
                size = "s"
            )
        })


    } else if(type=="register") {    
        
        output$showLogin <- renderUI({
            modalDialog(
                id = "login",
                title = "Login to Omics Playground",
                tagList( 
                    p("Please enter your email address below to login."),
                    textInput(ns("login_username"), "E-mail:"),
                    div( actionButton(ns("login_btn"), "Login"),
                        style="text-align: center;")
                ),
                footer = div(textOutput(ns("login_warning")),style="color: red;"),
                size = "s"
            )
        })
        
        observeEvent( input$register_link, {
            ##install.packages("countrycode")
            require(countrycode)
            all.countries = countrycode::codelist[,"country.name.en"]
            
            showModal( modalDialog(
                id = "register",
                title = "Omics Playground registration",
                tagList( 
                    p("Fill in the form below. Registration data is used to measure usage only. This helps us track and better serve our user community."),
                    ## textInput("register_name", "Name:"),
                    textInput("register_email", "E-mail:"),
                    ## textInput("register_organization", "Organization:"),
                    selectInput("register_country", "Country:", choices=c("",all.countries)),
                    div( ##actionButton(ns("register_btn_skip"), "Skip"),
                        actionButton(ns("register_btn"), "Register"),
                        style="text-align: center;")
                ),
                footer = div(textOutput("register_warning"),style="color: red;"),
                size = "s"
            ))
        })           

        observeEvent( input$register_btn, {

            email.OK = ( !is.null(input$register_email) &&
                         input$register_email != "" &&
                         grepl("@",input$register_email))  ## valid email
            country.OK = (!is.null(input$register_country) &&
                          input$register_country != "" )
            register.OK = ( email.OK && country.OK )

            if(register.OK) {
                rdata <- paste("tom@acme.com","Tom Cruise","ACME Inc.","USA",date(),sep=",")
                rdata <- paste(input$register_email,
                               input$register_name,
                               input$register_organization,
                               input$register_country,
                               date(), sep=",")
                write( rdata, file="../logs/registered.csv", append=TRUE )

                USER$name   <- input$register_email
                USER$logged <- TRUE  ## global reactive
                show("register_warning")
                output$register_warning = renderText(paste("Welcome",input$register_email,"!"))
                shinyjs::delay(3000, hide("register_warning", anim = TRUE, animType = "fade"))
                removeModal()
                
            } else {
                ##show("register_warning")
                output$register_warning = renderText("Invalid email or country")
                ##shinyjs::delay(2000, hide("register_warning", anim = TRUE, animType = "fade"))
                shinyjs::delay(2000, {output$register_warning <- renderText("")})
            }
        }) 

        ## observe-event end-if-REGISTER
    } else {
        output$showLogin <- renderUI({})
    }

    output$login_warning = renderText("")

    observeEvent( input$logout, {
        ##updateTextInput(session, ".username", value=NULL)
        reset(ns("login_username"))
        reset(ns("login_password"))
        USER$logged <- FALSE
    })

    observeEvent( input$login_btn, {           

        cat("type=",type,"\n")
        cat("username=",input$login_username,"\n")
        cat("password=",input$login_password,"\n")
        cat("Logged=",USER$logged,"\n")

        login.OK   = FALSE
        valid.date = FALSE
        valid.user = FALSE
        
        if(type=="register") {
            ##if( is.null(input$login_username) || is.null(input$login_password)) return(NULL)
            ##if( input$login_username=="" || input$login_password=="") return(NULL)    
            ##if( is.null(input$login_name) || input$login_name == "") return(NULL)
            if( is.null(input$login_username) || input$login_username=="") return(NULL)
            if(file.exists("../logs/registered.csv")) {
                registered <- read.csv(file="../logs/registered.csv",header=FALSE,
                                       stringsAsFactors=FALSE)
                registered.users <- unique(registered[,1]) ## emails
                cat("registered.users=",registered.users,"\n")
                registered.users <- c(registered.users,"demo@bigomics.ch")
                login.OK = (input$login_username %in% registered.users)
            }
        }
        
        if(type=="password") {
            if( is.null(input$login_username) || is.null(input$login_password)) return(NULL)
            if( input$login_username=="" || input$login_password=="") return(NULL)    
            username <- input$login_username
            valid.user <- isTRUE(username %in% rownames(credentials))
            valid.pw   <- isTRUE(credentials[username,"password"]==input$login_password)
            valid.date <- isTRUE( Sys.Date() < as.Date(credentials[username,"expiry"]) )
            login.OK = (valid.user && valid.pw && valid.date)
            message("--------- password login ---------")
            message("input.username=",input$login_username)
            message("input.password=",input$login_password)
            message("user.password=",credentials[username,"password"])
            message("user.expiry=",credentials[username,"expiry"])
            message("valid.user=",valid.user)
            message("valid.date=",valid.date)
            message("valid.pw=",valid.pw)
            message("----------------------------------")
        }
        
        if (login.OK) {
            output$login_warning = renderText("")
            removeModal()
            USER$name   <- input$login_username
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
            ##show("login_warning")
            if(type=="password") {
                if(!valid.date) {
                    output$login_warning = renderText("Registration expired")
                } else {
                    output$login_warning = renderText("Invalid username or password")
                }
            }
            if(type=="register") {
                output$login_warning = renderText("Email address not recognized")
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
