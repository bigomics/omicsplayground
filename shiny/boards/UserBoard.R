##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing UserBoard")

UserInputs <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::uiOutput(ns("description"))
    )
}

UserUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    
    shiny::div(
        style = "padding-left:1rem!important;",
        shiny::h1("User Profile"),
        shiny::div(
            class = "row",
            shiny::div(
                class = "col-md-3",
                uiOutput(ns("username")),
                uiOutput(ns("plan"))
            )
        ),
        h3("Subscriptions"),
        shiny::actionButton(
            ns("manage"),
            "Manage Subscription"
        ),
        shiny::div(
            id = "user-subs"
        )
    )
}

UserBoard <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE

    user <- env[["load"]][["auth"]]

    dbg("[UserBoard] >>> initializing UserBoard...")

    observeEvent(user$logged(), {
        if(!user$logged())
            return()
        
        session$sendCustomMessage(
            "get-subs",
            list(
                ns = ns(NULL)
            )
        )
        user.name  <- user$name()
        user.level <- user$level()
        user.email <- user$email()
        dbg("[UserBoard::description] names(user) = ",names(user) )
        dbg("[UserBoard::description] user.name = ",user.name )
        dbg("[UserBoard::description] user.email = ",user.email )
        dbg("[UserBoard::description] user.level = ",user.level )        

        if(is.null(user.name))  user.name  <- ""
        if(is.null(user.email)) user.email <- ""
        
        description = "Signed in as<h2><b>NAME<b></h2><h4><b>EMAIL<b></h4><br><h4>LEVEL</h4>"
        description = "Signed in as<h2><b>NAME<b></h2><h4><b>EMAIL<b></h4>"
        ##description = "Signed in as<h4><b>EMAIL<b></h4>"
        dbg("[UserBoard::description] 1a : " )
        description <- sub("EMAIL", as.character(user.email), description)
        dbg("[UserBoard::description] 1b : " )        
        description <- sub("NAME", as.character(user.name), description)
        dbg("[UserBoard::description] 1c : " )        
        description <- sub("LEVEL", as.character(user.level), description)
        dbg("[UserBoard::description] 1d : " )        
        shiny::HTML(description)
    })
    
    ##-----------------------------------------------------------------------------
    ## User interface
    ##-----------------------------------------------------------------------------
    
    output$inputsUI <- shiny::renderUI({        
    })

    output$userdata <- renderTable({
        dbg("[UserBoard::userdata]  renderDataTable")
        values <- c(
            name   = user$name(),
            email  = user$email()
        )
        values[which(values=="")] <- "(not set)"
        data.frame(' '=names(values), '  '=values, check.names=FALSE)
    })

    output$userdata2 <- renderTable({
        dbg("[UserBoard::userdata]  renderDataTable")
        values <- c(
            plan   = user$level(),
            ##logged = user$logged(),
            limit  = paste(user$limit(),collapse=';')
        )
        values[which(values=="")] <- "(not set)"
        data.frame(' '=names(values), '  '=values, check.names=FALSE)
    })
        
    output$userinfo_UI <- shiny::renderUI({
        tagList(
            shiny::HTML("<h4>Personal</h4>"),
            shiny::tableOutput(ns("userdata")),
            shiny::br(),
            shiny::HTML("<h4>Account</h4>"),            
            shiny::tableOutput(ns("userdata2")),
            shiny::br(),            
            shiny::HTML("<h4>Settings</h4>"),            
            shinyWidgets::prettySwitch(ns("enable_beta"),"enable beta features")
            ##shinyWidgets::prettySwitch(ns("enable_alpha"),"enable alpha features")
        )
    })
    shiny::outputOptions(output, "userinfo_UI", suspendWhenHidden=FALSE) ## important!

    
    ##---------------------------------------------------------------
    ##--------------------- modules for UsersMap --------------------
    ##---------------------------------------------------------------
    
    usersmap.RENDER %<a-% shiny::reactive({
        
        df <- ACCESS.LOG$visitors        
        ## sPDF <- rworldmap::getMap()  
        ## rworldmap::mapCountryData(sPDF, nameColumnToPlot='continent')
        sPDF <- rworldmap::joinCountryData2Map(
            df,
            joinCode = "ISO2",
            nameJoinColumn = "country_code")
        
        par(mai=c(0,0.4,0.2,1),xaxs="i",yaxs="i")
        mapParams <- rworldmap::mapCountryData(
            sPDF, nameColumnToPlot="count",
            ##mapTitle = "Number of unique IPs",
            mapTitle = "", addLegend='FALSE',
            colourPalette = RColorBrewer::brewer.pal(9,"Blues"),
            numCats=9, catMethod="logFixedWidth")   
                   
        ##add a modified legend using the same initial parameters as mapCountryData
        do.call( rworldmap::addMapLegend,
                c(mapParams, labelFontSize = 0.85, legendWidth = 1.2, legendShrink = 0.5,
                  legendMar = 4, horizontal = FALSE, legendArgs = NULL, tcl = -0.5,
                  sigFigs = 4, digits = 3)
                )
        
    })

    output$username <- renderUI({
        div(
            h3(user$name()),
            p(user$email())
        )
    })

    output$plan <- renderUI({
        plan_class <- "info"
        if(user$level() == "premium")
            plan_class <- "success"

        cl <- sprintf("badge badge-%s", plan_class)
        p(
            span("Subscription level", style="color:grey;"),
            span(class = cl, tools::toTitleCase(user$level()))
        )
    })

    observeEvent(input$manage, {
        response <- httr::POST(
            "https://api.stripe.com/v1/billing_portal/sessions",
            body = list(
                customer = user$stripe_id(),
                return_url = user$href()
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
    
    ##------------------------------------------------
    ## Board return object
    ##------------------------------------------------
    res <- list(
        enable_beta = reactive({ as.logical(input$enable_beta) })
    )
    return(res)
}
