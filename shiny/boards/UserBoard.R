##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing UserBoard")

UserInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description"))
        ## shiny::uiOutput(ns("inputsUI"))
    )
}

UserUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("User profile",uiOutput(ns("userinfo_UI"))),
            shiny::tabPanel("Visitors map",uiOutput(ns("usersmap_UI")))
            ## shiny::tabPanel("Community forum",uiOutput(ns("forum_UI")))
        )
    )
}

UserBoard <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE

    user <- env[["load"]][["auth"]]

    dbg("[UserBoard] >>> initializing UserBoard...")
    
    ##-----------------------------------------------------------------------------
    ## Description
    ##-----------------------------------------------------------------------------
    
    output$description <- shiny::renderUI({

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
    
    usersmap_info = "<strong>Visitors map.</strong> The world map shows the number of users visiting this site by unique IP."
    
    shiny::callModule(
        plotModule,
        id = "usersmap", ## label="a", 
        plotlib = "baseplot",
        func = usersmap.RENDER,
        func2 = usersmap.RENDER, 
        info.text = usersmap_info,
        ##options = usersmap_options,
        pdf.width=12, pdf.height=7, pdf.pointsize=13,
        height = c(450,600), width = c('auto',1000), res=72,
        ##datacsv = enrich_getWordFreq,
        title = "Number of visitors by country",
        add.watermark = WATERMARK
    )

    ##usersmap_caption = "<b>(a)</b> <b>Geo locate.</b>"
    output$usersmapInfo <- shiny::renderUI({

        u <- ACCESS.LOG
        df <- u$visitors
        rownames(df) <-  df$country_name
        tot.users <- sum(df$count)
        freq <- df$count
        names(freq) <- df$country_name
        top.countries <- head(sort(freq,dec=TRUE),10)
        top.countriesTT <- paste("<li>",names(top.countries),top.countries,collapse=" ")
        
        shiny::HTML(
            "<b>Total visitors:</b>",tot.users,"<br><br>",
            "<b>Top 10 countries:</b><br><ol>",top.countriesTT,"</ol><br>",
            "<b>Period:</b><br>",u$period,"<br><br>"
        )
    })
    
    output$usersmap_UI <- shiny::renderUI({
        shiny::fillCol(
            height = 600,
            shiny::fillRow(
                flex = c(1,4.5),
                shiny::wellPanel( shiny::uiOutput(ns("usersmapInfo"))),
                plotWidget(ns("usersmap"))
            )
        )
    })


    ##---------------------------------------------------------------
    ##----------------- modules for Forum ---------------------------
    ##---------------------------------------------------------------
    
    output$forum <- shiny::renderUI({
        parenturl <- paste0(session$clientData$url_protocol,
                            "//",session$clientData$url_hostname,
                            ":",session$clientData$url_port,
                            session$clientData$url_pathname)
        ## parenturl <- gsub("localhost","127.0.0.1",parenturl)
        parenturl <- URLencode(parenturl, TRUE)
        cat("[LoadingBoard:forum] parenturl =",parenturl,"\n")
        src = paste0('https://groups.google.com/forum/embed/?place=forum/omicsplayground',
                     '&showsearch=true&showpopout=true&parenturl=',parenturl)
        cat("src = ",src,"\n")
        shiny::tags$iframe(id="forum_embed", src=src, height=600, width='100%',
                    ##seamless="seamless",
                    frameborder='no')
        ##HTML(src)
    })
         
    output$tweet <- shiny::renderUI({
        ## NOT WORKING YET...
        shiny::tags$a(class="twitter-timeline",
               href="https://twitter.com/bigomics?ref_src=twsrc%5Etfw")
        ##shiny::tags$script('twttr.widgets.load(document.getElementById("tweet"));')
    })
            
    output$forum_UI <- shiny::renderUI({
        shiny::fillCol(
            height = 550,
            shiny::fillRow(
                flex=c(4,0),
                shiny::htmlOutput(ns("forum"))
                ##uiOutput("tweet")
            )
        )
    })
    
    ##------------------------------------------------
    ## Board return object
    ##------------------------------------------------
    res <- list(
        enable_beta = reactive({ as.logical(input$enable_beta) })
    )
    return(res)
}
