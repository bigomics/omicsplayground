##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2021 BigOmics Analytics Sagl. All rights reserved.
##

#' The main application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @export
app_ui <- function(request) {

    TABVIEWS <- list(
        "load"   = tabView("Home",LoadingInputs("load"),LoadingUI("load")),
        "view"   = tabView("DataView",DataViewInputs("view"),DataViewUI("view")),
        "clust"  = tabView("Cluster samples",ClusteringInputs("clust"),ClusteringUI("clust")),
        "ftmap"  = tabView("Cluster features",FeatureMapInputs("ftmap"),FeatureMapUI("ftmap")),    
        "wgcna"  = tabView("WGCNA (beta)",WgcnaInputs("wgcna"),WgcnaUI("wgcna")),
        "expr"   = tabView("Differential expression",ExpressionInputs("expr"),ExpressionUI("expr")),
        "cor"    = tabView("Correlation analysis", CorrelationInputs("cor"), CorrelationUI("cor")),
        "enrich" = tabView("Geneset enrichment",EnrichmentInputs("enrich"), EnrichmentUI("enrich")),
        "func"   = tabView("Pathway analysis", FunctionalInputs("func"), FunctionalUI("func")),
        "word"   = tabView("Word cloud", WordCloudInputs("word"), WordCloudUI("word")),
        "drug"   = tabView("Drug connectivity", DrugConnectivityInputs("drug"), DrugConnectivityUI("drug")),
        "isect"  = tabView("Compare signatures", IntersectionInputs("isect"), IntersectionUI("isect")),
        "sig"    = tabView("Test signatures", SignatureInputs("sig"), SignatureUI("sig")),
        "bio"    = tabView("Find biomarkers", BiomarkerInputs("bio"), BiomarkerUI("bio")),
        "cmap"   = tabView("Find similar experiments", ConnectivityInputs("cmap"), ConnectivityUI("cmap")),
        "scell"  = tabView("CellProfiling", SingleCellInputs("scell"), SingleCellUI("scell")),
        "tcga"   = tabView("TCGA survival (beta)", TcgaInputs("tcga"), TcgaUI("tcga")),
        "comp"   = tabView("Compare datasets (beta)", CompareInputs("comp"), CompareUI("comp"))
    )

    ##names(TABVIEWS)
    ##TABVIEWS <- TABVIEWS[names(TABVIEWS) %in% names(which(ENABLED))]
    ##names(TABVIEWS)

    ##-------------------------------------------------------
    ## Build USERMENU
    ##-------------------------------------------------------
    user.tab <-  tabView(title = "Settings", id="user", UserInputs("user"), UserUI("user"))    
    ##title = shiny::HTML("<span class='label label-info' id='authentication-user'></span>"),
    logout.tab  <- shiny::tabPanel(shiny::HTML("<a onClick='logout()' id='authentication-logout'>Logout</a>"))

    ## conditionally add if firebase authentication is enabled
    stop.tab    <- shiny::tabPanel(shiny::HTML("<a onClick='logout();quit();'>Quit</a>"))
    if(opt$AUTHENTICATION == "shinyproxy") {
        ## For ShinyProxy we need to redirect to /logout for clean session
        ## logout. Then we need a redirect to the /login page.
        logout.tab  <- shiny::tabPanel(shiny::HTML("<a href='/login' onClick='shinyproxy_logout();' id='authentication-logout'>Logout</a>"))    
    }

    upgrade.tab <- NULL
    if(opt$AUTHENTICATION == "firebase") {
        upgrade.tab <- shiny::tabPanel(shiny::HTML("<a onClick='show_plans()' style='font-weight:bold;color:#2a9d8f;cursor:pointer;' id='authentication-upgrade'>Upgrade</a>"))
    }

    user.menu <- shiny::navbarMenu(
         ##title="User",
         title=icon("user-circle","fa"),                     
         user.tab,
         upgrade.tab,
         "----",
         shiny::tabPanel(title=shiny::HTML("<a href='https://omicsplayground.readthedocs.io' target='_blank'>Documentation</a>")),
         shiny::tabPanel(title=shiny::HTML("<a href='https://www.youtube.com/watch?v=_Q2LJmb2ihU&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-' target='_blank'>Video tutorials</a>")),
         shiny::tabPanel(title=shiny::HTML("<a href='https://groups.google.com/d/forum/omicsplayground' target='_blank'>Community Forum</a>")),
         shiny::tabPanel(title=shiny::HTML("<a href='https://github.com/bigomics/omicsplayground' target='_blank'>GitHub</a>")),
         "----",         
         logout.tab,
         stop.tab
    )

    createUI <- function(tabs)
    {
        message("\n=======================================================================")
        message("================================ UI ===================================")
        message("=======================================================================\n")

        version <- scan("../VERSION", character())[1]
        TITLE = paste(opt$TITLE,version)
        LOGO = shiny::div(shiny::img(src="bigomics-logo-white-48px.png", height="48px"),
                          TITLE, id="navbar-logo", style="margin-top:-13px;")    
        title = shiny::tagList(LOGO)
        windowTitle = TITLE
        theme = shinythemes::shinytheme("cerulean")
        id = "maintabs"


        ## Add Google Tag manager header code    
        gtag <- NULL
        if(Sys.getenv("OMICS_GOOGLE_TAG")!="") {
            gtag.html <- htmltools::includeHTML("www/google-tags.html")
            gtag.html <- sub("GTM-0000000",Sys.getenv("OMICS_GOOGLE_TAG"),gtag.html)
            gtag <- shiny::tags$head(gtag.html)
        }
        
        header = shiny::tagList(
                            shiny::tags$head(shiny::tags$script(src="temp.js")),
                            shiny::tags$head(shiny::tags$script(src="bigomics-extra.js")),  ## chatra,clarity
                            gtag,   ## Google Tags???
                            shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "playground.css")),
                            shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "fonts.css")),
                            shiny::tags$head(shiny::tags$link(rel="shortcut icon", href="favicon.ico")),
                            shinyjs::useShinyjs(),
                            sever::useSever(),
                            shinylogs::use_tracking(),
                            shinyalert::useShinyalert(),  # Set up shinyalert
                            firebase::useFirebase(firestore = TRUE),
                            ##TAGS.JSSCRIPT,  ## window size
                            shiny::tags$script(async=NA, src="https://platform.twitter.com/widgets.js"),
                            shiny::div(shiny::textOutput("current_dataset"), class='current-data'),
                            shiny::div(class='label label-info current-user',id='authentication-user')   
                            ##QuestionBoard_UI("qa")
                        )
        names(header) <- NULL
        
        footer.gif = shiny::tagList(
                                shinybusy::busy_start_up(
                                               text = "\nPrepping your Omics Playground...", mode = "auto",
                                               background="#2780e3", color="#ffffff",
                                               ##loader = shiny::img(src=base64enc::dataURI(file="www/ready.png"))
                                               loader = shiny::img(src=base64enc::dataURI(file="www/monster-hi.png"))            
                                           )
                            )
        footer = footer.gif
        
        ##-------------------------------------
        ## create TAB list
        ##-------------------------------------
        createNavbarMenu <- function(title, tabs, icon=NULL) {
            tablist <- TABVIEWS[tabs]
            names(tablist) <- NULL
            do.call( navbarMenu, c(tablist, title=title, icon=icon) )
        }
        
        ##tablist <- TABVIEWS[tabs]
        tablist <- list()
        i=1
        for(i in 1:length(tabs)) {
            itab <- tabs[[i]]
            itab <- itab[which(ENABLED[itab])] ## only enabled
            if(length(itab)>1) {
                message("[MAIN] creating menu items for: ",paste(itab,collapse=" "))
                m <- createNavbarMenu( names(tabs)[i], itab )
                tablist[[i]] <- m
            } else if(length(itab)==1) {
                message("[MAIN] creating menu item for: ",itab)
                tablist[[i]] <- TABVIEWS[[itab]] 
            } else {
                
            }
        }
        tablist <- tablist[!sapply(tablist,is.null)]
        
        ## add user menu (profile, help + logout)
        tablist[["usermenu"]] <- user.menu
        
        ##-------------------------------------
        ## create navbarPage
        ##-------------------------------------
        selected = "Home"    
        names(tablist) <- NULL
        do.call( navbarPage, c(tablist,
                               title=title, id=id,
                               selected=selected,
                               windowTitle = windowTitle,
                               header = shiny::tagList(header),
                               footer = shiny::tagList(footer),
                               theme = theme))
    }

    tabs = list(
        "Home" = c("load"),
        "DataView" = "view",
        "Clustering" = c("clust","ftmap","wgcna"),
        "Expression" = c("expr","cor"),
        "Enrichment" = c("enrich","func","word","drug"),
        "Signature" = c("isect","sig","bio","cmap","comp","tcga"),
        "CellProfiling" = "scell",
        "DEV" = c("corsa","system","multi")
    )

    gtag2 <- NULL
    if(Sys.getenv("OMICS_GOOGLE_TAG")!="") {
        ## Add Google Tag manager body code
        gtag2 <- htmltools::includeHTML("www/google-tags-noscript.html")
        gtag2 <- sub("GTM-0000000",Sys.getenv("OMICS_GOOGLE_TAG"),gtag2)
    }

    tagList(
        gtag2,
        createUI(tabs)
    )

}
