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

if(DEV) {
    if(ENABLED["corsa"]) TABVIEWS$corsa = tabView("CORSA (dev)",CorsaInputs("corsa"),
                                                  CorsaUI("corsa"))
    if(ENABLED["system"]) TABVIEWS$system = tabView("Systems analysis (dev)",
                                                    SystemInputs("system"),SystemUI("system"))
    if(ENABLED["multi"]) TABVIEWS$multi = tabView("Multi-level (dev)", MultiLevelInputs("multi"),
                                                  MultiLevelUI("multi"))
}

names(TABVIEWS)
TABVIEWS <- TABVIEWS[names(TABVIEWS) %in% names(which(ENABLED))]
names(TABVIEWS)

#-------------------------------------------------------
## Build USERMENU
#-------------------------------------------------------
user.tab <-  tabView(title = "Settings", id="user", UserInputs("user"), UserUI("user"))    
##title = shiny::HTML("<span class='label label-info' id='authentication-user'></span>"),
# logout.tab  <- shiny::tabPanel(shiny::HTML("<a onClick='logout()' id='authentication-logout'>Logout</a>"))
logout.tab <- bigdash::navbarDropdownItem(
    "Logout",
    onClick = "logout"
)

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
    upgrade.tab <- bigdash::navbarDropdownItem(
        "Upgrade",
        onClick = "show_plans()"
    )
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
    message("\n======================================================")
    message("======================= UI ===========================")
    message("======================================================\n")

    version <- scan("../VERSION", character())[1]
    id = "maintabs"
    gtag <- NULL
    if(0 && file.exists("www/google-tags.html")) {
        gtag <- shiny::tags$head(includeHTML("www/google-tags.html"))
    }
        
    header = shiny::tagList(
        shiny::tags$head(shiny::tags$script(src="temp.js")),
        shiny::tags$head(shiny::tags$script(src="bigomics-extra.js")),  ## chatra,clarity
        gtag,   ## Google Tags???
        shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "playground.css")),
        shiny::tags$head(shiny::tags$link(rel="shortcut icon", href="favicon.ico")),
        shinyjs::useShinyjs(),
        sever::useSever(),
        shinylogs::use_tracking(),
        shinyalert::useShinyalert(),  # Set up shinyalert
        firebase::useFirebase(firestore = TRUE),
        shiny::tags$script(async=NA, src="https://platform.twitter.com/widgets.js"),
        # shiny::div(shiny::textOutput("current_dataset"), class='current-data'),
        shiny::div(class='label label-info current-user',id='authentication-user')        
    )
    
    footer <- shiny::tagList(
        shinybusy::busy_start_up(
            text = "\nPrepping your Omics Playground...", mode = "auto",
            background="#2780e3", color="#ffffff",
            ##loader = shiny::img(src=base64enc::dataURI(file="www/ready.png"))
            loader = shiny::img(src=base64enc::dataURI(file="www/monster-hi.png"))            
        )
    )
    
    #return(ui)

    logout.tab <- bigdash::navbarDropdownItem(
        "Logout",
        onClick = "logout()"
    )

    if(opt$AUTHENTICATION == "shinyproxy") {
        ## For ShinyProxy we need to redirect to /logout for clean session
        ## logout. Then we need a redirect to the /login page.
        logout.tab <- bigdash::navbarDropdownItem(
            "Logout",
            onClick = "shinyproxy_logout();",
            link = "/login"
        )
    }
    bigdash::bigPage(
        header,
        navbar = bigdash::navbar(
            tags$img(
                src = "assets/img/bigomics.png",
                width = "110",
            ),
            bigdash::navbarDropdown(
                "User",
                bigdash::navbarDropdownTab(
                    "Settings",
                    "userSettings"
                ),
                upgrade.tab,
                bigdash::navbarDropdownItem(
                    "Documentation",
                    link = "https://omicsplayground.readthedocs.io"
                ),
                bigdash::navbarDropdownItem(
                    "Video tutorials",
                    link = "https://www.youtube.com/watch?v=_Q2LJmb2ihU&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-"
                ),
                bigdash::navbarDropdownItem(
                    "Community Forum",
                    link = "https://groups.google.com/d/forum/omicsplayground"
                ),
                bigdash::navbarDropdownItem(
                    "Github",
                    link = "https://github.com/bigomics/omicsplayground"
                ),
                logout.tab
            )
        ),
        sidebar = bigdash::sidebar(
            bigdash::sidebarItem(
                "Home",
                "home-tab"
            ),
            bigdash::sidebarItem(
                "DataView",
                "dataview-tab"
            ),
            bigdash::sidebarMenu(
                "Clustering",
                bigdash::sidebarMenuItem(
                    "Cluster samples",
                    "clustersamples-tab"
                ),
                bigdash::sidebarMenuItem(
                    "Cluster features",
                    "clusterfeaters-tab"
                ),
                bigdash::sidebarMenuItem(
                    "WGCNA (beta)",
                    "wgcna-tab"
                )
            ),
            bigdash::sidebarMenu(
                "Expression",
                bigdash::sidebarMenuItem(
                    "Differential expression",
                    "diffexpr-tab"
                ),
                bigdash::sidebarMenuItem(
                    "Correlation analysis",
                    "corr-tab"
                )
            ),
            bigdash::sidebarMenu(
                "Enrichment",
                bigdash::sidebarMenuItem(
                    "Geneset enrichment",
                    "enrich-tab"
                ),
                bigdash::sidebarMenuItem(
                    "Pathway analysis",
                    "pathway-tab"
                ),
                bigdash::sidebarMenuItem(
                    "Word cloud",
                    "cloud-tab"
                ),
                bigdash::sidebarMenuItem(
                    "Drug connectivity",
                    "drug-tab"
                )
            ),
            bigdash::sidebarMenu(
                "Signature",
                bigdash::sidebarMenuItem(
                    "Compare signatures",
                    "isect-tab"
                ),
                bigdash::sidebarMenuItem(
                    "Test signatures",
                    "sig-tab"
                ),
                bigdash::sidebarMenuItem(
                    "Find biomarkers",
                    "bio-tab"
                ),
                bigdash::sidebarMenuItem(
                    "Find similar experiments",
                    "cmap-tab"
                ),
                bigdash::sidebarMenuItem(
                    "Compare datasets (beta)",
                    "comp-tab"
                ),
                bigdash::sidebarMenuItem(
                    "TCGA survival (beta)",
                    "tcga-tab"
                )
            ),
            bigdash::sidebarItem(
                "Cell profiling",
                "cell-tab"
            )
        ),
        bigdash::bigTabs(
            bigdash::bigTabItem(
                "home-tab",
                LoadingInputs("load"),
                LoadingUI("load")
            ),
            bigdash::bigTabItem(
                "dataview-tab",
                DataViewInputs("view"),
                DataViewUI("view")
            ),
            bigdash::bigTabItem(
                "clustersamples-tab",
                ClusteringInputs("clust"),
                ClusteringUI("clust")
            ),
            bigdash::bigTabItem(
                "wgcna-tab",
                WgcnaInputs("wgcna"),
                WgcnaUI("wgcna")
            ),
            bigdash::bigTabItem(
                "diffexpr-tab",
                ExpressionInputs("expr"),
                ExpressionUI("expr")
            ),
            bigdash::bigTabItem(
                "corr-tab",
                CorrelationInputs("cor"),
                CorrelationUI("cor")
            ),
            bigdash::bigTabItem(
                "enrich-tab",
                EnrichmentInputs("enrich"),
                EnrichmentUI("enrich")
            ),
            bigdash::bigTabItem(
                "pathway-tab",
                FunctionalInputs("func"),
                FunctionalUI("func")
            ),
            bigdash::bigTabItem(
                "cloud-tab",
                WordCloudInputs("word"),
                WordCloudUI("word")
            ),
            bigdash::bigTabItem(
                "drug-tab",
                DrugConnectivityInputs("drug"),
                DrugConnectivityUI("drug")
            ),
            bigdash::bigTabItem(
                "isect-tab",
                IntersectionInputs("isect"),
                IntersectionUI("isect")
            ),
            bigdash::bigTabItem(
                "sig-tab",
                SignatureInputs("sig"),
                SignatureUI("sig")
            ),
            bigdash::bigTabItem(
                "bio-tab",
                BiomarkerInputs("bio"),
                BiomarkerUI("bio")
            ),
            bigdash::bigTabItem(
                "cmap-tab",
                ConnectivityInputs("cmap"),
                ConnectivityUI("cmap")
            ),
            bigdash::bigTabItem(
                "comp-tab",
                CompareInputs("comp"),
                CompareUI("comp")
            ),
            bigdash::bigTabItem(
                "tcga-tab",
                TcgaInputs("tcga"),
                TcgaUI("tcga")
            ),
            bigdash::bigTabItem(
                "cell-tab",
                SingleCellInputs("scell"),
                SingleCellUI("scell")
            ),
            bigdash::bigTabItem(
                "userSettings",
                UserInputs("user"), 
                UserUI("user")
            )
        ),
        tagList(footer)
    )
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

ui = createUI(tabs)
