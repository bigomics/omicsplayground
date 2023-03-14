app_ui <- function() {
    #-------------------------------------------------------
    ## Build USERMENU
    #-------------------------------------------------------
    user.tab <-  tabView(title = "Settings", id="user", UserInputs("user"), UserUI("user"))

    upgrade.tab <- NULL
    if(opt$AUTHENTICATION == "firebase") {
        upgrade.tab <- bigdash::navbarDropdownItem(
            "Upgrade",
            onClick = "show_plans()"
        )
    }

    gtag2 <- NULL
    if(Sys.getenv("OMICS_GOOGLE_TAG")!="") {
        ## Add Google Tag manager body code
        gtag2 <- htmltools::includeHTML("www/google-tags-noscript.html")
        gtag2 <- sub("GTM-0000000",Sys.getenv("OMICS_GOOGLE_TAG"),gtag2)
    }

    createUI <- function()
    {
        message("\n======================================================")
        message("======================= UI ===========================")
        message("======================================================\n")

        version <- scan(file.path(OPG,"VERSION"), character())[1]
        id = "maintabs"
        header = shiny::tagList(
            shiny::tags$head(shiny::tags$script(src="temp.js")),
            shiny::tags$head(shiny::tags$script(src="bigomics-extra.js")),  ## chatra,clarity
            gtag2,   ## Google Tags???
            shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "styles.min.css")),
            shiny::tags$head(shiny::tags$link(rel="shortcut icon", href="favicon.ico")),
            shinyjs::useShinyjs(),
            sever::useSever(),
            # shinylogs::use_tracking(),
            firebase::useFirebase(firestore = TRUE, analytics = TRUE),
            ##shiny::div(class='label label-info current-user',id='authentication-user'),
            shiny::tags$script(async=NA, src="https://platform.twitter.com/widgets.js")
        )

        footer <- shiny::tagList(
            SocialMediaModuleUI("socialmodal"),
            SendReferralModuleUI("sendreferral"),
            shinybusy::busy_start_up(
                text = "\nPrepping your personal playground...", mode = "auto",
                background="#2780e3", color="#ffffff",
                loader = shiny::img(
                    src="static/ready.png"
                )
            )
        )

        logout.tab <- bigdash::navbarDropdownItem(
            "Logout",
            onClick = "logoutInApp()"
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

        
        menu_tree = list(
          ## "Home" = c("load"),
          "Load" = c(
            welcome = "Welcome",
            load    = "Load dataset",
            upload  = "Upload data"                        
          ),          
          "DataView" = c(
            dataview = "DataView"
          ),
          "Clustering" = c(
            clustersamples  = "Samples",
            clusterfeatures = "Features",
            wgcna = "WGCNA"),
          "Expression" = c(
            diffexpr = "Differential expression",
            corr = "Correlation analysis"
          ),
          "Enrichment" = c(
            enrich = "Geneset Enrichment",
            pathway = "Pathway analysis",
            wordcloud = "Word cloud",
            drug = "Drug connectivity"
          ),
          "Signature" = c(
            isect = "Compare signatures",
            sig = "Test signatures",
            bio = "Find biomarkers",
            cmap = "Similar experiments",
            comp = "Compare datasets (beta)",
            tcga = "TCGA survival (beta)"
          ),
          "CellProfiling" = c(
            cell = "Cell profiling"
          )
          ## "DEV" = c("corsa","system","multi")
        )

        ## filter disabled modules
        ENABLED['welcome'] <<- TRUE
        ENABLED['load'] <<- TRUE
        #ENABLED['upload'] <- TRUE
        dbg("[ui.R] sum.enabled = ",sum(ENABLED))
        dbg("[ui.R] names.enabled = ",names(ENABLED))
        menu_tree <- lapply(menu_tree, function(m) m[which(ENABLED[names(m)])])
        
        populateSidebar <- function(menu_tree) {

          sidebar_item <- function(title, name) {
            div(class='sidebar-item',  bigdash::sidebarItem( title, paste0(name,"-tab")))
          }
          sidebar_menu_item <- function(title, name) {
            bigdash::sidebarMenuItem( title, paste0(name,"-tab"))
          }
          sidebar_menu_with_items <- function(tabs, title) {
            ee <- list()
            for(i in 1:length(tabs)) {
              tab.name  <- names(tabs)[i]
              tab.title <- tabs[i]              
              ee[[i]] <- sidebar_menu_item(tab.title, tab.name)
            }
            bigdash::sidebarMenu(title, !!!ee)
          }
          
          ## This creates the menu from a menu_tree
          menu <- list()
          i=3
          for(i in 1:length(menu_tree)) {
            tab.names <- names(menu_tree[[i]])
            tab.titles <- menu_tree[[i]]
            menu.id <- names(menu_tree)[i]
            if(length(tab.names)==0) {
              ## 
            } else if(length(tab.names)==1) {              
              menu[[menu.id]] <- sidebar_item(tab.titles, tab.names)
            } else {
              menu[[menu.id]] <- sidebar_menu_with_items(menu_tree[[i]], menu.id)
            }
          }
          ##lapply(menu,cat)
          return(menu)
        }
        
        info("[ui.R] creating sidebar menu")
        mm <- populateSidebar(menu_tree)
        mm <- lapply(mm, as.character)
        mm <- HTML(unlist(mm))
        sidebar <- bigdash::sidebar("Menu", mm)
        
        sidebar.save = bigdash::sidebar(
          "Menu",
          bigdash::sidebarMenu(
            "Load",
            bigdash::sidebarMenuItem(
              "Welcome",
              "welcome-tab"
            ),
            withTooltip(div(
              bigdash::sidebarMenuItem(
                "Load dataset",
                "load-tab"
              )),
              "This panel shows the available datasets within the platform. These data sets
              have been pre-computed and are ready to be used. Select a dataset in the table
              and load the data set by clicking the 'load' button.",
              placement = "top"             
            ),
            bigdash::sidebarMenuItem(
              "Upload data",
              "upload-tab"
            )
          ),
          withTooltip(
             div(class="sidebar-item",
               bigdash::sidebarItem(
                 "DataView",
                 "dataview-tab"
               )),
             "Information and descriptive statistics to quickly lookup a gene,
              check your experiment QC, view the raw data, sample or contrast tables.",
             placement = "top"
           ),
             bigdash::sidebarMenu(
                 "Clustering",
                 bigdash::sidebarMenuItem(
                     "Samples",
                     "clustersamples-tab"
                 ),
                 bigdash::sidebarMenuItem(
                     "Features",
                     "clusterfeatures-tab"
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
                     "wordcloud-tab"
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
                     "Similar experiments",
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
             div(class="sidebar-item",
                 bigdash::sidebarItem(
                     "Cell profiling",
                     "cell-tab"
                 )
             )
        )
        
        bigdash::bigPage(
            header,
            sidebar = sidebar,
            navbar = bigdash::navbar(
                tags$img(
                    src = "assets/img/bigomics.png",
                    width = "110",
                ),
                shiny::div(shiny::textOutput("current_section"), class='current-section'),
                shiny::div(shiny::textOutput("current_dataset"), class='current-dataset'),
                bigdash::navbarDropdown(
                    "Support",
                    bigdash::navbarDropdownItem(
                        "Community Forum",
                        link = "https://groups.google.com/d/forum/omicsplayground",
                        target = "_blank"
                    ),
                    bigdash::navbarDropdownItem(
                        "Github issues",
                        link = "https://github.com/bigomics/omicsplayground/issues",
                        target = "_blank"
                    ),
                    bigdash::navbarDropdownItem(
                        "www.bigomics.ch",
                        link = "http://bigomics.ch",
                        target = "_blank"
                    )
                ),
                bigdash::navbarDropdown(
                    "Tutorials",
                    bigdash::navbarDropdownItem(
                        "Documentation",
                        link = "https://omicsplayground.readthedocs.io",
                        target = "_blank"
                    ),
                    bigdash::navbarDropdownItem(
                        "Video tutorials",
                        link = "https://www.youtube.com/watch?v=_Q2LJmb2ihU&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-",
                        target = "_blank"
                    ),
                    bigdash::navbarDropdownItem(
                        "Case studies",
                        link = "https://bigomics.ch/category/case-study/",
                        target = "_blank"
                    )
                ),
                bigdash::navbarDropdown(
                    ##"User",
                    ##shiny::div(class='label label-info current-user',id='authentication-user'),
                    shiny::textOutput("current_user", inline = TRUE),
                    bigdash::navbarDropdownTab(
                        "Settings",
                        "userSettings"
                    ),
                    upgrade.tab,
                    logout.tab
                )
            ),
            settings = bigdash::settings(
                "Settings"
            ),
            bigdash::sidebarHelp(
                bigdash::sidebarTabHelp(
                    "welcome-tab",
                    "BigOmics Playground",
                    "is your self-service bioinformatics platform for interactive analysis,
                    visualization and interpretation of transcriptomics and proteomics data.
                    Perform complex data analysis and visualization easily without coding,
                    and significantly reduce the time-to-discovery."
                ),
                bigdash::sidebarTabHelp(
                    "load-tab",
                    "Load dataset",
                    "This panel shows the available datasets within the platform. These data sets
                     have been pre-computed and are ready to be used. Select a
                     dataset in the table and load the data set by clicking the 'load' button."
                ),
                bigdash::sidebarTabHelp(
                    "upload-tab",
                    "Upload data",
                    "Here you can upload your own transcriptomics and proteomics data into
                     the platform and perform computations for the Playground."
                ),
                bigdash::sidebarTabHelp(
                    "dataview-tab",
                    "DataView",
                    "Information and descriptive statistics to quickly lookup a gene,
                     check your experiment QC, view the raw data, sample or contrast tables."
                ),
                bigdash::sidebarTabHelp(
                    "clustersamples-tab",
                    "Clustering Analysis",
                    "Discover clusters of similar genes or samples using unsupervised
                    machine learning."
                ),
                bigdash::sidebarTabHelp(
                    "wgcna-tab",
                    "Weighted Correlation",
                    "Weighted correlation network analysis (WGCNA) is a gene-level cluster analysis
                    method based on pairwise correlations between genes. It allows one to define modules (clusters),
                    intramodular hubs, and network nodes with regard to module membership, to study the relationships
                    between co-expression modules."
                ),
                bigdash::sidebarTabHelp(
                    "diffexpr-tab",
                    "Expression Analysis",
                    "Compare expression between
                    two conditions. Determine which genes are significantly downregulated or overexpressed in one of the groups."
                ),
                bigdash::sidebarTabHelp(
                    "corr-tab",
                    "Correlation Analysis",
                    "Compute the correlation between genes and find coregulated modules."
                ),
                bigdash::sidebarTabHelp(
                    "enrich-tab",
                    "Geneset Enrichment",
                    "Perform differential expression analysis on a geneset level,
                    also called geneset enrichment analysis."
                ),
                bigdash::sidebarTabHelp(
                    "pathway-tab",
                    "Functional Analysis",
                    "Perform specialized functional analysis
                    to understand biological functions including GO, KEGG, and drug connectivity mapping."
                ),
                bigdash::sidebarTabHelp(
                    "wordcloud-tab",
                    "Wordcloud",
                    "WordCloud analysis or 'keyword enrichment' analysis computes the
                    enrichment of keywords for the contrasts. The set of words frequently appearing in the top ranked
                    gene sets form an unbiased description of the contrast."
                ),
                bigdash::sidebarTabHelp(
                    "drug-tab",
                    "Drug Connectivity",
                    "Perform drug connectivity analysis
                    to see if certain drug activity or drug sensitivity signatures matches your experimental signatures.
                    Matching drug signatures to your experiments may elicudate biological functions through
                    mechanism-of-action (MOA) and known drug molecular targets."
                ),
                bigdash::sidebarTabHelp(
                    "isect-tab",
                    "Compare Signatures",
                    "Find genes that are commonly up/down regulated
                    between two or more signatures. Compute similarity between contrasts."
                ),
                bigdash::sidebarTabHelp(
                    "sig-tab",
                    "Signature Analysis",
                    "Users can test their gene signature by
                    calculating an enrichment score. Upload your own gene list, or select
                    a contrast which then takes the top differentially expressed genes as
                    signature."
                ),
                bigdash::sidebarTabHelp(
                    "bio-tab",
                    "Biomarker Board",
                    "Select biomarkers that can be used for
                    classification or prediction purposes. The phenotype of interest can
                    be multiple categories (classes) or patient survival data."
                ),
                bigdash::sidebarTabHelp(
                    "cmap-tab",
                    "Similar Experiments",
                    "Find similar experiments by correlating their signatures.
                    The main goal is to identify experiments showing similar signatures and find genes
                    that are commonly up/down regulated between experiments."
                ),
                bigdash::sidebarTabHelp(
                    "comp-tab",
                    "Compare Datasets",
                    "Compare expression and signatures between two datasets,
                    from similar experiments or from different datatypes, e.g. transcriptomics and proteomics."
                ),
                bigdash::sidebarTabHelp(
                    "tcga-tab",
                    "TCGA Analysis",
                    "Correlate your signature with the survival in cancer patients from the TCGA database. Warning: EXPERIMENTAL."
                ),
                bigdash::sidebarTabHelp(
                    "cell-tab",
                    "Single-Cell Profiling",
                    "Visualize the distribution of (inferred)
                    immune cell types, expressed genes and pathway activation."
                )
            ),
            bigdash::bigTabs(
                bigdash::bigTabItem(
                    "welcome-tab",
                    WelcomeBoardInputs("welcome"),
                    WelcomeBoardUI("welcome", ENABLED['upload'])
                ),
                bigdash::bigTabItem(
                    "load-tab",
                    LoadingInputs("load"),
                    LoadingUI("load")
                ),
                bigdash::bigTabItem(
                    "upload-tab",
                    # UploadInputs("upload"),
                    UploadUI("upload")
                ),
                bigdash::bigTabItem(
                    "dataview-tab",
                    DataViewInputs("dataview"),
                    DataViewUI("dataview")
                ),
                bigdash::bigTabItem(
                    "clustersamples-tab",
                    ClusteringInputs("clustersamples"),
                    ClusteringUI("clustersamples")
                ),
                bigdash::bigTabItem(
                    "clusterfeatures-tab",
                    FeatureMapInputs("clusterfeatures"),
                    FeatureMapUI("clusterfeatures")
                ),
                bigdash::bigTabItem(
                    "wgcna-tab",
                    WgcnaInputs("wgcna"),
                    WgcnaUI("wgcna")
                ),
                bigdash::bigTabItem(
                    "diffexpr-tab",
                    ExpressionInputs("diffexpr"),
                    ExpressionUI("diffexpr")
                ),
                bigdash::bigTabItem(
                    "corr-tab",
                    CorrelationInputs("corr"),
                    CorrelationUI("corr")
                ),
                bigdash::bigTabItem(
                    "enrich-tab",
                    EnrichmentInputs("enrich"),
                    EnrichmentUI("enrich")
                ),
                bigdash::bigTabItem(
                    "pathway-tab",
                    FunctionalInputs("pathway"),
                    FunctionalUI("pathway")
                ),
                bigdash::bigTabItem(
                    "wordcloud-tab",
                    WordCloudInputs("wordcloud"),
                    WordCloudUI("wordcloud")
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
                    SingleCellInputs("cell"),
                    SingleCellUI("cell")
                ),
                bigdash::bigTabItem(
                    "userSettings",
                    UserInputs("user"),
                    UserUI("user")
                )
            ),
            tagList(footer),
            # This JS script is sourced at this point because if sourced before creating
            # the UI elements, the snippet to prevent `dropdown-menu` to propagate clicks
            # will not be applied
            tags$script(src = "dropdown-helper.js")
        )
    }

    ui <- createUI()

    return(ui)
}

