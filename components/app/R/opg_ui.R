##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

opg_ui <- function() {

  message("\n======================================================")
  message("======================= UI ===========================")
  message("======================================================\n")

  #-------------------------------------------------------
  ## Build USERMENU
  #-------------------------------------------------------
  VERSION <- scan(file.path(OPG, "VERSION"), character())[1]

  ## upgrade.tab <- NULL
  ## if (opt$AUTHENTICATION == "firebase") {
  ##   upgrade.tab <- bigdash::navbarDropdownItem(
  ##     "Upgrade",
  ##     onClick = "show_plans()"
  ##   )
  ## }
    
  createUI <- function(menu_tree) {
    
    version <- scan(file.path(OPG, "VERSION"), character())[1]
    id <- "maintabs"
    
    logout.tab <- bigdash::navbarDropdownItem(
      "Logout",
      onClick = "logoutInApp(); setTimeout(() => window.location.reload(), 200);"
    )

    if (opt$AUTHENTICATION == "shinyproxy") {
      ## For ShinyProxy we need to redirect to /logout for clean session
      ## logout. Then we need a redirect to the /login page.
      logout.tab <- bigdash::navbarDropdownItem(
        "Exit",
        onClick = "shinyproxy_logout();",
        link = "/login"
      )
    } else if (opt$AUTHENTICATION == "apache-cookie") {
      ## For apache SSO we need to redirect to /mellon/logout for SSO logout
      logout.tab <- bigdash::navbarDropdownItem(
        "Logout",
        link = paste0(opt$APACHE_COOKIE_PATH, "mellon/logout?ReturnTo=#")
      )
    }

    ## filter disabled modules
    #ENABLED["welcome"] <<- TRUE
    #ENABLED["load"] <<- TRUE

    dbg("names(menu_tree) = ", names(menu_tree))
    dbg("names.ENABLED = ", names(ENABLED))
    menu_tree <- menu_tree[MODULES_ENABLED]
    ## menu_tree <- lapply(menu_tree, function(m) m[which(ENABLED[names(m)])])
    ENABLED <<- array(BOARDS %in% sapply(menu_tree, function(m) names(m)), dimnames = list(BOARDS))

    createMenu <- function(tree) {
      sidebar_item <- function(title, name) {
        div(class = "sidebar-item", bigdash::sidebarItem(title, paste0(name, "-tab")))
      }
      sidebar_menu_item <- function(title, name) {
        bigdash::sidebarMenuItem(title, paste0(name, "-tab"))
      }
      sidebar_menu_with_items <- function(tabs, title) {
        ee <- list()
        for (i in 1:length(tabs)) {
          tab.name <- names(tabs)[i]
          tab.title <- tabs[i]
          ee[[i]] <- sidebar_menu_item(tab.title, tab.name)
        }
        bigdash::sidebarMenu(title, !!!ee)
      }
      menu <- list()
      for (i in 1:length(tree)) {
        tab.names <- names(tree[[i]])
        tab.titles <- tree[[i]]
        menu.id <- names(tree)[i]
        if (length(tab.names) == 0) {} else if (length(tab.names) == 1) {
          menu[[menu.id]] <- sidebar_item(tab.titles, tab.names)
        } else {
          menu[[menu.id]] <- sidebar_menu_with_items(tree[[i]], menu.id)
        }
      }
      menu <- lapply(menu, as.character)
      HTML(unlist(menu))
    }

    basic_menu_tree <- list(
      "DataView"                = c(dataview       = "DataView"),
      "Cluster Samples"         = c(clustersamples = "Cluster Samples"),
      "Differential expression" = c(diffexpr       = "Differential expression"),
      "Geneset Enrichment"      = c(enrich         = "Geneset Enrichment"),
      "Pathway analysis"        = c(pathway        = "Pathway analysis")
    )

    initial_is_full <- (opt$USER_LEVEL != "BASIC")
    info("[opg_ui] creating sidebar menu")
    info("[opg_ui] initial_is_full = ", initial_is_full)
    
    sidebar <- bigdash::sidebar(
      "Menu",
      ## div(
      ##   id = "menu-mode-switch",
      ##   shinyWidgets::switchInput(
      ##     "menu_mode_toggle",
      ##     label    = NULL,
      ##     value    = initial_is_full,
      ##     onLabel  = "Advanced",
      ##     offLabel = "Basic",
      ##     size     = "mini"
      ##   )
      ## ),
      div(
        id = "menu-full",
        class = "nodisp", style = "diplay: none;",
        createMenu(menu_tree)
      ),
      div(
        id = "menu-basic",
        class = "nodisp", style = "diplay: none;",        
        createMenu(basic_menu_tree)
      )
    )
    
    big_theme2 <- bigdash::big_theme()
    big_theme2 <- bslib::bs_add_variables(big_theme2,
      "grid-breakpoints" = "map-merge($grid-breakpoints, ('xxxl': 2400px))",
      .where = "declarations"
    )

    ## ------------------------- bigPage ----------------------------------
    make_sidebarHelp <- function(...) {
      do.call(bigdash::sidebarHelp, rlang::list2(...))
    }

    ## empty navbar
    navbar <- div(
      ## THIS IS SO WEIRD. if we remove/comment out the
      ## prettySwitch, the header of all plotModules f*ck
      ## up... (IK). HELP!!! we do not need this button...      
      style = "visibility: hidden; display: none;",
      shinyWidgets::prettySwitch("I_AM_WEIRD_BUTTON", "remove me")
    )
    
    bigdash::bigPage(
      shiny.i18n::usei18n(i18n),
      # header,
      title = "Omics Playground 4",
      theme = big_theme2,
      navbar = navbar,
      sidebar = sidebar,
      settings = bigdash::settings(
        "Settings"
      ),
      ## bigdash::sidebarHelp(
      make_sidebarHelp(
        bigdash::sidebarTabHelp(
          "upload-tab",
          "Upload new",
          "Here you can upload your own transcriptomics and proteomics data into
                    the platform and perform computations for the Playground."
        ),
        bigdash::sidebarTabHelp(
          "dataview-tab",
          "DataView",
          tspan("Information and descriptive statistics to quickly lookup a gene,
                    check your experiment QC, view the raw data, sample or contrast tables.")
        ),
        bigdash::sidebarTabHelp(
          "clustersamples-tab",
          "Clustering Analysis",
          tspan("Discover clusters of similar genes or samples using unsupervised
                    machine learning.")
        ),
        bigdash::sidebarTabHelp(
          "wgcna-tab",
          "Weighted Correlation",
          tspan("Weighted correlation network analysis (WGCNA) is a gene-level cluster
                    analysis method based on pairwise correlations between genes. It
                    allows one to define modules (clusters), intramodular hubs, and
                    network nodes with regard to module membership, to study the
                    relationships between co-expression modules.")
        ),
        bigdash::sidebarTabHelp(
          "mofa-tab",
          "MOFA",
          tspan("Multi-omics Factor Analysis (MOFA) is a multi-omics
                  integration method based on multi-omcis factor analysis.")
        ),
        bigdash::sidebarTabHelp(
          "mgsea-tab",
          "multiGSEA",
          tspan("multiGSEA perform multi-omics integration on gene set level.")
        ),
        bigdash::sidebarTabHelp(
          "snf-tab",
          "SNF",
          tspan("SNF clustering")
        ),
        bigdash::sidebarTabHelp(
          "pcsf-tab",
          "PCSF Network Analysis",
          "PCSF performs fast and user-friendly network analysis using
                    interaction networks as a template, it determines high-confidence
                    subnetworks relevant to the data, which potentially leads to
                    predictions of functional units."
        ),
        bigdash::sidebarTabHelp(
          "diffexpr-tab",
          "Expression Analysis",
          tspan("Compare expression between two conditions. Determine which genes are
                    significantly downregulated or overexpressed in one of the groups.")
        ),
        bigdash::sidebarTabHelp(
          "corr-tab",
          "Correlation Analysis",
          tspan("Compute the correlation between genes and find coregulated modules.")
        ),
        bigdash::sidebarTabHelp(
          "enrich-tab",
          "Geneset Enrichment",
          tspan("Perform differential expression analysis on a geneset level,
                    also called geneset enrichment analysis.")
        ),
        bigdash::sidebarTabHelp(
          "pathway-tab",
          "Pathway Analysis",
          "Perform functional analysis to understand biological pathways using WikiPathways, Reactome and Gene Ontology."
        ),
        bigdash::sidebarTabHelp(
          "wordcloud-tab",
          "Wordcloud",
          tspan("WordCloud analysis or 'keyword enrichment' analysis computes the
                    enrichment of keywords for the contrasts. The set of words frequently appearing in the top ranked
                    genesets form an unbiased description of the contrast.")
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
          tspan("Find genes that are commonly up/down regulated
                    between two or more signatures. Compute similarity between contrasts.")
        ),
        bigdash::sidebarTabHelp(
          "sig-tab",
          "Signature Analysis",
          tspan("Users can test their gene signature by
                    calculating an enrichment score. Upload your own gene list, or select
                    a contrast which then takes the top differentially expressed genes as
                    signature.")
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
          tspan("Find similar experiments by correlating their signatures.
                 The main goal is to identify experiments showing similar signatures and find 
                 genes that are commonly up/down regulated between experiments.")
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
          tspan("Visualize the distribution of (inferred)
                    immune cell types, expressed genes and pathway activation.")
        ),
        bigdash::sidebarTabHelp(
          "consensus-tab",
          "Consensus WGCNA",
          tspan("Consensus analysis using the WGCNA framework")
        ),
        bigdash::sidebarTabHelp(
          "preservation-tab",
          "Preservation WGCNA",
          tspan("Preservation analysis using the WGCNA framework")
        ),
        !!!MODULE.multiomics$module_help() ### HELP!!! DOES NOT WORK!!!
      ),
      bigdash::bigTabs(
        ## bigdash::bigTabItem(
        ##   "welcome-tab",
        ##   WelcomeBoardInputs("welcome"),
        ##   WelcomeBoardUI("welcome")
        ## ),
        ## bigdash::bigTabItem(
        ##   "load-tab",
        ##   # LoadingInputs("load")
        ##   LoadingUI("load")
        ## ),
        ## bigdash::bigTabItem(
        ##   "upload-tab",
        ##   UploadUI("upload")
        ## ),
        ## bigdash::bigTabItem(
        ##   "userprofile-tab",
        ##   UserProfileUI("user_profile")
        ## ),
        ## bigdash::bigTabItem(
        ##   "usersettings-tab",
        ##   AppSettingsInputs("app_settings"),
        ##   AppSettingsUI("app_settings")
        ## ),
        ## if (isTRUE(opt$ENABLE_ADMIN)) {
        ##   bigdash::bigTabItem(
        ##     "admin-tab",
        ##     AdminPanelUI("admin_panel")
        ##   )
        ## }
        ## bigdash::bigTabItem(
        ##   "sharing-tab",
        ##   SharedDatasetsUI("load")
        ## )
      )
      ## UploadUI("upload")
    ) ## end of bigPage
  }


  full_menu_tree <- list(
    "DataView" = c(
      dataview = "DataView"
    ),
    "Clustering" = c(
      clustersamples = "Samples",
      clusterfeatures = "Features"
    ),
    "Expression" = c(
      diffexpr = "Differential expression",
      timeseries = "TimeSeries", ## here???
      corr = "Correlation analysis",
      bio = "Find biomarkers"
    ),
    "GeneSets" = c(
      enrich = "Geneset Enrichment",
      sig = "Test geneset",
      pathway = "Pathway analysis",
      wordcloud = "Word cloud"
    ),
    "Compare" = c(
      isect = "Compare signatures",
      comp = "Compare datasets",
      cmap = "Similar experiments"
    ),
    "SystemsBio" = c(
      drug = "Drug connectivity",
      cell = "Cell profiling",
      pcsf = "PCSF",
      tcga = "TCGA survival (beta)"
    ),
    "MultiOmics" = MODULE.multiomics$module_menu(),
    "WGCNA" = MODULE.wgcna$module_menu()
  )

  info("[opg_ui] >>> creating UI")
  ui <- createUI(full_menu_tree)
  info("[opg_ui] <<< finished UI!")

  return(ui)
}
