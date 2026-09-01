##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Full sidebar menu: module group -> named vector of boards (id = title).
#' Shared with the admin panel, which offers these as the basic-menu choices.
opg_menu_tree <- function() {
  list(
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
}


#' Board ids kept in BASIC mode, in full-menu order. The selection is chosen
#' in Admin panel > Basic menu and persisted per deploy in
#' etc/BASIC_MENU-<HOSTNAME> (read in global.R). There is only one sidebar
#' menu tree now -- opg_server.R's tab_control() intersects this with the set
#' it already computes from dataset/content availability and hands the result
#' to bigdash::bigdash.filterTabs(), so this only needs to express the admin's
#' picks, not what is actually available right now.
opg_basic_menu_boards <- function(menu_tree = opg_menu_tree(), boards = opt$BASIC_MENU) {
  all_boards <- names(unlist(unname(menu_tree)))
  boards <- intersect(all_boards, boards) ## drops unknown/disabled boards, keeps full-menu order
  if (!length(boards)) boards <- intersect(all_boards, BASIC_MENU_DEFAULT)
  if (!length(boards)) return(character(0)) ## paste0(character(0), "-tab") == "-tab", not character(0)
  paste0(boards, "-tab")
}


opg_ui <- function(id) {

  message("\n======================================================")
  message("======================= UI ===========================")
  message("======================================================\n")

  #-------------------------------------------------------
  ## Build USERMENU
  #-------------------------------------------------------
VERSION <- scan(file.path(OPG, "VERSION"), character())[1]

  createUI <- function(menu_tree) {
    
    version <- scan(file.path(OPG, "VERSION"), character())[1]
    ##id <- "maintabs"
    
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
    } else if (opt$AUTHENTICATION %in% c("shinyproxy-sso", "shinyproxy-sso-admin")) {
      ## Upstream-header auth (e.g. ShinyProxy + SAML). Hit /logout so
      ## ShinyProxy clears its session and triggers the IdP SLO via the
      ## configured saml.logout-url, then bounce back to /login.
      logout.tab <- bigdash::navbarDropdownItem(
        "Logout",
        link = "/logout"
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

    menu_tree <- menu_tree[MODULES_ENABLED]
    ENABLED <<- array(BOARDS %in% sapply(menu_tree, function(m) names(m)), dimnames = list(BOARDS))

    createMenu <- function(tree) {
      sidebar_item <- function(title, name) {
        #div(class = "sidebar-item", bigdash::sidebarItem(title, paste0(name, "-tab")))
        bigdash::sidebarItem(title, paste0(name, "-tab"))      
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
        ## promote_single: with one shared tree now driving both full and
        ## basic mode via bigdash.filterTabs(), a group filtered down to one
        ## visible board (e.g. Basic menu keeping only "Samples" out of
        ## Clustering) renders as a flat top-level item instead of a
        ## one-item group -- matching the old, separately-built flat basic
        ## menu's look without needing a second tree.
        bigdash::sidebarMenu(title, !!!ee, promote_single = TRUE)
      }
      ## a config can filter every board out (see MODULES_ENABLED below); an
      ## empty menu is survivable, `for (i in 1:0)` is not
      if (!length(tree)) {
        return(HTML(""))
      }
      menu <- list()
      for (i in 1:length(tree)) {
        tab.names <- names(tree[[i]])
        tab.titles <- tree[[i]]
        menu.id <- names(tree)[i]
        if (length(tab.names) == 0) {} else if (length(tab.names) == 1 && tolower(tab.names) == tolower(menu.id)) {
          menu[[menu.id]] <- sidebar_item(tab.titles, tab.names)
        } else {
          menu[[menu.id]] <- sidebar_menu_with_items(tree[[i]], menu.id)
        }
      }
      menu <- lapply(menu, as.character)
      HTML(unlist(menu))
    }

    info("[opg_ui] creating sidebar menu")

    ## One tree now for both full and basic mode -- opg_server.R's
    ## tab_control() drives which items are visible via
    ## bigdash::bigdash.filterTabs(), which hides/shows the matching
    ## sidebar-item/sidebarMenuItem elements (and auto-promotes a group left
    ## with a single visible item, see bigdash's refreshMenuPromotion()).
    sidebar <- bigdash::sidebar(
      "Menu",
      createMenu(menu_tree)
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
      id = id,  ## default was 'app'
      shiny.i18n::usei18n(i18n),
      ## shiny.i18n's subscribe() hands jQuery's Event object straight to Shiny's
      ## callback, which expects a boolean, so every update_lang() logs
      ## "Unexpected input value mode: '[object Object]'". Must come after
      ## usei18n() so the binding is already registered.
      shiny::tags$script(shiny::HTML(
        "Shiny.inputBindings.bindingNames['shiny.shinyi18n'].binding.subscribe =
           function (el, callback) { $(el).on('change.shinyi18n', function () { callback(false); }); };"
      )),
      # header,
      title = "Omics Playground",
      theme = big_theme2,
      navbar = navbar,
      sidebar = sidebar,
      settings = bigdash::settings(
        "Settings"
      ),
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
        ## One shell per tab: that tab's own inputs -- which settings.js moves
        ## into the settings drawer at boot -- plus a spinner. bigTabsLazy() in
        ## the server fills in the board itself on first visit.
        ##
        ## Taken from each group's module_ui() rather than spelled out here, so
        ## a tab's shell and the module_lazy() entry that fills it stay in one
        ## file and cannot drift apart. DataView is the one tab with no group
        ## registry of its own.
        ##
        ## lock_advanced() greys a locked board's settings accordions for BASIC
        ## users. The accordions live in this shell (each board's *Inputs()),
        ## not in the lazily-loaded body, so it belongs here, at boot, rather
        ## than wherever bigTabsLazy() later fills the tab in.
        lock_advanced(bigdash::bigTabItem("dataview-tab", DataViewInputs("dataview"), create_loader("dataview-loader"))),
        lock_advanced(MODULE.clustering$module_ui()),
        lock_advanced(MODULE.expression$module_ui()),
        lock_advanced(MODULE.enrichment$module_ui()),
        lock_advanced(MODULE.compare$module_ui()),
        lock_advanced(MODULE.systems$module_ui()),
        lock_advanced(MODULE.multiomics$module_ui()),
        lock_advanced(MODULE.wgcna$module_ui())
      )
    ) ## end of bigPage
  }


  full_menu_tree <- opg_menu_tree()

  info("[opg_ui] >>> creating UI")
  ui <- createUI(full_menu_tree)
  info("[opg_ui] <<< finished UI!")

  return(ui)
}
