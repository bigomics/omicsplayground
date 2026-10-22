##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
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
    "WGCNA" = MODULE.wgcna$module_menu(),
    "Epigenomics" = MODULE.epigenomics$module_menu()
  )
}


#' Bare board ids (e.g. "pcsf", "dataview") a menu_tree allows, flattened
#' across all its groups regardless of group membership. Used to filter UI
#' shells/server registrations down to exactly the boards a menu_tree lists,
#' even when only some boards of an otherwise-included group are wanted
#' (see opg_multiomics_menu_tree()).
opg_menu_boards <- function(menu_tree) {
  unique(unlist(lapply(menu_tree, names)))
}

#' Split a menu_tree's "MultiOmics" group into flat, single-board top-level
#' entries, leaving every other group untouched.
#'
#' The single "app" OPG instance's sidebar DOM is built once per session by
#' opg_ui()'s createMenu() -- opg_server.R's tab_control() only shows/hides
#' those existing nodes at runtime (bigdash.filterTabs()), it cannot turn a
#' collapsible group into flat items or back. So "flat vs grouped" has to be
#' decided here, before the tree ever reaches opg_ui() -- not left to
#' whichever menu_tree opg_server.R happens to be reactively filtering with.
#' createMenu() renders any tree entry with exactly one board as a flat
#' top-level sidebar item (like DataView), so splitting MultiOmics's boards
#' (MOFA, multiGSEA, SNF, ...) into one-board entries promotes them out of
#' the collapsible "MultiOmics" submenu. One split entry keeps the
#' "MultiOmics" key (whichever board holds it doesn't matter -- the key is
#' never shown) purely so tab_control()'s "MultiOmics" %in% allowed_groups()
#' gate still recognizes the group as present.
opg_promote_multiomics <- function(menu_tree) {
  mo_boards <- menu_tree[["MultiOmics"]]
  if (is.null(mo_boards) || length(mo_boards) <= 1) {
    return(menu_tree)
  }
  mo_entries <- lapply(seq_along(mo_boards), function(i) mo_boards[i])
  names(mo_entries) <- c("MultiOmics", names(mo_boards)[-1])

  i <- match("MultiOmics", names(menu_tree))
  c(
    menu_tree[seq_len(i - 1)],
    mo_entries,
    menu_tree[seq_len(length(menu_tree) - i) + i]
  )
}

#' MultiOmics-only menu tree for the single "app" OPG instance's MultiOmics
#' view (server.R's opg_view()): the full MultiOmics group, plus DataView and
#' all Expression tabs, plus single-board additions from groups that
#' otherwise don't belong here (PCSF from SystemsBio, multiWGCNA from WGCNA).
#' MultiOmics's own boards are promoted flat, same as the main sidebar --
#' see opg_promote_multiomics().
opg_multiomics_menu_tree <- function(menu_tree = opg_menu_tree()) {
  opg_promote_multiomics(list(
    DataView   = menu_tree[["DataView"]],
    Expression = menu_tree[["Expression"]],
    SystemsBio = menu_tree[["SystemsBio"]]["pcsf"],
    MultiOmics = menu_tree[["MultiOmics"]],
    WGCNA      = menu_tree[["WGCNA"]]["mwgcna"]
  ))
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


#' Namespace helper matching bigdash's own `BIGDASH_DEFAULT_ID` ("app")
#' convention: `id = "app"` -- the pre-existing, single-instance OPG -- stays
#' completely unscoped (bare ids, byte-identical markup/inputs to before this
#' module conversion); any other `id` is a real nested `bigdash::bigPage()`
#' instance and gets every id `"<id>-"`-prefixed, matching what
#' `bigdash:::scoped_id()` already does for the container ids it generates
#' and what `shiny::moduleServer(id, ...)`'s own `session$ns()` produces for
#' the boards nested inside. Kept local (rather than reaching into bigdash's
#' unexported `scoped_id()`) so this file has no dependency on bigdash
#' internals.
opg_ns <- function(id) {
  function(suffix) {
    if (is.null(id) || identical(id, "app")) suffix else paste0(id, "-", suffix)
  }
}

opg_ui <- function(id, menu_tree = opg_menu_tree()) {

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
    ns <- opg_ns(id)

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

    ## Fully-namespaced tab ids this instance's menu_tree allows -- narrows
    ## a group's module_ui() shells down to just the boards menu_tree lists
    ## for that group (e.g. only "pcsf-tab" out of all of SystemsBio's),
    ## not just whole groups. See filter_boards() below and opg_server.R's
    ## matching filter_lazy() for the server-side half.
    allowed_tab_ids <- ns(paste0(opg_menu_boards(menu_tree), "-tab"))
    filter_boards <- function(tabs) {
      if (inherits(tabs, "shiny.tag")) tabs <- list(tabs)
      Filter(function(tag) {
        nm <- tag$attribs[["data-name"]]
        is.null(nm) || nm %in% allowed_tab_ids
      }, tabs)
    }

    createMenu <- function(tree) {
      sidebar_item <- function(title, name) {
        #div(class = "sidebar-item", bigdash::sidebarItem(title, ns(paste0(name, "-tab"))))
        bigdash::sidebarItem(title, ns(paste0(name, "-tab")))
      }
      sidebar_menu_item <- function(title, name) {
        bigdash::sidebarMenuItem(title, ns(paste0(name, "-tab")))
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
        if (length(tab.names) == 0) {} else if (length(tab.names) == 1) {
          ## A group with exactly one board renders as a flat top-level
          ## item, not a one-item collapsible group -- matching what
          ## bigdash's client-side promote_single gives a group *runtime*
          ## filtered down to one, but built in from the start server-side
          ## here, since promote_single only reacts to a change and leaves
          ## an empty "<group> >" header behind for a group that's single-
          ## item from the very first render (e.g. Epigenomics, or a
          ## menu_tree -- opg_multiomics_menu_tree() -- that only wants one
          ## board out of an otherwise multi-board group like SystemsBio).
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
      createMenu(menu_tree),
      id = id
    )
    
    big_theme2 <- bigdash::big_theme()
    big_theme2 <- bslib::bs_add_variables(big_theme2,
      "grid-breakpoints" = "map-merge($grid-breakpoints, ('xxxl': 2400px))",
      .where = "declarations"
    )

    ## ------------------------- bigPage ----------------------------------
    make_sidebarHelp <- function(...) {
      do.call(bigdash::sidebarHelp, rlang::list2(..., id = id))
    }

    ## empty navbar
    navbar <- div(
      ## THIS IS SO WEIRD. if we remove/comment out the
      ## prettySwitch, the header of all plotModules f*ck
      ## up... (IK). HELP!!! we do not need this button...
      ## ns()-wrapped: navbar sits outside bigPage()'s own data-bigdash-id
      ## wrapper, so a bare id here would collide with a second instance's.
      style = "visibility: hidden; display: none;",
      shinyWidgets::prettySwitch(ns("I_AM_WEIRD_BUTTON"), "remove me")
    )
    
    bigdash::bigPage(
      id = id,  ## default was 'app'
      ## shiny.i18n::usei18n() renders its own state div with a fixed id
      ## ("i18n-state"), unaware of bigdash's per-instance namespacing --
      ## calling it once per opg_ui() instance would duplicate that id (the
      ## second instance's div becomes unreachable to shiny.i18n's own JS,
      ## since getElementById()/jQuery only ever find the first match). Only
      ## the primary instance sets it up; language is app-wide anyway, not
      ## per-OPG-instance.
      if (identical(id, "app")) {
        list(
          shiny.i18n::usei18n(i18n),
          ## shiny.i18n's subscribe() hands jQuery's Event object straight to Shiny's
          ## callback, which expects a boolean, so every update_lang() logs
          ## "Unexpected input value mode: '[object Object]'". Must come after
          ## usei18n() so the binding is already registered.
          shiny::tags$script(shiny::HTML(
            "Shiny.inputBindings.bindingNames['shiny.shinyi18n'].binding.subscribe =
               function (el, callback) { $(el).on('change.shinyi18n', function () { callback(false); }); };"
          ))
        )
      },
      # header,
      title = "Omics Playground",
      theme = big_theme2,
      navbar = navbar,
      sidebar = sidebar,
      settings = bigdash::settings(
        "Settings",
        id = id
      ),
      make_sidebarHelp(
        ## "upload-tab" is the top-level Upload nav panel, outside this
        ## bigPage()'s own bigTabs() -- left un-namespaced, same as before.
        bigdash::sidebarTabHelp(
          "upload-tab",
          "Upload new",
          "Here you can upload your own transcriptomics and proteomics data into
                    the platform and perform computations for the Playground."
        ),
        bigdash::sidebarTabHelp(
          ns("dataview-tab"),
          "DataView",
          tspan("Information and descriptive statistics to quickly lookup a gene,
                    check your experiment QC, view the raw data, sample or contrast tables.")
        ),
        bigdash::sidebarTabHelp(
          ns("clustersamples-tab"),
          "Clustering Analysis",
          tspan("Discover clusters of similar genes or samples using unsupervised
                    machine learning.")
        ),
        bigdash::sidebarTabHelp(
          ns("wgcna-tab"),
          "Weighted Correlation",
          tspan("Weighted correlation network analysis (WGCNA) is a gene-level cluster
                    analysis method based on pairwise correlations between genes. It
                    allows one to define modules (clusters), intramodular hubs, and
                    network nodes with regard to module membership, to study the
                    relationships between co-expression modules.")
        ),
        bigdash::sidebarTabHelp(
          ns("mofa-tab"),
          "MOFA",
          tspan("Multi-omics Factor Analysis (MOFA) is a multi-omics
                  integration method based on multi-omcis factor analysis.")
        ),
        bigdash::sidebarTabHelp(
          ns("mgsea-tab"),
          "multiGSEA",
          tspan("multiGSEA perform multi-omics integration on gene set level.")
        ),
        bigdash::sidebarTabHelp(
          ns("snf-tab"),
          "SNF",
          tspan("SNF clustering")
        ),
        bigdash::sidebarTabHelp(
          ns("pcsf-tab"),
          "PCSF Network Analysis",
          "PCSF performs fast and user-friendly network analysis using
                    interaction networks as a template, it determines high-confidence
                    subnetworks relevant to the data, which potentially leads to
                    predictions of functional units."
        ),
        bigdash::sidebarTabHelp(
          ns("diffexpr-tab"),
          "Expression Analysis",
          tspan("Compare expression between two conditions. Determine which genes are
                    significantly downregulated or overexpressed in one of the groups.")
        ),
        bigdash::sidebarTabHelp(
          ns("corr-tab"),
          "Correlation Analysis",
          tspan("Compute the correlation between genes and find coregulated modules.")
        ),
        bigdash::sidebarTabHelp(
          ns("enrich-tab"),
          "Geneset Enrichment",
          tspan("Perform differential expression analysis on a geneset level,
                    also called geneset enrichment analysis.")
        ),
        bigdash::sidebarTabHelp(
          ns("pathway-tab"),
          "Pathway Analysis",
          "Perform functional analysis to understand biological pathways using WikiPathways, Reactome and Gene Ontology."
        ),
        bigdash::sidebarTabHelp(
          ns("wordcloud-tab"),
          "Wordcloud",
          tspan("WordCloud analysis or 'keyword enrichment' analysis computes the
                    enrichment of keywords for the contrasts. The set of words frequently appearing in the top ranked
                    genesets form an unbiased description of the contrast.")
        ),
        bigdash::sidebarTabHelp(
          ns("drug-tab"),
          "Drug Connectivity",
          "Perform drug connectivity analysis
                    to see if certain drug activity or drug sensitivity signatures matches your experimental signatures.
                    Matching drug signatures to your experiments may elicudate biological functions through
                    mechanism-of-action (MOA) and known drug molecular targets."
        ),
        bigdash::sidebarTabHelp(
          ns("isect-tab"),
          "Compare Signatures",
          tspan("Find genes that are commonly up/down regulated
                    between two or more signatures. Compute similarity between contrasts.")
        ),
        bigdash::sidebarTabHelp(
          ns("sig-tab"),
          "Signature Analysis",
          tspan("Users can test their gene signature by
                    calculating an enrichment score. Upload your own gene list, or select
                    a contrast which then takes the top differentially expressed genes as
                    signature.")
        ),
        bigdash::sidebarTabHelp(
          ns("bio-tab"),
          "Biomarker Board",
          "Select biomarkers that can be used for
                    classification or prediction purposes. The phenotype of interest can
                    be multiple categories (classes) or patient survival data."
        ),
        bigdash::sidebarTabHelp(
          ns("cmap-tab"),
          "Similar Experiments",
          tspan("Find similar experiments by correlating their signatures.
                 The main goal is to identify experiments showing similar signatures and find
                 genes that are commonly up/down regulated between experiments.")
        ),
        bigdash::sidebarTabHelp(
          ns("comp-tab"),
          "Compare Datasets",
          "Compare expression and signatures between two datasets,
                    from similar experiments or from different datatypes, e.g. transcriptomics and proteomics."
        ),
        bigdash::sidebarTabHelp(
          ns("tcga-tab"),
          "TCGA Analysis",
          "Correlate your signature with the survival in cancer patients from the TCGA database. Warning: EXPERIMENTAL."
        ),
        bigdash::sidebarTabHelp(
          ns("cell-tab"),
          "Single-Cell Profiling",
          tspan("Visualize the distribution of (inferred)
                    immune cell types, expressed genes and pathway activation.")
        ),
        bigdash::sidebarTabHelp(
          ns("consensus-tab"),
          "Consensus WGCNA",
          tspan("Consensus analysis using the WGCNA framework")
        ),
        bigdash::sidebarTabHelp(
          ns("preservation-tab"),
          "Preservation WGCNA",
          tspan("Preservation analysis using the WGCNA framework")
        ),
        bigdash::sidebarTabHelp(
          ns("ideograms-tab"),
          "Beta Ideograms",
          tspan("Epigenomics visualizations and analyses for methylomics data.")
        ),
        !!!MODULE.multiomics$module_help(ns) ### HELP!!! DOES NOT WORK!!!
      ),
      bigdash::bigTabs(
        id = id,
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
        lock_advanced(filter_boards(bigdash::bigTabItem(ns("dataview-tab"), DataViewInputs(ns("dataview")), create_loader(ns("dataview-loader"))))),
        lock_advanced(filter_boards(MODULE.clustering$module_ui(ns))),
        lock_advanced(filter_boards(MODULE.expression$module_ui(ns))),
        lock_advanced(filter_boards(MODULE.enrichment$module_ui(ns))),
        lock_advanced(filter_boards(MODULE.compare$module_ui(ns))),
        lock_advanced(filter_boards(MODULE.systems$module_ui(ns))),
        lock_advanced(filter_boards(MODULE.multiomics$module_ui(ns))),
        lock_advanced(filter_boards(MODULE.wgcna$module_ui(ns))),
        lock_advanced(filter_boards(MODULE.epigenomics$module_ui(ns)))
      )
    ) ## end of bigPage
  }


  info("[opg_ui] >>> creating UI")
  ui <- createUI(menu_tree)
  info("[opg_ui] <<< finished UI!")

  return(ui)
}
