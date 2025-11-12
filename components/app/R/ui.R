##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

app_ui <- function(x) {
  if (identical("/cookie", x$PATH_INFO)) {
    value <- x$HTTP_HEADER_USER_COOKIE
    return(cookies::set_cookie_response(
      cookie_name = "persistentOPG",
      cookie_value = value,
      http_only = TRUE,
      secure_only = TRUE,
      redirect = "/close",
      same_site = "Strict",
      path = "/"
    ))
  } else if (identical("/cookie_nonce", x$PATH_INFO)) {
    value <- x$HTTP_HEADER_USER_COOKIE
    return(cookies::set_cookie_response(
      cookie_name = "persistentOPG_nonce",
      cookie_value = value,
      http_only = TRUE,
      secure_only = TRUE,
      redirect = "/close",
      same_site = "Strict",
      path = "/"
    ))
  } else if (identical("/cookie_remove", x$PATH_INFO)) {
    return(cookies::set_cookie_response(
      cookie_name = "persistentOPG",
      cookie_value = "",
      expiration = -1,
      http_only = TRUE,
      secure_only = TRUE,
      redirect = "/close",
      same_site = "Strict",
      path = "/"
    ))
  } else if (identical("/close", x$PATH_INFO)) {} else if (identical("/", x$PATH_INFO)) {
    #-------------------------------------------------------
    ## Build USERMENU
    #-------------------------------------------------------
    VERSION <- scan(file.path(OPG, "VERSION"), character())[1]

    upgrade.tab <- NULL
    if (opt$AUTHENTICATION == "firebase") {
      upgrade.tab <- bigdash::navbarDropdownItem(
        "Upgrade",
        onClick = "show_plans()"
      )
    }

    gtag2 <- NULL
    if (Sys.getenv("OMICS_GOOGLE_TAG") != "") {
      ## Add Google Tag manager body code
      gtag2 <- htmltools::includeHTML("www/google-tags-noscript.html")
      gtag2 <- sub("GTM-0000000", Sys.getenv("OMICS_GOOGLE_TAG"), gtag2)
    }

    createUI <- function() {
      message("\n======================================================")
      message("======================= UI ===========================")
      message("======================================================\n")

      version <- scan(file.path(OPG, "VERSION"), character())[1]
      id <- "maintabs"
      header <- shiny::tagList(
        shiny::tags$head(htmltools::includeHTML("www/hubspot-embed.html")),
        ##    gtag2, ## Google Tag Manager???
        shiny::tags$head(shiny::tags$script(src = "custom/temp.js")),
        shiny::tags$head(shiny::tags$script(src = "static/copy-info-helper.js")),
        shiny::tags$script(src = "custom/close-message.js"),
        shiny::tags$head(shiny::tags$script(src = "static/add-tick-helper.js")),
        shiny::tags$head(shiny::tags$script(src = "custom/dropdown-helper.js")),
        shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "custom/styles.min.css")),
        shiny::tags$head(shiny::tags$link(rel = "shortcut icon", href = "custom/favicon.ico")),
        visnetwork = visNetwork::visNetworkOutput("a", height = "0px"),
        shinyjs::useShinyjs(),
        waiter::use_waiter(),
        sever::useSever(),
        bigLoaders::addBigLoaderDeps(),
        firebase::useFirebase(firestore = TRUE, analytics = TRUE),
        shinybrowser::detect(),
        shinybusy::busy_start_up(
          text = tags$h2("\nPrepping your personal playground..."), mode = "auto",
          background = "#2780e3", color = "#ffffff",
          loader = shinybusy::spin_epic("hollow-dots", color = "#FFF")
        )
      )

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

      ## ============================================================================
      ## MENU TREE CONFIGURATION
      ## ============================================================================
      ##
      ## Each menu has two parts:
      ##   • items: The tabs/pages in that menu section
      ##   • icon:  FontAwesome icon shown next to the menu title
      ##
      ## BASIC FORMAT (for simple menus):
      ##   "Menu Name" = list(
      ##     items = c(tab_id = "Display Name"),
      ##     icon = "icon-name"
      ##   )
      ##
      ## SINGLE-ITEM MENU WITH BADGE (badge shows on main menu item):
      ##   "Menu Name" = list(
      ##     items = list(
      ##       list(id = "tab_id", title = "Display Name", badge = "new", badge_color = "success")
      ##     ),
      ##     icon = "icon-name"
      ##   )
      ##
      ## MULTI-ITEM MENU WITH BADGES (badges show on sub-items only):
      ##   "Menu Name" = list(
      ##     items = list(
      ##       list(id = "tab1", title = "Tab 1"),
      ##       list(id = "tab2", title = "Tab 2", badge = "new", badge_color = "success")
      ##     ),
      ##     icon = "icon-name"
      ##   )
      ##
      ## Badge colors: "success" (green), "warning" (orange), "danger" (red), "info" (blue)
      ## Icons: Search FontAwesome 6 for icon names (use without "fa-" prefix)
      ## ============================================================================
      menu_tree <- list(
        "Welcome" = list(
          items = c(welcome = "Welcome"),
          icon = "house"
        ),
        "Datasets" = list(
          items = c(load = "Home"),
          icon = "database"
        ),
        "DataView" = list(
          items = c(dataview = "DataView"),
          icon = "table"
        ),
        "Clustering" = list(
          items = c(
            clustersamples = "Samples",
            clusterfeatures = "Features"
          ),
          icon = "circle-nodes"
        ),
        "Expression" = list(
          items = list(
            list(id = "diffexpr", title = "Differential expression"),
            list(id = "timeseries", title = "TimeSeries", badge = "new", badge_color = "success"),
            list(id = "corr", title = "Correlation analysis"),
            list(id = "bio", title = "Find biomarkers")
          ),
          icon = "dna"
        ),
        "GeneSets" = list(
          items = c(
            enrich = "Geneset Enrichment",
            sig = "Test geneset",
            pathway = "Pathway analysis",
            wordcloud = "Word cloud"
          ),
          icon = "layer-group"
        ),
        "Compare" = list(
          items = c(
            isect = "Compare signatures",
            comp = "Compare datasets",
            cmap = "Similar experiments"
          ),
          icon = "code-compare"
        ),
        "SystemsBio" = list(
          items = c(
            drug = "Drug connectivity",
            cell = "Cell profiling",
            pcsf = "PCSF",
            wgcna = "WGCNA",
            tcga = "TCGA survival (beta)"
          ),
          icon = "diagram-project"
        ),
        "MultiOmics (beta)" = MODULE.multiomics$module_menu()
      )

      ## filter disabled modules
      ENABLED["welcome"] <<- TRUE
      ENABLED["load"] <<- TRUE

      dbg("names(menu_tree) = ", names(menu_tree))
      dbg("names.ENABLED = ", names(ENABLED))
      menu_tree <- menu_tree[MODULES_ENABLED]
      ## menu_tree <- lapply(menu_tree, function(m) m[which(ENABLED[names(m)])])
      ENABLED <<- array(BOARDS %in% sapply(menu_tree, function(m) {
        if (is.list(m) && "items" %in% names(m)) {
          items <- m$items
          ## Handle per-item icon structure (list of lists)
          if (is.list(items) && length(items) > 0 && is.list(items[[1]])) {
            sapply(items, function(x) x$id)
          } else {
            names(items)
          }
        } else {
          ## Fallback for old structure
          names(m)
        }
      }), dimnames = list(BOARDS))

      populateSidebar <- function(menu_tree) {
        sidebar_item <- function(title, name, icon, badge = NULL, badge_color = "success") {
          div(class = "sidebar-item", bigdash::sidebarItem(title, paste0(name, "-tab"), icon = icon, badge = badge, badge_color = badge_color))
        }
        sidebar_menu_item <- function(title, name, icon, badge = NULL, badge_color = "success") {
          bigdash::sidebarMenuItem(title, paste0(name, "-tab"), icon, badge = badge, badge_color = badge_color)
        }
        sidebar_menu_with_items <- function(items, title, menu_icon) {
          ee <- list()

          ## Check if items is a list of lists (with badges) or named vector
          if (is.list(items) && length(items) > 0 && is.list(items[[1]])) {
            ## List structure: list(list(id, title, badge, badge_color), ...)
            for (i in seq_along(items)) {
              item <- items[[i]]
              tab.name <- item$id
              tab.title <- item$title
              tab.badge <- item$badge
              tab.badge_color <- if (!is.null(item$badge_color)) item$badge_color else "success"
              ee[[i]] <- sidebar_menu_item(tab.title, tab.name, NULL, tab.badge, tab.badge_color)
            }
          } else {
            ## Simple named vector: no icons or badges for sub-items
            for (i in seq_along(items)) {
              tab.name <- names(items)[i]
              tab.title <- items[i]
              ee[[i]] <- sidebar_menu_item(tab.title, tab.name, NULL, NULL, "success")
            }
          }

          ## Add icon to the main menu title if provided
          menu_title <- if (!is.null(menu_icon)) {
            shiny::tagList(shiny::icon(menu_icon), " ", title)
          } else {
            title
          }

          bigdash::sidebarMenu(menu_title, !!!ee)
        }

        ## This creates the menu from a menu_tree
        menu <- list()
        for (i in seq_along(menu_tree)) {
          menu.id <- names(menu_tree)[i]
          menu.item <- menu_tree[[i]]

          ## Handle new structure with icon field
          if (is.list(menu.item) && "items" %in% names(menu.item)) {
            items <- menu.item$items
            tab.icon <- menu.item$icon

            ## Extract tab names for length check
            if (is.list(items) && length(items) > 0 && is.list(items[[1]])) {
              tab.names <- sapply(items, function(x) x$id)
            } else {
              tab.names <- names(items)
            }
          } else {
            ## Fallback for old structure (e.g., MODULE.multiomics)
            items <- menu.item
            tab.names <- names(menu.item)
            tab.icon <- NULL
          }

          if (length(tab.names) == 0) {
            ## Empty menu
          } else if (length(tab.names) == 1) {
            ## Single item - extract title, icon, badge properly
            if (is.list(items) && !is.null(items[[1]]$title)) {
              tab.title <- items[[1]]$title
              tab.icon <- if (!is.null(items[[1]]$icon)) items[[1]]$icon else tab.icon
              tab.badge <- items[[1]]$badge
              tab.badge_color <- if (!is.null(items[[1]]$badge_color)) items[[1]]$badge_color else "success"
            } else {
              tab.title <- items[1]
              tab.badge <- NULL
              tab.badge_color <- "success"
            }
            menu[[menu.id]] <- sidebar_item(tab.title, tab.names, tab.icon, tab.badge, tab.badge_color)
          } else {
            ## Multiple items
            menu[[menu.id]] <- sidebar_menu_with_items(items, menu.id, tab.icon)
          }
        }
        return(menu)
      }

      info("[ui.R] creating sidebar menu")
      mm <- populateSidebar(menu_tree)
      mm <- lapply(mm, as.character)
      mm <- HTML(unlist(mm))
      sidebar <- bigdash::sidebar("Menu", mm)

      big_theme2 <- bigdash::big_theme()
      big_theme2 <- bslib::bs_add_variables(big_theme2,
        "grid-breakpoints" = "map-merge($grid-breakpoints, ('xxxl': 2400px))",
        .where = "declarations"
      )

      ## offcanvas chatbox
      div.chirpbutton <- NULL
      if (opt$ENABLE_CHIRP) {
        div.chirpbutton <- shiny::actionButton("chirp_button", "Discuss!",
          width = "auto", class = "quick-button",
          onclick = "window.open('https://www.reddit.com/r/omicsplayground', '_blank')"
        )
      }

      div.invitebutton <- InviteFriendUI("invite")
      div.upgradebutton <- if (opt$ENABLE_UPGRADE) {
        UpgradeModuleUI("upgrade")
      } else {
        NULL
      }

      ## ------------------------- bigPage ----------------------------------
      bigdash.sidebarHelp2 <- function(...) {
        do.call(bigdash::sidebarHelp, rlang::list2(...))
      }

      bigdash::bigPage(
        shiny.i18n::usei18n(i18n),
        header,
        title = "Omics Playground 4",
        theme = big_theme2,
        sidebar = sidebar,
        navbar = bigdash::navbar(
          title = tags$img(
            id = "logo-bigomics",
            ## src = "assets/img/bigomics.png",
            src = "static/bigomics-logo.png",
            height = "30"
            # width = "110",
          ),
          center = tags$div(
            shiny::div(shiny::uiOutput("current_dataset"), class = "current-dataset")
          ),
          left = tags$div(
            style = "padding: 0 0 0 20px;",
            div(
              style = "display: inline-block; ",
              bigdash::navbarDropdown(
                "Datasets",
                style = "border: 1px; padding: 2px 6px;",
                bigdash::navbarDropdownTab(
                  "Upload new",
                  "upload-tab"
                ),
                bigdash::navbarDropdownTab(
                  "Load from library",
                  "load-tab"
                ),
                bigdash::navbarDropdownTab(
                  "Shared datasets",
                  "sharing-tab"
                )
              )
            )
          ),
          div.upgradebutton,
          div.invitebutton,
          div.chirpbutton,
          div(
            id = "mainmenu_help",
            bigdash::navbarDropdown(
              "Help",
              bigdash::navbarDropdownItem(
                "Documentation",
                link = "https://omicsplayground.readthedocs.io",
                target = "_blank"
              ),
              bigdash::navbarDropdownItem(
                "Video tutorials",
                link = "https://bigomics.ch/tutorials/",
                target = "_blank"
              ),
              bigdash::navbarDropdownItem(
                "Google forum",
                link = "https://groups.google.com/d/forum/omicsplayground",
                target = "_blank"
              ),
              bigdash::navbarDropdownItem(
                "Submit a support ticket",
                link = "https://share-eu1.hsforms.com/1glP7Cm6GQrWIGXgZrC0qrweva7t",
                target = "_blank"
              ),
              bigdash::navbarDropdownItem(
                "Github issues",
                link = "https://github.com/bigomics/omicsplayground/issues",
                target = "_blank"
              ),
              bigdash::navbarDropdownItem(
                "Case studies",
                link = "https://bigomics.ch/case-studies/",
                target = "_blank"
              )
            )
          ),
          div(
            id = "mainmenu_user",
            bigdash::navbarDropdown(
              ## "User",
              shiny::textOutput("current_user", inline = TRUE),
              bigdash::navbarDropdownTab(
                "User profile",
                "userprofile-tab"
              ),
              bigdash::navbarDropdownTab(
                "App settings",
                "usersettings-tab"
              ),
              upgrade.tab,
              tags$li(
                actionLink("navbar_about", "About")
              ),
              logout.tab
            )
          ),
          div(
            id = "mainmenu_appsettings",
            bigdash::navbarDropdown(
              auto_close = "outside",
              shiny::icon("cog"),
              div(
                class = "dropdown-items",
                bslib::input_switch("enable_beta", "Enable beta features"),
                bslib::input_switch("enable_info", "Show info boxes", value = TRUE),
                selector_switch(
                  class = "card-footer-checked",
                  label = "show captions",
                  is.checked = FALSE
                )
              ),
              bigdash::navbarDropdownItem(
                withTooltip(
                  shiny::selectInput(
                    inputId = "selected_labeltype",
                    label = "Label type:",
                    choices = c("feature", "symbol", "name"),
                    selected = "feature",
                    width = "100%"
                  ),
                  "Choose a label type to be displayed in the plots",
                  placement = "right", options = list(container = "body")
                )
              )
            )
          ),
          ## THIS IS SO WEIRD. if we remove/comment out the
          ## prettySwitch, the header of all plotModules f*ck
          ## up... (IK). HELP!!! we do not need this button...
          div(
            style = "visibility: hidden; display: none;",
            shinyWidgets::prettySwitch("I_AM_WEIRD_BUTTON", "remove me")
          )
        ),
        settings = bigdash::settings(
          "Settings"
        ),
        ## bigdash::sidebarHelp(
        bigdash.sidebarHelp2(
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
            "Analyze dataset",
            "This panel shows the available datasets within the platform. These data sets
                    have been pre-computed and are ready to be used. Select a
                    dataset in the table and load the data set by clicking the 'Analyze dataset' button."
          ),
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
                    The main goal is to identify experiments showing similar signatures and find genes
                    that are commonly up/down regulated between experiments.")
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
          !!!MODULE.multiomics$module_help() ### HELP!!! DOES NOT WORK!!!
        ),
        bigdash::bigTabs(
          bigdash::bigTabItem(
            "welcome-tab",
            WelcomeBoardInputs("welcome"),
            WelcomeBoardUI("welcome")
          ),
          bigdash::bigTabItem(
            "load-tab",
            # LoadingInputs("load")
            LoadingUI("load")
          ),
          bigdash::bigTabItem(
            "upload-tab",
            UploadUI("upload")
          ),
          bigdash::bigTabItem(
            "userprofile-tab",
            UserProfileUI("user_profile")
          ),
          bigdash::bigTabItem(
            "usersettings-tab",
            AppSettingsInputs("app_settings"),
            AppSettingsUI("app_settings")
          ),
          bigdash::bigTabItem(
            "sharing-tab",
            SharedDatasetsUI("load")
          )
        )
        ## UploadUI("upload")
      ) ## end of bigPage
    }

    info("[ui.R] >>> creating UI")
    ui <- createUI()
    info("[ui.R] <<< finished UI!")

    return(ui)
  }
}

