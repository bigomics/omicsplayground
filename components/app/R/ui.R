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
      same_site = "Strict"
    ))
  } else if (identical("/cookie_nonce", x$PATH_INFO)) {
    value <- x$HTTP_HEADER_USER_COOKIE
    return(cookies::set_cookie_response(
      cookie_name = "persistentOPG_nonce",
      cookie_value = value,
      http_only = TRUE,
      secure_only = TRUE,
      redirect = "/close",
      same_site = "Strict"
    ))
  } else if (identical("/cookie_remove", x$PATH_INFO)) {
    return(cookies::set_cookie_response(
      cookie_name = "persistentOPG",
      cookie_value = "",
      expiration = -1,
      http_only = TRUE,
      secure_only = TRUE,
      redirect = "/close"
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
        onClick = "logoutInApp()"
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

      menu_tree <- list(
        "Welcome" = c(
          welcome = "Welcome"
        ),
        "Datasets" = c(
          load = "Home"
        ),
        "DataView" = c(
          dataview = "DataView"
        ),
        "Clustering" = c(
          clustersamples = "Samples",
          clusterfeatures = "Features"
        ),
        "Expression" = c(
          diffexpr = "Differential expression",
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
          wgcna = "WGCNA",
          tcga = "TCGA survival (beta)"
        )
      )

      ## filter disabled modules
      ENABLED["welcome"] <<- TRUE
      ENABLED["load"] <<- TRUE

      menu_tree <- lapply(menu_tree, function(m) m[which(ENABLED[names(m)])])

      populateSidebar <- function(menu_tree) {
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

        ## This creates the menu from a menu_tree
        menu <- list()
        i <- 3
        for (i in 1:length(menu_tree)) {
          tab.names <- names(menu_tree[[i]])
          tab.titles <- menu_tree[[i]]
          menu.id <- names(menu_tree)[i]
          if (length(tab.names) == 0) {
          } else if (length(tab.names) == 1) {
            menu[[menu.id]] <- sidebar_item(tab.titles, tab.names)
          } else {
            menu[[menu.id]] <- sidebar_menu_with_items(menu_tree[[i]], menu.id)
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

      ## ## offcanvas chatbox
      ## div.chirpbutton <- NULL
      ## if (opt$ENABLE_CHIRP) {
      ##   div.chirpbutton <- shiny::actionButton("chirp_button", "Discuss!",
      ##     width = "auto", class = "quick-button",
      ##     onclick = "window.open('https://www.reddit.com/r/omicsplayground', '_blank')"
      ##   )
      ## }

      div.invitebutton <- InviteFriendUI("invite")

      div.copilotbutton <- NULL
      if(opt$DEVMODE) {
        div.copilotbutton <- uiOutput("copilot_button")
      }
      
      ## ------------------------- bigPage ----------------------------------
      bigdash::bigPage(
        shiny.i18n::usei18n(i18n),
        header,
        title = "Omics Playground v3",
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
            ## div(style = "display: inline-block; ",
            ## bigdash::navbarDropdown(
            ##   "Actions",
            ##   style = "display: inline-block; border: 1px; padding: 2px 6px;",
            ##   tags$li(
            ##     actionLink("menu_createreport", "Create report")
            ##   ),
            ##   tags$li(
            ##     actionLink("menu_reanalyze", "Reanalyze")
            ##   )
            ## ))
          ),
          div.copilotbutton,
          div.invitebutton,          
##        div.chirpbutton,
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
                link = "https://www.youtube.com/watch?v=_Q2LJmb2ihU&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-",
                target = "_blank"
              ),
              bigdash::navbarDropdownItem(
                "Google forum",
                link = "https://groups.google.com/d/forum/omicsplayground",
                target = "_blank"
              ),
              bigdash::navbarDropdownItem(
                "Discuss on Reddit",
                link = "https://www.reddit.com/r/omicsplayground",
                target = "_blank"
              ),
              bigdash::navbarDropdownItem(
                "Github issues",
                link = "https://github.com/bigomics/omicsplayground/issues",
                target = "_blank"
              ),
              bigdash::navbarDropdownItem(
                "Case studies",
                link = "https://bigomics.ch/blog/category/case-study/",
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
                    selected = NULL,
                    width = "100%"
                  ),
                  "Choose a label type to be displayed in the heatmap.",
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
            "Perform specialized functional analysis
                    to understand biological pathways including GO, KEGG, and drug connectivity mapping."
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
          )
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
