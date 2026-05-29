
## NOTE: This is not a real shiny module (yet...). We should move as
## much as possible OPG server related code here.

opg_server <- function(input, output, session, PGX, env, auth) {

  labeltype <- reactiveVal("feature") # can be feature (rownames counts), symbol or name

  ## Hide/show tabs. Open sidebar and settings
  shiny::observeEvent(
    {
      list(
        auth$logged,
        env$user_settings$enable_beta(),
        env$trigger_on_change_dataset()
      )
    },
    {
      ## trigger on change dataset
      shiny::req(PGX$X)
      info("[SERVER] trigger on change dataset")

      tab_control()

      ## hide all main tabs until we have an object
      if (is.null(PGX) || is.null(PGX$name) || !auth$logged) {
        warning("[SERVER] !!! no data. hiding menu.")
        shinyjs::runjs("sidebarClose()")
        shinyjs::runjs("settingsClose()")
        bigdash.selectTab(session, selected = "dataview-tab")        
        return(NULL)
      }

      ## show all main tabs
      shinyjs::runjs("sidebarOpen()")
      shinyjs::runjs("settingsOpen()")

    }
  )

  ## Hide/show basic or full menu
  observeEvent(input$menu_basic, {
    if (isTRUE(input$menu_basic)) {
      shinyjs::runjs("$('#menu-basic').removeClass('nodisp').show(); $('#menu-full').addClass('nodisp').hide();")
    } else {
      shinyjs::runjs("$('#menu-full').removeClass('nodisp').show(); $('#menu-basic').addClass('nodisp').hide();")
    }
  }, ignoreInit = FALSE)

  ## Modules needed after dataset is loaded (deferred) --------------
  observeEvent(env$load$is_data_loaded(), {
    # depending on datatpye, subset modules enabled and create modules active,
    if (tolower(PGX$datatype) == "multi-omics") {
      MODULES_ACTIVE <- MODULES_MULTIOMICS
    } else {
      MODULES_ACTIVE <- MODULES_TRANSCRIPTOMICS
    }
    if (env$load$is_data_loaded() == 1) {
      bigdash.hideMenuElement(session, "Clustering")
      bigdash.hideMenuElement(session, "Expression")
      bigdash.hideMenuElement(session, "GeneSets")
      bigdash.hideMenuElement(session, "Compare")
      bigdash.hideMenuElement(session, "SystemsBio")
      bigdash.hideMenuElement(session, "MultiOmics")
      bigdash.hideMenuElement(session, "WGCNA")
    }
    # ###################### I STILL HAVE TO REMOVE THE UI!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MODULES_TO_REMOVE <- xor(MODULES_LOADED, MODULES_ACTIVE) & MODULES_LOADED
    MODULES_TO_LOAD <- xor(MODULES_LOADED, MODULES_ACTIVE) & MODULES_ACTIVE

    lapply(names(MODULES_TO_REMOVE[MODULES_TO_REMOVE]), function(x) {
      if (x == "DataView") {
        bigdash.removeTab(session, "dataview-tab")
        bigdash.hideMenuElement(session, "DataView")
      }
      if (x == "Clustering") {
        lapply(names(MODULE.clustering$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "Clustering")
        loaded$clustering <- 0
      }
      if (x == "Expression") {
        lapply(names(MODULE.expression$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "Expression")
        loaded$expression <- 0
      }
      if (x == "GeneSets") {
        lapply(names(MODULE.enrichment$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "GeneSets")
        loaded$enrichment <- 0
      }
      if (x == "Compare") {
        lapply(names(MODULE.compare$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "Compare")
        loaded$compare <- 0
      }
      if (x == "SystemsBio") {
        lapply(names(MODULE.systems$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "SystemsBio")
        loaded$systems <- 0
      }
      if (x == "MultiOmics") {
        lapply(names(MODULE.multiomics$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "MultiOmics")
        loaded$multiomics <- 0
      }
      if (x == "WGCNA") {
        lapply(names(MODULE.wgcna$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "WGCNA")
        loaded$wgcna <- 0
      }
    })

    if (env$load$is_data_loaded()) { # == 1) {

      additional_ui_tabs <- list(
        dataview = bigdash::bigTabItem(
          "dataview-tab",
          DataViewInputs("dataview"),
          DataViewUI("dataview")
        )
      )

      insertBigTabUI <- function(ui) {
        for (i in 1:length(ui)) {
          shiny::insertUI(
            selector = "#big-tabs",
            where = "beforeEnd",
            ui = ui[[i]],
            immediate = TRUE
          )
        }
        bigdash.openSettings()
      }
      insertBigTabItem <- function(tab) {
        insertBigTabUI(additional_ui_tabs[tab])
      }

      ## show the dashboard item in app sidebar
      #bslib::nav_show("app-sidebar", "Dashboard")
      shinyjs::enable(selector = "a[data-value='Dashboard']")
      shinyjs::enable(selector = "a[data-value='Studio']")
      shinyjs::enable(selector = "a[data-value='Copilot']")
      shinyjs::enable(selector = "a[data-value='Copilot2']")      
      
      shiny::withProgress(
        message = "Preparing your dashboard (server)...",
        value = 0,
        {
          if (MODULES_TO_LOAD["DataView"]) {
            info("[SERVER] calling DataView module")
            insertBigTabItem("dataview")
            DataViewBoard("dataview",
              pgx = PGX, labeltype = labeltype
            )
            bslib::nav_select("app-sidebar", selected = "Dashboard")
            bigdash.selectTab(session, "dataview-tab")            
          }
          shiny::incProgress(0.1)

          if (MODULES_TO_LOAD["Clustering"]) {
            mod <- MODULE.clustering
            insertBigTabUI(mod$module_ui())
            info("[SERVER:UI:1] calling Clustering module")
            bigdash.showMenuElement(session, "Clustering")
            lapply(names(MODULE.clustering$module_menu()), function(x) {
              bigdash.showTab(session, paste0(x, "-tab"))
            })
          }
          shiny::incProgress(0.1)

          if (MODULES_TO_LOAD["Expression"]) {
            mod <- MODULE.expression
            insertBigTabUI(mod$module_ui())
            info("[SERVER:UI:1] calling Expression module")
            info("[SERVER:UI:1] calling DiffExprBoard module")
            bigdash.showMenuElement(session, "Expression")
            ExpressionBoard("diffexpr",
              pgx = PGX, labeltype = labeltype
            ) ->> env$diffexpr
            lapply(names(MODULE.expression$module_menu()), function(x) {
              bigdash.showTab(session, paste0(x, "-tab"))
            })
          }
          shiny::incProgress(0.1)

          if (MODULES_TO_LOAD["GeneSets"]) {
            mod <- MODULE.enrichment
            insertBigTabUI(mod$module_ui())
            info("[SERVER:UI:1] calling GeneSets module")
            info("[SERVER:UI:1] calling EnrichmentBoard module")
            bigdash.showMenuElement(session, "GeneSets")
            EnrichmentBoard("enrich",
              pgx = PGX,
              selected_gxmethods = env$diffexpr$selected_gxmethods
            ) ->> env$enrich
            lapply(names(MODULE.enrichment$module_menu()), function(x) {
              bigdash.showTab(session, paste0(x, "-tab"))
            })
          }
          shiny::incProgress(0.1)

          if (MODULES_TO_LOAD["Compare"]) {
            mod <- MODULE.compare
            insertBigTabUI(mod$module_ui())
            info("[SERVER:UI:1] calling Compare module")
            bigdash.showMenuElement(session, "Compare")
            lapply(names(MODULE.compare$module_menu()), function(x) {
              bigdash.showTab(session, paste0(x, "-tab"))
            })
          }
          shiny::incProgress(0.1)

          if (MODULES_TO_LOAD["SystemsBio"]) {
            mod <- MODULE.systems
            insertBigTabUI(mod$module_ui())
            info("[SERVER:UI:1] calling SystemsBio module")
            bigdash.showMenuElement(session, "SystemsBio")
            lapply(names(mod$module_menu()), function(x) {
              bigdash.showTab(session, paste0(x, "-tab"))
            })
            bigdash.toggleTab(session, "tcga-tab", env$user_settings$enable_beta() && dir.exists(file.path(OPG, "libx")))
          }
          shiny::incProgress(0.1)

          if (MODULES_TO_LOAD["MultiOmics"] && exists("MODULE.multiomics")) {
            info("[SERVER:UI:1] initializing MultiOmics module")
            mod <- MODULE.multiomics
            insertBigTabUI(mod$module_ui())
            bigdash.showMenuElement(session, "MultiOmics")
            lapply(names(MODULE.multiomics$module_menu()), function(x) {
              bigdash.showTab(session, paste0(x, "-tab"))
            })
          }

          if (MODULES_TO_LOAD["WGCNA"] && exists("MODULE.wgcna")) {
            info("[SERVER:UI:1] initializing WGCNA module")
            mod <- MODULE.wgcna
            insertBigTabUI(mod$module_ui())
            bigdash.showMenuElement(session, "WGCNA")
            lapply(names(MODULE.wgcna$module_menu()), function(x) {
              bigdash.showTab(session, paste0(x, "-tab"))
            })
          }

          MODULES_LOADED <<- MODULES_ACTIVE

          if (env$load$is_data_loaded() > 0) {
            env$trigger_on_change_dataset(runif(1))
          }
          info("[SERVER:UI:1] calling modules done!")
        }
      )
    }

    if (env$load$is_data_loaded() == 1) {
      # this is a function - like "handleSettings()" in bigdash- needed to
      # make the settings sidebar show up for the inserted tabs
      shinyjs::runjs(
        "  $('.big-tab')
          .each((index, el) => {
            let settings = $(el)
              .find('.tab-settings')
              .first();
            $(settings).data('target', $(el).data('name'));
            $(settings).appendTo('#settings-content');
          });"
      )
    }

    ## Goto dataview
    bslib::nav_select("app-sidebar", selected = "Dashboard")
    bigdash.openSettings(lock = TRUE)
    bigdash.openSidebar()
    bigdash.showTabs(session)
    bigdash.selectTab(session, "dataview-tab")

    ## remove loading modal from LoadingBoard
    shinyjs::delay(2000, {
      shiny::removeModal()
    })
    
  })

  insertBigTabUI2 <- function(ui, menu) {
    for (i in 1:length(ui)) {
      for (j in 2:length(ui[[i]])) {
        shiny::insertUI(
          selector = paste0("div.big-tab[data-name='", ui[[i]][[1]], "']"),
          where = "beforeEnd",
          ui = ui[[i]][[j]],
          immediate = TRUE
        )
      }
    }
    shinyjs::runjs(
      "  $('.big-tab')
    .each((index, el) => {
      let settings = $(el)
        .find('.tab-settings')
        .first();
      $(settings).data('target', $(el).data('name'));
      $(settings).appendTo('#settings-content');
    });"
    )
    bigdash.openSettings()
    shinyjs::hide(selector = paste0("[id='", names(menu), "-loader']"))
  }

  loaded <- shiny::reactiveValues(
    clustering = 0,
    expression = 0,
    enrichment = 0,
    compare = 0,
    systems = 0,
    multiomics = 0,
    wgcna = 0
  )

  observeEvent(input$nav, {
    dbg("[SERVER] input$nav =", input$nav)

    if (input$nav %in% c("clustersamples-tab", "clusterfeatures-tab") &&
      loaded$clustering == 0) {
      info("[SERVER:UI:2] reacted: calling Clustering module")
      mod <- MODULE.clustering
      insertBigTabUI2(mod$module_ui2(), mod$module_menu())
      mod$module_server(PGX, labeltype = labeltype)
      loaded$clustering <- 1
      tab_control()
    }
    if (input$nav %in% c("diffexpr-tab", "corr-tab", "bio-tab", "timeseries-tab") &&
      loaded$expression == 0) {
      info("[SERVER:UI:2] reacted: calling Expression module")
      mod <- MODULE.expression
      insertBigTabUI2(mod$module_ui2(), mod$module_menu())
      mod$module_server(PGX, labeltype = labeltype)
      loaded$expression <- 1
      tab_control()
    }
    if (input$nav %in% c("enrich-tab", "sig-tab", "pathway-tab", "wordcloud-tab") &&
      loaded$enrichment == 0) {
      info("[SERVER:UI:2] reacted: calling Enrichment module")
      mod <- MODULE.enrichment
      insertBigTabUI2(mod$module_ui2(), mod$module_menu())
      mod$module_server(PGX, labeltype = labeltype, env = env)
      loaded$enrichment <- 1
      tab_control()
    }
    if (input$nav %in% c("isect-tab", "comp-tab", "cmap-tab") && loaded$compare == 0) {
      info("[SERVER:UI] reacted: calling Compare module")
      mod <- MODULE.compare
      insertBigTabUI2(mod$module_ui2(), mod$module_menu())
      mod$module_server(PGX, labeltype = labeltype, auth = auth, env = env, reload_pgxdir = reload_pgxdir)
      loaded$compare <- 1
      tab_control()
    }
    if (input$nav %in% c("drug-tab", "tcga-tab", "cell-tab", "pcsf-tab") &&
      loaded$systems == 0) {
      info("[SERVER:UI:2] reacted: calling Systems module")
      mod <- MODULE.systems
      insertBigTabUI2(mod$module_ui2(), mod$module_menu())
      mod$module_server(PGX)
      loaded$systems <- 1
      tab_control()
    }
    if (input$nav %in% c(
      "mofa-tab", "mgsea-tab", "snf-tab", "lasagna-tab",
      "deepnet-tab"
    ) && loaded$multiomics == 0) {
      info("[SERVER:UI:2] reacted: calling Multi-Omics module")
      mod <- MODULE.multiomics
      insertBigTabUI2(mod$module_ui2(), mod$module_menu())
      mod$module_server(PGX)
      loaded$multiomics <- 1
      tab_control()
    }
    if (input$nav %in% c(
      "wgcna-tab", "mwgcna-tab", "consensus-tab",
      "preservation-tab"
    ) && loaded$wgcna == 0) {
      info("[SERVER:UI:2] reacted: calling WGCNA module")
      mod <- MODULE.wgcna
      insertBigTabUI2(mod$module_ui2(), mod$module_menu())
      mod$module_server(PGX)
      loaded$wgcna <- 1
      tab_control()
    }
  })


  tab_control <- function() {

    ## show beta feauture
    show.beta <- env$user_settings$enable_beta()
    if (is.null(show.beta) || length(show.beta) == 0) show.beta <- FALSE

    has.libx <- dir.exists(file.path(OPG, "libx"))
    is.multiomics <- tolower(PGX$datatype) == "multi-omics"
    
    ## Hide beta main tabs
    bigdash.toggleTab(session, "tcga-tab", show.beta && has.libx)
    bigdash.toggleTab(session, "consensus-tab", show.beta)
    bigdash.toggleTab(session, "preservation-tab", opt$DEVMODE && show.beta)
    #bigdash.toggleTab(session, "mwgcna-tab", show.beta && is.multiomics)
    bigdash.toggleTab(session, "mwgcna-tab", is.multiomics)

    ## hide beta subtabs..
    toggleTab("drug-tabs", "Connectivity map (beta)", show.beta) ## too slow
    toggleTab("pathway-tabs", "Enrichment Map (beta)", show.beta) ## too slow
    toggleTab("wgcna-tabs", "AI Report✨", show.beta)
    toggleTab("drug-tabs", "AI Summary✨", show.beta)     

    ## Control tab to only be displayed if there is custom fc + baseline fc
    has.customfc <- "custom" %in% colnames(PGX$gx.meta$meta[[1]]$fc) &&
      length(colnames(PGX$gx.meta$meta[[1]]$fc)) > 1
    toggleTab("diffexpr-tabs1", "FC-FC comparison", has.customfc)

    ## Dynamically show upon availability in pgx object
    tabRequire(PGX, session, "drug-tab", "drugs", TRUE)
    tabRequire(PGX, session, "wordcloud-tab", "wordcloud", TRUE)
    tabRequire(PGX, session, "cell-tab", "deconv", TRUE)
    tabRequireTS(PGX, session, "timeseries-tab", TRUE)
    tabRequire(PGX, session, "cmap-tab", "connectivity", TRUE)
    gset_tabs <- c("enrich-tab", "pathway-tab", "isect-tab", "sig-tab")
    for (tab_i in gset_tabs) {
      tabRequire(PGX, session, tab_i, "gsetX", TRUE)
      tabRequire(PGX, session, tab_i, "gset.meta", TRUE)
    }

    ## Hide PCSF and WGCNA for metabolomics.
    # WGCNA will be available upon gmt refactoring
    if (PGX$datatype == "metabolomics") {
      info("[SERVER] disabling modules for metabolomics data")
      bigdash.hideTab(session, "cmap-tab")
    }

    if (PGX$datatype == "multi-omics") {
      info("[SERVER] disabling modules for multi-omics data")
      bigdash.hideTab(session, "drug-tab")
      bigdash.hideTab(session, "cell-tab")
      bigdash.hideTab(session, "wordcloud-tab")
      bigdash.hideTab(session, "cmap-tab")
    }
    
  }

  ## -------------------------------------------------------------
  ## Labeltype stuff
  ## -------------------------------------------------------------

  # populate labeltype selector based on pgx$genes
  observeEvent(
    {
      list(PGX$X, PGX$name)
    },
    {
      req(PGX$genes)
      genes_mat <- PGX$genes
      # remove NA columns and columns with only 1 unique value
      genes_mat <- genes_mat[, colMeans(is.na(genes_mat)) < 1, drop = FALSE]
      genes_mat <- genes_mat[, sapply(genes_mat, function(x) length(unique(x)) > 1), drop = FALSE]
      genes_mat <- genes_mat[, !duplicated(t(genes_mat)), drop = FALSE]
      label_types <- colnames(genes_mat)
      names(label_types) <- label_types
      names(label_types)[names(label_types) == "gene_title"] <- "title"
      label_types <- label_types[!grepl("pos|map|tx_len|source", names(label_types))]
      names(label_types) <- sub("^chr$", "chromosome", names(label_types))
      # if one of the label_types unique values amounts for less than 10% of total genes, remove it
      n_genes <- nrow(PGX$genes)
      keep_types <- sapply(label_types, function(col) {
        if (col %in% colnames(PGX$genes)) {
          n_unique <- length(unique(PGX$genes[, col]))
          (n_unique / n_genes) >= 0.10
        } else {
          TRUE
        }
      })
      label_types <- label_types[keep_types]

      sel.labeltype <- "feature"
      shiny::updateSelectInput(
        session,
        "selected_labeltype",
        choices = label_types,
        selected = sel.labeltype
      )
    }
  )

  # change label type based on selected input
  shiny::observeEvent(
    {
      input$selected_labeltype
    },
    {
      labeltype(input$selected_labeltype)
      if (!is.null(PGX$genes)) {
        lab <- input$selected_labeltype
        if (lab == "gene_title") {
          tt <- paste0(PGX$genes[, "gene_title"], " (", PGX$genes[, "symbol"], ")")
          PGX$genes$gene_name <- tt
        } else if (lab %in% colnames(PGX$genes)) {
          PGX$genes$gene_name <- PGX$genes[, lab]
        } else {
          PGX$genes$gene_name <- rownames(PGX$genes)
        }
      }
    }
  )

  
  
}
