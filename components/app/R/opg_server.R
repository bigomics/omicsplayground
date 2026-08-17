
## NOTE: This is not a real shiny module (yet...). We should move as
## much as possible OPG server related code here.

opg_server <- function(input, output, session, PGX, env, auth, reload_pgxdir) {

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
    } else if (tolower(PGX$datatype) == "methylomics") {
      MODULES_ACTIVE <- MODULES_METHYLOMICS
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
      bigdash.hideMenuElement(session, "Epigenomics")
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
        reset_group_tabs(MODULE.clustering)
      }
      if (x == "Expression") {
        lapply(names(MODULE.expression$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "Expression")
        loaded$expression <- 0
        reset_group_tabs(MODULE.expression)
      }
      if (x == "GeneSets") {
        lapply(names(MODULE.enrichment$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "GeneSets")
        loaded$enrichment <- 0
        reset_group_tabs(MODULE.enrichment)
      }
      if (x == "Compare") {
        lapply(names(MODULE.compare$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "Compare")
        loaded$compare <- 0
        reset_group_tabs(MODULE.compare)
      }
      if (x == "SystemsBio") {
        lapply(names(MODULE.systems$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "SystemsBio")
        loaded$systems <- 0
        reset_group_tabs(MODULE.systems)
      }
      if (x == "MultiOmics") {
        lapply(names(MODULE.multiomics$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "MultiOmics")
        loaded$multiomics <- 0
        reset_group_tabs(MODULE.multiomics)
      }
      if (x == "WGCNA") {
        lapply(names(MODULE.wgcna$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "WGCNA")
        loaded$wgcna <- 0
        reset_group_tabs(MODULE.wgcna)
      }
      if (x == "Epigenomics") {
        lapply(names(MODULE.epigenomics$module_menu()), function(x) {
          bigdash.removeTab(session, paste0(x, "-tab"))
        })
        bigdash.hideMenuElement(session, "Epigenomics")
        loaded$epigenomics <- 0
        reset_group_tabs(MODULE.epigenomics)
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
            ## preload expression board
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
            ## preload enrichment board
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

          if (MODULES_TO_LOAD["Epigenomics"] && exists("MODULE.epigenomics")) {
            info("[SERVER:UI:1] initializing Epigenomics module")
            mod <- MODULE.epigenomics
            insertBigTabUI(mod$module_ui())
            bigdash.showMenuElement(session, "Epigenomics")
            lapply(names(MODULE.epigenomics$module_menu()), function(x) {
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

  ## `only`: insert just this tab's UI instead of the whole group's. Building a
  ## group's UI is the larger half of opening it -- measured 3.9-6.0s of a
  ## 5.7-9.3s first click, against 1.8-3.3s for starting its servers -- and
  ## most of it is for boards the user has not asked for.
  insertBigTabUI2 <- function(ui, menu, only = NULL) {
    if (!is.null(only)) {
      keep <- vapply(ui, function(x) identical(x[[1]], only), logical(1))
      ui <- ui[keep]
      if (!length(ui)) {
        return(invisible(NULL))
      }
    }
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
    wgcna = 0,
    epigenomics = 0
  )

  ## ---------------------------------------------------------------
  ## Per-board loading
  ## ---------------------------------------------------------------
  ##
  ## A board group used to be materialised whole: clicking one tab built and
  ## started every board in its group. Open Differential Expression and you
  ## also paid for Correlation, Biomarker and TimeSeries.
  ##
  ## load_board() does one board. Both halves have to move together -- a board
  ## server whose UI is not in the DOM loses the updateSelectInput() calls it
  ## makes while initialising -- so the UI is inserted first, then that board's
  ## server is started, then the tab is latched.
  ##
  ## Which boards a group can load this way is given by its `module_servers`,
  ## a list of server functions keyed by tab name. A group without one is left
  ## to the eager path.
  loaded_tabs <- new.env(parent = emptyenv())

  load_board <- function(tab, mod, run_server) {
    if (is.null(tab) || isTRUE(loaded_tabs[[tab]])) {
      return(invisible(FALSE))
    }
    fn <- mod$module_servers[[tab]]
    if (is.null(fn)) {
      return(invisible(FALSE))
    }
    ## Latch before running: a board that errors must not re-insert its UI on
    ## every later click.
    loaded_tabs[[tab]] <- TRUE
    info("[SERVER:UI:2] loading board ", tab)
    insertBigTabUI2(mod$module_ui2(), mod$module_menu(), only = tab)
    run_server(fn)
    tab_control()
    invisible(TRUE)
  }

  ## Called wherever a group is torn down (datatype change), so its boards load
  ## again next time they are opened.
  reset_group_tabs <- function(mod) {
    for (tab in names(mod$module_servers)) {
      if (!is.null(loaded_tabs[[tab]])) rm(list = tab, envir = loaded_tabs)
    }
  }

  observeEvent(input$nav, {
    dbg("[SERVER] input$nav =", input$nav)

    ## Phase 1 (at dataset load) inserts the tab shells; before that there is
    ## nothing for insertUI() to target. Acting earlier would mark the board
    ## loaded while its UI silently went nowhere, and it would never appear --
    ## reachable by clicking a sidebar item while the dashboard is still
    ## initialising.
    shiny::req(PGX$name)

    if (input$nav %in% names(MODULE.clustering$module_servers)) {
      load_board(input$nav, MODULE.clustering, function(f) f(PGX, labeltype = labeltype))
    }
    if (input$nav %in% names(MODULE.expression$module_servers)) {
      load_board(input$nav, MODULE.expression, function(f) f(PGX, labeltype = labeltype))
    }
    if (input$nav %in% names(MODULE.enrichment$module_servers)) {
      load_board(input$nav, MODULE.enrichment, function(f) f(PGX, labeltype = labeltype, env = env))
    }
    if (input$nav %in% names(MODULE.compare$module_servers)) {
      load_board(input$nav, MODULE.compare, function(f) f(PGX, labeltype = labeltype, auth = auth, env = env, reload_pgxdir = reload_pgxdir))
    }
    if (input$nav %in% names(MODULE.systems$module_servers)) {
      load_board(input$nav, MODULE.systems, function(f) f(PGX))
    }
    if (input$nav %in% names(MODULE.multiomics$module_servers)) {
      load_board(input$nav, MODULE.multiomics, function(f) f(PGX))
    }
    if (input$nav %in% names(MODULE.wgcna$module_servers)) {
      load_board(input$nav, MODULE.wgcna, function(f) f(PGX, save_pgx = env$save_pgx))
    }
    if (input$nav %in% names(MODULE.epigenomics$module_servers)) {
      load_board(input$nav, MODULE.epigenomics, function(f) f(PGX))
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
    bigdash.toggleTab(session, "wgcna-tab", !is.multiomics || show.beta)

    ## hide beta subtabs..
    toggleTab("drug-tabs", "Connectivity map (beta)", show.beta) ## too slow
    toggleTab("pathway-tabs", "Enrichment Map (beta)", show.beta) ## too slow
    toggleTab("wgcna-tabs", "AI Report✨", show.beta)
    toggleTab("mwgcna-tabs", "AI Report✨", show.beta)
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

    ## Show Epigenomics only for methylomics data
    if (!is.null(PGX$datatype) && tolower(PGX$datatype) != "methylomics") {
      bigdash.hideTab(session, "ideograms-tab")
      bigdash.hideMenuElement(session, "Epigenomics")
    }

    ## Hide PCSF for methylomics DMP (CpG probe level — no meaningful PPI matching)
    if (!is.null(PGX$datatype) && tolower(PGX$datatype) == "methylomics") {
      is_dmp <- if (!is.null(PGX$dma)) {
        PGX$dma == "Differentially methylated positions"
      } else {
        ## fallback for old pgx files without dma field: CpG probe IDs start with "cg"
        mean(grepl("^cg[0-9]+", rownames(PGX$X))) > 0.5
      }
      if (is_dmp) {
        bigdash.hideTab(session, "pcsf-tab")
      } else {
        bigdash.showTab(session, "pcsf-tab")
      }
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
