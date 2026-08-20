
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

  ## Register all lazy tabs at startup.  bigTabsLazy() will
  ## materialise each tab's heavy UI (via module_lazy closures)
  ## the first time the user clicks on it.
  ## ------------------------------------------------------------
  lazy_tabs <- list()
  
  lazy_tabs[["dataview-tab"]] <- list(
    ui      = function() DataViewUI("dataview"),
    server  = function() DataViewBoard("dataview", pgx = PGX, labeltype = labeltype),
    preload = TRUE
  )

  lazy_tabs[["diffexpr-tab"]] <- list(
    ui      = function() ExpressionUI("diffexpr"),
    server  = function() ExpressionBoard("diffexpr", pgx = PGX, labeltype = labeltype) ->>
                           env$diffexpr,
    preload = TRUE
  )

  lazy_tabs[["enrich-tab"]] = list(
    ui     = function() EnrichmentUI("enrich"),
    server = function() EnrichmentBoard("enrich", pgx = PGX, selected_gxmethods = env$diffexpr$selected_gxmethods) ->>
                          env$enrich,
    preload = TRUE
  )
  
  lazy_tabs <- c(lazy_tabs, MODULE.clustering$module_lazy(PGX = PGX, labeltype = labeltype))
  lazy_tabs <- c(lazy_tabs, MODULE.expression$module_lazy(PGX = PGX, labeltype = labeltype))
  lazy_tabs <- c(lazy_tabs, MODULE.enrichment$module_lazy(PGX = PGX, labeltype = labeltype, env = env))
  lazy_tabs <- c(lazy_tabs, MODULE.compare$module_lazy(PGX = PGX, labeltype = labeltype, auth = auth,
    env = env, reload_pgxdir = reload_pgxdir))
  lazy_tabs <- c(lazy_tabs, MODULE.systems$module_lazy(PGX = PGX))

  if (exists("MODULE.multiomics")) {
    lazy_tabs <- c(lazy_tabs, MODULE.multiomics$module_lazy(PGX = PGX))
  }
  if (exists("MODULE.wgcna")) {
    lazy_tabs <- c(lazy_tabs, MODULE.wgcna$module_lazy(PGX = PGX, save_pgx = env$save_pgx))
  }
  if (exists("MODULE.epigenomics")) {
    lazy_tabs <- c(lazy_tabs, MODULE.epigenomics$module_lazy(PGX = PGX))
  }

  ## Preload enrich-tab so selected_gsetmethods is available
  lazy_tabs[["enrich-tab"]]$preload <- TRUE

#  bigdash::bigTabsLazy(lazy_tabs)

  ## Modules needed after dataset is loaded (deferred) --------------
  observeEvent(env$load$is_data_loaded(), {

    # depending on datatpye, subset modules enabled and create modules active,
    ## if (tolower(PGX$datatype) == "multi-omics") {
    ##   MODULES_ACTIVE <- MODULES_MULTIOMICS
    ## } else if (tolower(PGX$datatype) == "methylomics") {
    ##   MODULES_ACTIVE <- MODULES_METHYLOMICS
    ## } else {
    ##   MODULES_ACTIVE <- MODULES_TRANSCRIPTOMICS
    ## }

    ## On the first dataset upload, we initialize the modules and
    ## UI. Preload modules get materialized here immediately, all
    ## others will only load/show when visited (i.e. lazy)
    if (env$load$is_data_loaded() == 1) {
      
      info("[SERVER] calling bigTabsLazy")
      shiny::withProgress(
        message = "Preparing your dashboard...",
        value = 0,
        {
          bigdash::bigTabsLazy(lazy_tabs)

          shinyjs::enable(selector = "a[data-value='Dashboard']")
          shinyjs::enable(selector = "a[data-value='Studio']")
          shinyjs::enable(selector = "a[data-value='Copilot']")
        })
        
      ## bigdash.hideMenuElement(session, "Clustering")
      ## bigdash.hideMenuElement(session, "Expression")
      ## bigdash.hideMenuElement(session, "GeneSets")
      ## bigdash.hideMenuElement(session, "Compare")
      ## bigdash.hideMenuElement(session, "SystemsBio")
      ## bigdash.hideMenuElement(session, "MultiOmics")
      ## bigdash.hideMenuElement(session, "WGCNA")
      ## bigdash.hideMenuElement(session, "Epigenomics")
      
    }
    
    ##MODULES_TO_REMOVE <- xor(MODULES_LOADED, MODULES_ACTIVE) & MODULES_LOADED
    ##MODULES_TO_LOAD <- xor(MODULES_LOADED, MODULES_ACTIVE) & MODULES_ACTIVE    
    ## lapply(names(MODULES_TO_REMOVE[MODULES_TO_REMOVE]), function(x) {
    ##   if (x == "DataView") {
    ##     bigdash.removeTab(session, "dataview-tab")
    ##     bigdash.hideMenuElement(session, "DataView")
    ##   }
    ##   if (x == "Clustering") {
    ##     lapply(names(MODULE.clustering$module_menu()), function(x) {
    ##       bigdash.removeTab(session, paste0(x, "-tab"))
    ##     })
    ##     bigdash.hideMenuElement(session, "Clustering")
    ##   }
    ##   if (x == "Expression") {
    ##     lapply(names(MODULE.expression$module_menu()), function(x) {
    ##       bigdash.removeTab(session, paste0(x, "-tab"))
    ##     })
    ##     bigdash.hideMenuElement(session, "Expression")
    ##   }
    ##   if (x == "GeneSets") {
    ##     lapply(names(MODULE.enrichment$module_menu()), function(x) {
    ##       bigdash.removeTab(session, paste0(x, "-tab"))
    ##     })
    ##     bigdash.hideMenuElement(session, "GeneSets")
    ##   }
    ##   if (x == "Compare") {
    ##     lapply(names(MODULE.compare$module_menu()), function(x) {
    ##       bigdash.removeTab(session, paste0(x, "-tab"))
    ##     })
    ##     bigdash.hideMenuElement(session, "Compare")
    ##   }
    ##   if (x == "SystemsBio") {
    ##     lapply(names(MODULE.systems$module_menu()), function(x) {
    ##       bigdash.removeTab(session, paste0(x, "-tab"))
    ##     })
    ##     bigdash.hideMenuElement(session, "SystemsBio")
    ##   }
    ##   if (x == "MultiOmics") {
    ##     lapply(names(MODULE.multiomics$module_menu()), function(x) {
    ##       bigdash.removeTab(session, paste0(x, "-tab"))
    ##     })
    ##     bigdash.hideMenuElement(session, "MultiOmics")
    ##   }
    ##   if (x == "WGCNA") {
    ##     lapply(names(MODULE.wgcna$module_menu()), function(x) {
    ##       bigdash.removeTab(session, paste0(x, "-tab"))
    ##     })
    ##     bigdash.hideMenuElement(session, "WGCNA")
    ##   }
    ##   if (x == "Epigenomics") {
    ##     lapply(names(MODULE.epigenomics$module_menu()), function(x) {
    ##       bigdash.removeTab(session, paste0(x, "-tab"))
    ##     })
    ##     bigdash.hideMenuElement(session, "Epigenomics")
    ##   }
    ## })
    ##
    ## if (env$load$is_data_loaded()) { # == 1) {
    ##   if (MODULES_TO_LOAD["DataView"]) {
    ##     bslib::nav_select("app-sidebar", selected = "Dashboard")
    ##   }
    ##   if (MODULES_TO_LOAD["Clustering"]) {
    ##     bigdash.showMenuElement(session, "Clustering")
    ##   }
    ##   if (MODULES_TO_LOAD["Expression"]) {
    ##     bigdash.showMenuElement(session, "Expression")
    ##   }
    ##   if (MODULES_TO_LOAD["GeneSets"]) {
    ##     bigdash.showMenuElement(session, "GeneSets")
    ##   }
    ##   if (MODULES_TO_LOAD["Compare"]) {
    ##     bigdash.showMenuElement(session, "Compare")
    ##   }
    ##   if (MODULES_TO_LOAD["SystemsBio"]) {
    ##     bigdash.showMenuElement(session, "SystemsBio")
    ##   }
    ##   if (MODULES_TO_LOAD["MultiOmics"] && exists("MODULE.multiomics")) {
    ##     bigdash.showMenuElement(session, "MultiOmics")
    ##   }
    ##   if (MODULES_TO_LOAD["WGCNA"] && exists("MODULE.wgcna")) {
    ##     bigdash.showMenuElement(session, "WGCNA")
    ##   }
    ##   if (MODULES_TO_LOAD["Epigenomics"] && exists("MODULE.epigenomics")) {
    ##     bigdash.showMenuElement(session, "Epigenomics")
    ##   }
    ## MODULES_LOADED <<- MODULES_ACTIVE

    ## initially shown. maybe better start with empty menu???
    all_tabs <- names(lazy_tabs)
    bigdash.filterTabs(session, all_tabs)
    
    if (env$load$is_data_loaded() > 0) {
      env$trigger_on_change_dataset(runif(1))
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

  ## tab_control is defined below; called on data load and after each tab
  ## navigation via bigdash::bigTabsLazy's internal observer.

  tab_control <- function() {

    show.beta <- env$user_settings$enable_beta()
    if (is.null(show.beta) || length(show.beta) == 0) show.beta <- FALSE

    has.libx <- dir.exists(file.path(OPG, "libx"))
    is.multiomics <- tolower(PGX$datatype) == "multi-omics"

    ## Collect all tab names from the currently active module set
    ## MODULES_MULTIOMICS/TRANSCRIPTOMICS/METHYLOMICS are logical arrays
    ## named by module name (e.g. "DataView", "Clustering", ...).
    if (is.multiomics) {
      active_mods <- names(MODULES_MULTIOMICS)[MODULES_MULTIOMICS]
    } else if (tolower(PGX$datatype) == "methylomics") {
      active_mods <- names(MODULES_METHYLOMICS)[MODULES_METHYLOMICS]
    } else {
      active_mods <- names(MODULES_TRANSCRIPTOMICS)[MODULES_TRANSCRIPTOMICS]
    }

    module_map <- list(
      DataView      = "dataview-tab",
      Clustering    = paste0(names(MODULE.clustering$module_menu()), "-tab"),
      Expression    = paste0(names(MODULE.expression$module_menu()), "-tab"),
      GeneSets      = paste0(names(MODULE.enrichment$module_menu()), "-tab"),
      Compare       = paste0(names(MODULE.compare$module_menu()), "-tab"),
      SystemsBio    = paste0(names(MODULE.systems$module_menu()), "-tab"),
      MultiOmics    = if (exists("MODULE.multiomics")) paste0(names(MODULE.multiomics$module_menu()), "-tab"),
      WGCNA         = if (exists("MODULE.wgcna")) paste0(names(MODULE.wgcna$module_menu()), "-tab"),
      Epigenomics   = if (exists("MODULE.epigenomics")) paste0(names(MODULE.epigenomics$module_menu()), "-tab")
    )

    tabs <- character(0)
    for (mod_name in names(module_map)) {
      if (mod_name %in% active_mods) {
        tabs <- c(tabs, module_map[[mod_name]])
      }
    }
    tabs <- unique(tabs)

    ## Beta visibility rules
    if (!(show.beta && has.libx))             tabs <- setdiff(tabs, "tcga-tab")
    if (!show.beta)                           tabs <- setdiff(tabs, "consensus-tab")
    if (!(isTRUE(opt$DEVMODE) && show.beta))  tabs <- setdiff(tabs, "preservation-tab")
    if (is.multiomics && !show.beta)          tabs <- setdiff(tabs, "wgcna-tab")
    if (!is.multiomics)                       tabs <- setdiff(tabs, "mwgcna-tab")

    ## Content requirements (replaces tabRequire)
    if (!rv_has_value(PGX, "drugs"))          tabs <- setdiff(tabs, "drug-tab")
    if (!rv_has_value(PGX, "wordcloud"))      tabs <- setdiff(tabs, "wordcloud-tab")
    if (!rv_has_value(PGX, "deconv"))         tabs <- setdiff(tabs, "cell-tab")
    if (!rv_has_value(PGX, "connectivity"))   tabs <- setdiff(tabs, "cmap-tab")
    time.vars <- playbase::get_timevars()
    found.time.var <- any(grepl(paste(time.vars, collapse = "|"), colnames(PGX$samples), ignore.case = TRUE))
    valid.contrasts <- any(grepl("IA:*", colnames(PGX$contrasts)))
    if (!(found.time.var && valid.contrasts))             tabs <- setdiff(tabs, "timeseries-tab")
    gset_avail <- rv_has_value(PGX, "gsetX") && rv_has_value(PGX, "gset.meta")
    if (!gset_avail) {
      tabs <- setdiff(tabs, c("enrich-tab", "pathway-tab", "isect-tab", "sig-tab"))
    }

    ## Datatype-specific filters
    if (PGX$datatype == "metabolomics") {
      tabs <- setdiff(tabs, "cmap-tab")
    }
    if (PGX$datatype == "multi-omics") {
      tabs <- setdiff(tabs, c("drug-tab", "cell-tab", "wordcloud-tab", "cmap-tab"))
    }
    if (!is.null(PGX$datatype) && tolower(PGX$datatype) != "methylomics") {
      tabs <- setdiff(tabs, "ideograms-tab")
    }

    ## Hide PCSF for methylomics DMP (CpG probe level — no meaningful PPI matching)
    if (!is.null(PGX$datatype) && tolower(PGX$datatype) == "methylomics") {
      is_dmp <- if (!is.null(PGX$dma)) {
        isTRUE(PGX$dma == "Differentially methylated positions")
      } else {
        mean(grepl("^cg[0-9]+", rownames(PGX$X))) > 0.5
      }
      if (is_dmp) {
        tabs <- setdiff(tabs, "pcsf-tab")
      }
    }

    if (is.null(PGX$datatype) || tolower(PGX$datatype) != "methylomics") {
      tabs <- setdiff(tabs, module_map[["Epigenomics"]] )
    }
    
    ## Apply filtering — hides tabs + sidebar items not in the allowed set
    bigdash.filterTabs(session, tabs)

    ## Hide Epigenomics sidebar menu group when not methylomics
    ## if (is.null(PGX$datatype) || tolower(PGX$datatype) != "methylomics") {
    ##   bigdash.hideMenuElement(session, "Epigenomics")
    ## }

    ## Subtab toggles inside boards — use try() since tabsetPanels may
    ## not exist yet. NOTE: should move to bigdash??
    try(toggleTab("drug-tabs", "Connectivity map (beta)", show.beta), silent = TRUE)
    try(toggleTab("pathway-tabs", "Enrichment Map (beta)", show.beta), silent = TRUE)
    has.customfc <- "custom" %in% colnames(PGX$gx.meta$meta[[1]]$fc) &&
      length(colnames(PGX$gx.meta$meta[[1]]$fc)) > 1
    try(toggleTab("diffexpr-tabs1", "FC-FC comparison", has.customfc), silent = TRUE)

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
