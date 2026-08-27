
## NOTE: This is not a real shiny module (yet...). We should move as
## much as possible OPG server related code here.

opg_server <- function(id, input, output, session, PGX, env, auth, reload_pgxdir,
                       load_example = NULL) {

  labeltype <- reactiveVal("feature") # can be feature (rownames counts), symbol or name
  
  if(id != "app") {
    stop("FATAL: opg_server is not a proper ShinyModule yet")
  }
  
  ## -------------------------------------------------------------
  ## No dataset loaded: offer the example dataset (like Qsee)
  ## -------------------------------------------------------------

  ## The Dashboard can be entered without any dataset (e.g. from the app
  ## launcher, which selects the nav panel programmatically and so bypasses
  ## the disabled nav link). Then all boards are empty, so ask the user
  ## whether to load the example dataset instead.
  shiny::observeEvent(input[["app-sidebar"]], {
    shiny::req(input[["app-sidebar"]] == "Dashboard")
    shiny::req(isTRUE(auth$logged))
    noX <- is.null(PGX$X) || length(PGX$X) == 0
    if (!noX) return(NULL)
    info("[SERVER] no dataset loaded: asking for example data")
    shiny::showModal(
      shiny::modalDialog(
        title = "No dataset loaded",
        shiny::p(
          "No dataset has been loaded yet. What would you like to do?"
        ),
        div(
          style = "text-align:center;",
          shiny::actionButton(
            "opg_load_example_from_popup",
            "Load example dataset",
            class = "btn btn-outline-info welcome-btn-sm"
          ),
          shiny::actionButton(
            "opg_upload_new_from_popup",
            "Upload new data",
            class = "btn btn-outline-info welcome-btn-sm"
          ),
          shiny::actionButton(
            "opg_load_library_from_popup",
            "Load from library",
            class = "btn btn-outline-primary welcome-btn-sm"
          )
        ),
        ##footer = shiny::modalButton("Cancel"),
        footer = NULL,
        size = "s",
        easyClose = FALSE
      )
    )
  })

  shiny::observeEvent(input$opg_load_example_from_popup, {
    shiny::removeModal()
    if (is.null(load_example)) {
      warning("[SERVER] !!! no load_example trigger available")
      return(NULL)
    }
    if (is.null(load_example())) {
      load_example(1)
    } else {
      load_example(load_example() + 1)
    }
    info("[SERVER] loading example data from popup")
  })

  shiny::observeEvent(input$opg_upload_new_from_popup, {
    shiny::removeModal()
    info("[SERVER] opening upload panel from popup")
    bslib::nav_select("app-sidebar", selected = "Upload", session = session)
  })

  shiny::observeEvent(input$opg_load_library_from_popup, {
    shiny::removeModal()
    info("[SERVER] opening library panel from popup")
    bslib::nav_select("app-sidebar", selected = "Library", session = session)
  })

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

  ## Hide/show basic or full menu. The body class drives the board-level
  ## simplifications (greyed .advanced-option blocks, hidden extra tabs) in
  ## scss/components/_app.scss -- CSS so it also applies to boards that are
  ## inserted lazily, after this observer has already run.
  observeEvent(input$menu_basic, {
    if (isTRUE(input$menu_basic)) {
      shinyjs::runjs("$('#menu-basic').removeClass('nodisp').show(); $('#menu-full').addClass('nodisp').hide();")
      shinyjs::addClass(selector = "body", class = "basic-mode")
      ## The active tab may be one we are about to hide. Only once the boards
      ## exist -- they are inserted when a dataset loads, and addressing a
      ## tabsetPanel that is not in the DOM yet just logs "there is no
      ## tabsetPanel with id ..." in the browser console.
      if (env$load$is_data_loaded() > 0) {
        shiny::updateTabsetPanel(session, "dataview-tabs", selected = "Overview")
        shiny::updateTabsetPanel(session, "diffexpr-tabs1", selected = "Overview")
      }
    } else {
      shinyjs::runjs("$('#menu-full').removeClass('nodisp').show(); $('#menu-basic').addClass('nodisp').hide();")
      shinyjs::removeClass(selector = "body", class = "basic-mode")
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

  ## Set by bigTabsLazy() below; reports which tabs have been materialised so
  ## far, so tab_control() only drives widgets that actually exist yet.
  loaded_tabs <- function() character(0)

  ## Each tab's shell (opg_ui) shows a spinner until its board arrives, and
  ## tab_control()'s subtab toggles can only reach a board's tabsetPanel once
  ## that board exists. Both are per-load concerns, so they hang off the server
  ## closure: bigTabsLazy() inserts the UI first, then calls this.
  with_tab_loaded <- function(name, server_fun) {
    force(name)
    force(server_fun)
    function() {
      if (is.function(server_fun)) server_fun()
      shinyjs::hide(selector = paste0("[id='", sub("-tab$", "-loader", name), "']"))
      ## try(): a board that fails to build must not take the whole tab set
      ## with it, and tab_control() reads a lot of PGX.
      try(tab_control(), silent = TRUE)
    }
  }
  lazy_tabs <- Map(
    function(name, spec) {
      spec$server <- with_tab_loaded(name, spec$server)
      spec
    },
    names(lazy_tabs), lazy_tabs
  )

  ## Modules needed after dataset is loaded (deferred) --------------
  observeEvent(env$load$is_data_loaded(), {

    ## On the first dataset upload, we initialize the modules and
    ## UI. Preload modules get materialized here immediately, all
    ## others will only load/show when visited (i.e. lazy)
    if (env$load$is_data_loaded() == 1) {
      
      info("[SERVER] calling bigTabsLazy")
      shiny::withProgress(
        message = "Preparing your dashboard...",
        value = 0,
        {
          loaded_tabs <<- bigdash::bigTabsLazy(
            lazy_tabs
          )

          shinyjs::enable(selector = "a[data-value='Dashboard']")
          shinyjs::enable(selector = "a[data-value='Studio']")
          shinyjs::enable(selector = "a[data-value='Copilot']")
        })
    }

    ## No "show everything" pass here: tab_control(), reached through the
    ## trigger below, computes the allowed set from scratch and filterTabs()
    ## hides whatever is not in it. Showing all 28 tabs first only flashed the
    ## ones about to be hidden again.
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

    ## Subtab toggles inside boards. showTab()/hideTab() fail on the *client*
    ## ("there is no tabsetPanel with id ..."), where a server-side try() cannot
    ## catch them -- so ask which boards exist before driving their tabsets
    ## rather than firing at all of them and swallowing the fallout.
    loaded <- loaded_tabs()
    if ("drug-tab" %in% loaded) {
      toggleTab("drug-tabs", "Connectivity map (beta)", show.beta)
    }
    if ("pathway-tab" %in% loaded) {
      toggleTab("pathway-tabs", "Enrichment Map (beta)", show.beta)
    }
    if ("diffexpr-tab" %in% loaded) {
      has.customfc <- "custom" %in% colnames(PGX$gx.meta$meta[[1]]$fc) &&
        length(colnames(PGX$gx.meta$meta[[1]]$fc)) > 1
      toggleTab("diffexpr-tabs1", "FC-FC comparison", has.customfc)
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
