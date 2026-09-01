
## NOTE: This is not a real shiny module (yet...). We should move as
## much as possible OPG server related code here.

#' @param id Namespace id. Leave as the default `"app"` for the pre-existing,
#'   single OPG instance -- it runs completely unscoped (bare ids, no
#'   `moduleServer()` wrap), byte-identical to before this module conversion.
#'   Pass a distinct `id` (matching the same `id` given to [opg_ui()]) to run
#'   a second, independent OPG dashboard in the same session -- it is then a
#'   real nested `shiny::moduleServer()`, and every id it touches (its own
#'   tabs, and the board modules it mounts) is `"<id>-"`-prefixed via
#'   `session$ns()`, matching what `bigdash::bigPage(id = id, ...)` already
#'   generates in the UI. See `opg_ns()` in opg_ui.R.
#' @param menu_tree Which menu groups (and boards within them) this instance
#'   offers, in [opg_menu_tree()]'s shape. Must be the *same* tree passed to
#'   [opg_ui()] for this `id` -- only the groups named here have their boards
#'   registered/loaded at all (not just hidden), so passing a restricted tree
#'   (e.g. `opg_menu_tree()["MultiOmics"]`) both narrows the sidebar and
#'   avoids mounting boards this instance never shows.
#' @param parent_session The *root* app session. A handful of things are
#'   genuinely global rather than per-instance -- the outer top-level
#'   "app-sidebar" navset (Dashboard/Upload/Library panels) and the single
#'   app-wide "Label type" selector in the Settings panel (see
#'   board.user/appsettings_ui.R, deliberately un-namespaced) -- and must be
#'   read/written against the true root session even when this OPG instance
#'   is itself a nested module. Defaults to the current reactive domain, so
#'   existing `id = "app"` callers (where this *is* already the root session)
#'   need not pass it.
#' @param dashboard_nav_value The `parent_session`-level "app-sidebar" value
#'   that hosts this OPG instance, e.g. `"Dashboard"` for the primary
#'   instance. Used both to detect when this instance becomes visible (to
#'   offer the "no dataset loaded" popup) and to jump back to it once a
#'   dataset finishes loading.
#' @param load_example_dataset Companion reactiveVal to `load_example`
#'   (`LoadingBoard`'s own `load_example_dataset` param) -- set to
#'   `example_dataset` right before bumping `load_example()`, so
#'   `LoadingBoard` knows *which* dataset this instance's "Load example
#'   dataset" popup button means. `NULL` (the default) just bumps
#'   `load_example()` as before, and `LoadingBoard` falls back to
#'   "example-data".
#' @param example_dataset Dataset name (no `.pgx`) this instance's popup
#'   loads, e.g. `"mox-brca"` for a multiomics-only instance. Only takes
#'   effect when `load_example_dataset` is also given.
opg_server <- function(id, PGX, env, auth, reload_pgxdir, load_example = NULL,
                        load_example_dataset = NULL,
                        example_dataset = "example-data",
                        menu_tree = opg_menu_tree(),
                        parent_session = shiny::getDefaultReactiveDomain(),
                        dashboard_nav_value = "Dashboard") {

  allowed_groups <- names(menu_tree)
  ## Bare board ids (e.g. "pcsf") menu_tree lists, regardless of group --
  ## narrows a partly-included group (e.g. only "pcsf" of SystemsBio) down
  ## to just those boards. See opg_ui.R's opg_menu_boards()/filter_boards().
  allowed_boards <- opg_menu_boards(menu_tree)

  body <- function(input, output, session) {

  labeltype <- reactiveVal("feature") # can be feature (rownames counts), symbol or name

  ## Admin-picked boards for 'Basic menu' mode (opg_ui.R:opg_basic_menu_boards).
  ## A reactiveVal, not a plain vector: Admin Panel > Basic Menu
  ## (admin_server.R) updates it live via getUserOption() so the admin's own
  ## session re-filters immediately after Save, without needing a reload.
  ## Published through session$userData -- this app's existing cross-module
  ## channel (see setUserOption()/getUserOption() in utils.R, used the same
  ## way by AppSettingsBoard/copilot for AI settings) -- because
  ## AdminPanelBoard runs in its own module-scoped session, not this one.
  ##
  ## session$userData is shared across the *whole* browser session, not
  ## scoped per module -- so this wiring only ever targets the one, primary
  ## "app" instance the admin panel actually controls. A second instance's
  ## basic_mode_boards/locked_boards below still work (each gets its own
  ## reactiveVal, seeded from opg_basic_menu_boards()/opt$BASIC_LOCKED), they
  ## just are not admin-adjustable and would collide on these userData keys
  ## if also published.
  basic_mode_boards <- shiny::reactiveVal(opg_basic_menu_boards())
  locked_boards <- shiny::reactiveVal(opt$BASIC_LOCKED)
  if (identical(id, "app")) {
    setUserOption(session, "basic_menu_boards", basic_mode_boards) ## stores the reactiveVal itself, not its value
    setUserOption(session, "locked_boards", locked_boards)
  }

  shiny::observe({
    ids <- locked_boards()
    ## Clear everywhere first: a board can be *un*-locked as well.
    shinyjs::removeClass(selector = ".advanced-option-candidate", class = "advanced-option")
    for (b in ids) {
      shinyjs::addClass(selector = paste0(".advanced-option-candidate[data-board='", b, "']"),
                        class = "advanced-option")
    }
  })

  ## -------------------------------------------------------------
  ## No dataset loaded: offer the example dataset (like Qsee)
  ## -------------------------------------------------------------

  ## The Dashboard can be entered without any dataset (e.g. from the app
  ## launcher, which selects the nav panel programmatically and so bypasses
  ## the disabled nav link). Then all boards are empty, so ask the user
  ## whether to load the example dataset instead.
  ##
  ## "app-sidebar" is the *outer* app's own top-level navset (Dashboard /
  ## Upload / Library / ...), not part of this OPG instance's own bigPage()
  ## -- always read/write it against parent_session (== session itself for
  ## the default id = "app").
  shiny::observeEvent(parent_session$input[["app-sidebar"]], {
    shiny::req(parent_session$input[["app-sidebar"]] == dashboard_nav_value)
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
            session$ns("opg_load_example_from_popup"),
            "Load example dataset",
            class = "btn btn-outline-info welcome-btn-sm"
          ),
          shiny::actionButton(
            session$ns("opg_upload_new_from_popup"),
            "Upload new data",
            class = "btn btn-outline-info welcome-btn-sm"
          ),
          shiny::actionButton(
            session$ns("opg_load_library_from_popup"),
            "Load from library",
            class = "btn btn-outline-primary welcome-btn-sm"
          )
        ),
        ##footer = shiny::modalButton("Cancel"),
        footer = NULL,
        size = "s",
        easyClose = TRUE ## lets Escape (and click-outside) dismiss it
      )
    )
  })

  shiny::observeEvent(input$opg_load_example_from_popup, {
    shiny::removeModal()
    if (is.null(load_example)) {
      warning("[SERVER] !!! no load_example trigger available")
      return(NULL)
    }
    if (!is.null(load_example_dataset)) {
      load_example_dataset(example_dataset)
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
    bslib::nav_select("app-sidebar", selected = "Upload", session = parent_session)
  })

  shiny::observeEvent(input$opg_load_library_from_popup, {
    shiny::removeModal()
    info("[SERVER] opening library panel from popup")
    bslib::nav_select("app-sidebar", selected = "Library", session = parent_session)
  })

  ## Hide/show tabs. Open sidebar and settings
  shiny::observeEvent(
    {
      list(
        auth$logged,
        env$user_settings$enable_beta(),
        env$trigger_on_change_dataset(),
        parent_session$input$menu_basic,
        basic_mode_boards()
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
        bigdash.closeSidebar(session)
        bigdash.closeSettings(session)
        bigdash.selectTab(session, selected = session$ns("dataview-tab"))
        return(NULL)
      }

      ## show all main tabs
      bigdash.openSidebar(session)
      bigdash.openSettings(session = session)

    }
  )

  ## Toggle the basic-mode body class. The sidebar's own menu items are
  ## handled by tab_control() -> bigdash.filterTabs() (see the "Hide/show
  ## tabs" observer above, which lists input$menu_basic as a dependency). The
  ## body class drives the board-level simplifications (greyed
  ## .advanced-option blocks, hidden extra tabs) in scss/components/_app.scss
  ## -- CSS so it also applies to boards that are inserted lazily, after this
  ## observer has already run.
  observeEvent(parent_session$input$menu_basic, {
    if (isTRUE(parent_session$input$menu_basic)) {
      shinyjs::addClass(selector = "body", class = "basic-mode")
      ## The active subtab may be one we are about to hide. Only once the
      ## boards exist -- they are inserted when a dataset loads, and
      ## addressing a tabsetPanel that is not in the DOM yet just logs "there
      ## is no tabsetPanel with id ..." in the browser console.
      if (isTRUE(env$load$is_data_loaded() > 0)) {
        ## updateTabsetPanel() auto-applies session$ns() to inputId itself
        ## (like showTab()/hideTab() below) -- pass the bare local id, not a
        ## manually-namespaced one, or a nested instance double-prefixes.
        shiny::updateTabsetPanel(session, "dataview-tabs", selected = "Overview")
        shiny::updateTabsetPanel(session, "diffexpr-tabs1", selected = "Overview")
      }
    } else {
      shinyjs::removeClass(selector = "body", class = "basic-mode")
    }
  }, ignoreInit = FALSE)

  ## Register all lazy tabs at startup.  bigTabsLazy() will
  ## materialise each tab's heavy UI (via module_lazy closures)
  ## the first time the user clicks on it.
  ## ------------------------------------------------------------
  ## Keep only entries (bare-named, e.g. "clustersamples-tab") whose board
  ## is one this instance's menu_tree actually lists -- so a group that's
  ## only partly wanted (e.g. just "pcsf" out of all of SystemsBio) doesn't
  ## get its other boards mounted as Shiny modules at all, not just hidden.
  ## Matches opg_ui.R's filter_boards() for the UI-shell half.
  filter_lazy <- function(lazy_list) {
    lazy_list[sub("-tab$", "", names(lazy_list)) %in% allowed_boards]
  }

  lazy_tabs <- list()

  if ("dataview" %in% allowed_boards) {
    lazy_tabs[["dataview-tab"]] <- list(
      ui      = function() DataViewUI(session$ns("dataview")),
      server  = function() DataViewBoard("dataview", pgx = PGX, labeltype = labeltype),
      preload = TRUE
    )
  }

  if ("diffexpr" %in% allowed_boards) {
    lazy_tabs[["diffexpr-tab"]] <- list(
      ui      = function() ExpressionUI(session$ns("diffexpr")),
      server  = function() ExpressionBoard("diffexpr", pgx = PGX, labeltype = labeltype) ->>
                             env$diffexpr,
      preload = TRUE
    )
  }

  if ("enrich" %in% allowed_boards) {
    lazy_tabs[["enrich-tab"]] = list(
      ui     = function() EnrichmentUI(session$ns("enrich")),
      server = function() EnrichmentBoard("enrich", pgx = PGX, selected_gxmethods = env$diffexpr$selected_gxmethods) ->>
                            env$enrich,
      preload = TRUE
    )
  }

  lazy_tabs <- c(lazy_tabs, filter_lazy(MODULE.clustering$module_lazy(PGX = PGX, labeltype = labeltype, ns = session$ns)))
  lazy_tabs <- c(lazy_tabs, filter_lazy(MODULE.expression$module_lazy(PGX = PGX, labeltype = labeltype, ns = session$ns)))
  lazy_tabs <- c(lazy_tabs, filter_lazy(MODULE.enrichment$module_lazy(PGX = PGX, labeltype = labeltype, env = env, ns = session$ns)))
  lazy_tabs <- c(lazy_tabs, filter_lazy(MODULE.compare$module_lazy(PGX = PGX, labeltype = labeltype, auth = auth,
    env = env, reload_pgxdir = reload_pgxdir, ns = session$ns)))
  lazy_tabs <- c(lazy_tabs, filter_lazy(MODULE.systems$module_lazy(PGX = PGX, ns = session$ns)))

  if (exists("MODULE.multiomics")) {
    lazy_tabs <- c(lazy_tabs, filter_lazy(MODULE.multiomics$module_lazy(PGX = PGX, ns = session$ns)))
  }
  if (exists("MODULE.wgcna")) {
    lazy_tabs <- c(lazy_tabs, filter_lazy(MODULE.wgcna$module_lazy(PGX = PGX, save_pgx = env$save_pgx, ns = session$ns)))
  }
  if (exists("MODULE.epigenomics")) {
    lazy_tabs <- c(lazy_tabs, filter_lazy(MODULE.epigenomics$module_lazy(PGX = PGX, ns = session$ns)))
  }

  ## Every group above returns its lazy_tabs entries bare-named (e.g.
  ## "clustersamples-tab"); namespace them all here, once, exactly like
  ## qsee_server.R does -- matches the fully-namespaced bigTabItem() names
  ## opg_ui.R generated for this same `id` (via its own opg_ns()/session$ns()
  ## equivalent), which bigTabsLazy() requires.
  names(lazy_tabs) <- session$ns(names(lazy_tabs))

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
            lazy_tabs,
            id = id
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
    bslib::nav_select("app-sidebar", selected = dashboard_nav_value, session = parent_session)
    bigdash.openSettings(lock = TRUE, session = session)
    bigdash.openSidebar(session)
    bigdash.showTabs(session)
    bigdash.selectTab(session, session$ns("dataview-tab"))

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

    ## Restrict to whatever menu_tree/allowed_groups this OPG instance was
    ## built with (opg_server(..., menu_tree = )) -- e.g. a second instance
    ## mounted with only the "MultiOmics" group never shows (or loads) any
    ## other board, regardless of dataset type.
    active_mods <- intersect(active_mods, allowed_groups)

    module_map <- list(
      DataView      = session$ns("dataview-tab"),
      Clustering    = session$ns(paste0(names(MODULE.clustering$module_menu()), "-tab")),
      Expression    = session$ns(paste0(names(MODULE.expression$module_menu()), "-tab")),
      GeneSets      = session$ns(paste0(names(MODULE.enrichment$module_menu()), "-tab")),
      Compare       = session$ns(paste0(names(MODULE.compare$module_menu()), "-tab")),
      SystemsBio    = session$ns(paste0(names(MODULE.systems$module_menu()), "-tab")),
      MultiOmics    = if (exists("MODULE.multiomics")) session$ns(paste0(names(MODULE.multiomics$module_menu()), "-tab")),
      WGCNA         = if (exists("MODULE.wgcna")) session$ns(paste0(names(MODULE.wgcna$module_menu()), "-tab")),
      Epigenomics   = if (exists("MODULE.epigenomics")) session$ns(paste0(names(MODULE.epigenomics$module_menu()), "-tab"))
    )

    tabs <- character(0)
    for (mod_name in names(module_map)) {
      if (mod_name %in% active_mods) {
        tabs <- c(tabs, module_map[[mod_name]])
      }
    }
    tabs <- unique(tabs)

    ## Restrict to exactly the boards this instance's menu_tree lists --
    ## module_map above is built per whole group, so this is what narrows
    ## e.g. an "active" SystemsBio group down to just "pcsf" when that's
    ## all menu_tree included for it. See allowed_boards/filter_lazy() for
    ## the matching server-side board registration restriction.
    tabs <- intersect(tabs, session$ns(paste0(allowed_boards, "-tab")))

    ## Beta visibility rules
    if (!(show.beta && has.libx))             tabs <- setdiff(tabs, session$ns("tcga-tab"))
    if (!show.beta)                           tabs <- setdiff(tabs, session$ns("consensus-tab"))
    if (!(isTRUE(opt$DEVMODE) && show.beta))  tabs <- setdiff(tabs, session$ns("preservation-tab"))
    if (is.multiomics && !show.beta)          tabs <- setdiff(tabs, session$ns("wgcna-tab"))
    if (!is.multiomics)                       tabs <- setdiff(tabs, session$ns("mwgcna-tab"))

    ## Content requirements (replaces tabRequire)
    if (!rv_has_value(PGX, "drugs"))          tabs <- setdiff(tabs, session$ns("drug-tab"))
    if (!rv_has_value(PGX, "wordcloud"))      tabs <- setdiff(tabs, session$ns("wordcloud-tab"))
    if (!rv_has_value(PGX, "deconv"))         tabs <- setdiff(tabs, session$ns("cell-tab"))
    if (!rv_has_value(PGX, "connectivity"))   tabs <- setdiff(tabs, session$ns("cmap-tab"))
    time.vars <- playbase::get_timevars()
    found.time.var <- any(grepl(paste(time.vars, collapse = "|"), colnames(PGX$samples), ignore.case = TRUE))
    valid.contrasts <- any(grepl("IA:*", colnames(PGX$contrasts)))
    if (!(found.time.var && valid.contrasts))             tabs <- setdiff(tabs, session$ns("timeseries-tab"))
    gset_avail <- rv_has_value(PGX, "gsetX") && rv_has_value(PGX, "gset.meta")
    if (!gset_avail) {
      tabs <- setdiff(tabs, session$ns(c("enrich-tab", "pathway-tab", "isect-tab", "sig-tab")))
    }

    ## Datatype-specific filters
    if (PGX$datatype == "metabolomics") {
      tabs <- setdiff(tabs, session$ns("cmap-tab"))
    }
    if (PGX$datatype == "multi-omics") {
      tabs <- setdiff(tabs, session$ns(c("drug-tab", "cell-tab", "wordcloud-tab", "cmap-tab")))
    }
    if (!is.null(PGX$datatype) && tolower(PGX$datatype) != "methylomics") {
      tabs <- setdiff(tabs, session$ns("ideograms-tab"))
    }

    ## Hide PCSF for methylomics DMP (CpG probe level — no meaningful PPI matching)
    if (!is.null(PGX$datatype) && tolower(PGX$datatype) == "methylomics") {
      is_dmp <- if (!is.null(PGX$dma)) {
        isTRUE(PGX$dma == "Differentially methylated positions")
      } else {
        mean(grepl("^cg[0-9]+", rownames(PGX$X))) > 0.5
      }
      if (is_dmp) {
        tabs <- setdiff(tabs, session$ns("pcsf-tab"))
      }
    }

    if (is.null(PGX$datatype) || tolower(PGX$datatype) != "methylomics") {
      tabs <- setdiff(tabs, module_map[["Epigenomics"]] )
    }

    ## Basic mode: further restrict to what the admin picked in
    ## Admin panel > Basic menu (opg_ui.R:opg_basic_menu_boards). Admin only
    ## ever controls the primary "app" instance (see basic_mode_boards()
    ## above) and its picks are bare board ids, not namespaced -- applying
    ## them to a second instance's already-namespaced `tabs` would just
    ## intersect to nothing and hide everything, so skip there.
    if (identical(id, "app") && isTRUE(parent_session$input$menu_basic)) {
      tabs <- intersect(tabs, basic_mode_boards())
    }

    ## Apply filtering — hides tabs + sidebar items not in the allowed set
    bigdash.filterTabs(session, tabs, id = id)

    ## Subtab toggles inside boards. showTab()/hideTab() fail on the *client*
    ## ("there is no tabsetPanel with id ..."), where a server-side try() cannot
    ## catch them -- so ask which boards exist before driving their tabsets
    ## rather than firing at all of them and swallowing the fallout.
    ##
    ## showTab()/hideTab() (inside toggleTab()) auto-apply session$ns() to
    ## their inputId themselves via the ambient reactive domain -- pass bare
    ## local ids here, matching the updateTabsetPanel() calls above, or a
    ## nested instance double-prefixes and neither call finds its tabsetPanel.
    loaded <- loaded_tabs()
    if (session$ns("drug-tab") %in% loaded) {
      toggleTab("drug-tabs", "Connectivity map (beta)", show.beta, session = session)
    }
    if (session$ns("pathway-tab") %in% loaded) {
      toggleTab("pathway-tabs", "Enrichment Map (beta)", show.beta, session = session)
    }
    if (session$ns("diffexpr-tab") %in% loaded) {
      has.customfc <- "custom" %in% colnames(PGX$gx.meta$meta[[1]]$fc) &&
        length(colnames(PGX$gx.meta$meta[[1]]$fc)) > 1
      toggleTab("diffexpr-tabs1", "FC-FC comparison", has.customfc, session = session)
    }

  }

  ## -------------------------------------------------------------
  ## Labeltype stuff
  ## -------------------------------------------------------------

  ## "selected_labeltype" is the single, app-wide "Label type" selector
  ## in the Settings panel (board.user/appsettings_ui.R) -- deliberately
  ## defined there without ns(), i.e. it is not per-OPG-instance. Read/write
  ## it against parent_session (== session itself for id = "app"); the local
  ## `labeltype` reactiveVal above is what's actually threaded down to this
  ## instance's own boards, so multiple instances still each get their own
  ## copy even though they all react to the one shared dropdown.

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
        parent_session,
        "selected_labeltype",
        choices = label_types,
        selected = sel.labeltype
      )
    }
  )

  # change label type based on selected input
  shiny::observeEvent(
    {
      parent_session$input$selected_labeltype
    },
    {
      labeltype(parent_session$input$selected_labeltype)
      if (!is.null(PGX$genes)) {
        lab <- parent_session$input$selected_labeltype
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

  } ## end of body()

  if (identical(id, "app")) {
    ## The pre-existing, single-instance OPG: run directly against the root
    ## session, exactly as before this file gained module support -- no
    ## shiny::moduleServer() wrap, so every id stays bare (opg_ns()/session$ns()
    ## are identity in this branch, matching bigdash's own BIGDASH_DEFAULT_ID
    ## convention for bigPage(id = "app", ...)).
    session <- parent_session
    input <- session$input
    output <- session$output
    body(input, output, session)
  } else {
    shiny::moduleServer(id, body)
  }
}
