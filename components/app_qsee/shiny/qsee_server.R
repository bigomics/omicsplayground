##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

qsee_server <- function(id, pgx=NULL, matx="counts", parent=NULL, purge=NULL, lazy=TRUE) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      message("[qsee_server:board] called ")

      ## Resolved once and threaded down to every board's own bd_is_visible()
      ## call (see qsee_resolve_purge()): NULL (the default) keeps whatever
      ## ENABLE_PLOTLY_PURGE in etc/OPTIONS says; an explicit TRUE/FALSE here
      ## overrides it for the whole module, e.g. to A/B compare purge on vs
      ## off from the test harness without touching etc/OPTIONS.
      purge_enabled <- qsee_resolve_purge(purge)

      ## Holds user-uploaded matrices from the upload modal (NULL = use defaults)
      uploaded <- shiny::reactiveValues(
        X = NULL,
        Y = NULL,
        K = NULL
      )

      ## Track if we've already shown the "no dataset" popup (to avoid spamming)
      has_shown_popup <- shiny::reactiveVal(FALSE)

      ## bigdash only auto-selects the first tab of a bigPage() that is
      ## already visible at the moment Shiny connects. `bigdash.openSidebar()`
      ## / `bigdash.openSettings()` now scope themselves to the calling
      ## module's own bigPage() (via `session`), so they no longer no-op or
      ## leak into other bigPage() instances on the page for a nested, scoped
      ## bigPage() like this one -- requires bigdash >= the version adding
      ## scoped nav helpers (PR bigomics/bigdash#feat/scoped-sidebar-settings-nav).
      ## Once Qsee is first shown: select its first tab and force its own
      ## left menu / right settings panel open.
      ## No plots are rendered directly under the top-level Qsee module, so
      ## nothing is ever purged on its behalf -- purge = FALSE.
      is_visible <- bigdash::bd_is_visible(input, purge = FALSE, label="BOARD")
      
      initial_tab_selected <- shiny::reactiveVal(FALSE)
      shiny::observeEvent(is_visible(), {
        shiny::req(isTRUE(is_visible()))
        if (!shiny::isolate(initial_tab_selected())) {
          bigdash.selectTab(session, session$ns("normalize-tab"))
          bigdash.openSidebar(session)
          bigdash.openSettings(session = session)
          initial_tab_selected(TRUE)
        }
      })

      shiny::observeEvent( list(pgx$X, pgx$samples), {
        if(matx == "counts") {
          prior <- quantile(pgx$counts[pgx$counts>0], probs=0.001, na.rm=TRUE)[1]
          uploaded$X <- log2(pgx$counts + prior)
        } else {
          uploaded$X <- pgx$X
        }
        uploaded$Y <- pgx$samples
      })

      show_upload_modal <- function() {
        shiny::showModal(shiny::modalDialog(
          title = "Upload CSV files",
          shiny::fileInput(
            session$ns("counts_csv"),
            "Counts (counts.csv)",
            accept = c(".csv", "text/csv", "text/comma-separated-values")
          ),
          shiny::fileInput(
            session$ns("samples_csv"),
            "Samples (samples.csv)",
            accept = c(".csv", "text/csv", "text/comma-separated-values")
          ),
          easyClose = TRUE,
          footer = shiny::modalButton("Close")
        ))
      }

      shiny::observeEvent(input$upload, {
        show_upload_modal()
      })

      shiny::observeEvent(input$load_example, {
        example_pgx <- playdata::GEIGER_PGX
        uploaded$X <- example_pgx$X
        uploaded$Y <- example_pgx$samples[,-1]
        message("[qsee_server:board] loaded example data")
      })

      shiny::observeEvent(input$load_pgx, {
        if(is.null(pgx$X)) {
          shinyalert::shinyalert(text = "no dataset uploaded", type = "error")
          return(NULL)
        }
        uploaded$X <- pgx$X
        uploaded$Y <- pgx$samples
        message("[qsee_server:board] loaded current pgx")
      })
      
      ## Handle popup buttons
      shiny::observeEvent(input$load_example_from_popup, {
        shiny::removeModal()
        pgx <- playdata::GEIGER_PGX
        if(matx == "counts") {
          prior <- quantile(pgx$counts[pgx$counts>0], probs=0.001, na.rm=TRUE)[1]
          uploaded$X <- log2(pgx$counts + prior)
        } else {
          uploaded$X <- pgx$X          
        }
        uploaded$Y <- pgx$samples[,-1]
        message("[qsee_server:board] loaded example data from popup")
      })

      shiny::observeEvent(input$counts_csv, {
        shiny::req(input$counts_csv)
        df <- playbase::read_counts(input$counts_csv$datapath)
        df <- as.matrix(df)
        if(!playbase::is_logged(df)) {
          prior <- quantile(df[df>0], probs=0.001, na.rm=TRUE)[1]          
          df <- log2(df + prior)
        }
        uploaded$X <- df
      })

      shiny::observeEvent(input$samples_csv, {
        shiny::req(input$samples_csv)
        df <- playbase::read_samples(input$samples_csv$datapath)        
        uploaded$Y <- df
      })
      
      ## Clean expression matrix. Noise injection for the Normalization
      ## board is controlled there via the "Noise amount" slider.
      getX <- reactive({
        X <- uploaded$X 
        X[which(X == 0)] <- NA
        if(!playbase::is_logged(X)) X <- log2(X)
        X
      })

      getY <- reactive({
        uploaded$Y
      })
      
      ## lazy = TRUE (default): bigTabsLazy() inserts each board's UI (built
      ## in qsee_ui.R's bigTabItem(), see comment there) and calls its server
      ## the first time that tab is opened, instead of paying for all six
      ## boards' moduleServer setup unconditionally at Qsee startup.
      ## lazy = FALSE calls every board server directly here instead, up
      ## front -- the pre-bigTabsLazy behaviour, kept switchable for an A/B
      ## comparison of the two. Must match the `lazy` given to qsee_ui() for
      ## the same `id`: qsee_ui(lazy=FALSE) already builds every board's UI
      ## up front, so bigTabsLazy() would find nothing left to insert.
      if (isTRUE(lazy)) {
        lazy_tabs <- list(
          "normalize-tab" = list(
            ui     = function() qsee_normalization_ui(session$ns("normalize")),
            server = function() qsee_normalization_server("normalize", rX=getX, rY=getY, purge=purge_enabled)
          ),
          "impute-tab" = list(
            ui     = function() qsee_imputation_ui(session$ns("impute")),
            server = function() qsee_imputation_server("impute", rX=getX, rY=getY, purge=purge_enabled)
          ),
          "bsee-tab" = list(
            ui     = function() qsee_bsee_ui(session$ns("bsee")),
            server = function() qsee_bsee_server("bsee", rX=getX, rY=getY, purge=purge_enabled)
          ),
          "outlier-tab" = list(
            ui     = function() qsee_outlier_ui(session$ns("outlier")),
            server = function() qsee_outlier_server("outlier", rX=getX, rY=getY, purge=purge_enabled)
          ),
          "filtering-tab" = list(
            ui     = function() qsee_filtering_ui(session$ns("filtering")),
            server = function() qsee_filtering_server("filtering", rX=getX, rY=getY, purge=purge_enabled)
          ),
          "pcaexplorer-tab" = list(
            ui     = function() qsee_pcaexplorer_ui(session$ns("pcaexplorer")),
            server = function() qsee_pcaexplorer_server("pcaexplorer", rX=getX, rY=getY, purge=purge_enabled)
          )
        )
        ## bigTabsLazy() matches tabs by the fully-namespaced name
        ## bigTabItem() was given in qsee_ui.R (e.g. "qsee-normalize-tab").
        ## Namespace these from the list's own names rather than a
        ## separately hand-kept vector, so reordering an entry above can't
        ## silently desync the two.
        names(lazy_tabs) <- session$ns(names(lazy_tabs))
        loaded <- bigdash::bigTabsLazy(lazy_tabs, id = id)
      } else {
        qsee_normalization_server("normalize", rX=getX, rY=getY, purge=purge_enabled)
        qsee_imputation_server("impute", rX=getX, rY=getY, purge=purge_enabled)
        qsee_bsee_server("bsee", rX=getX, rY=getY, purge=purge_enabled)
        qsee_outlier_server("outlier", rX=getX, rY=getY, purge=purge_enabled)
        qsee_filtering_server("filtering", rX=getX, rY=getY, purge=purge_enabled)
        qsee_pcaexplorer_server("pcaexplorer", rX=getX, rY=getY, purge=purge_enabled)
      }

      ## If no dataset is loaded (X is NULL/empty), show a helpful popup only
      ## after the Qsee module itself is visible.  In the main application Qsee
      ## is a hidden parent nav panel during startup.
      ##
      ## Depend on is_visible() alone, not on input$nav: nav only changes once
      ## bigdash.selectTab() above completes its own client round trip, and
      ## chaining the popup off *that* instead of the visibility probe just
      ## adds a redundant extra hop before the popup can appear.
      shiny::observeEvent(is_visible(), {

        X <- uploaded$X
        message("[qsee_server:board] isnullX = ", is.null(X))
        message("[qsee_server:board] isvisible = ", is_visible())
        message("[qsee_server:board] has_shown_popup = ", has_shown_popup())
        
        shiny::req(isTRUE(is_visible()))
        noX <- is.null(X) || length(X) == 0 || nrow(X) == 0 || ncol(X) == 0
        if (noX && !has_shown_popup()) {
          has_shown_popup(TRUE)
          shiny::showModal(
            shiny::modalDialog(
              title = "No dataset loaded",
              shiny::p(
                "No dataset has been loaded yet. What would you like to do?"
              ),
              div(
                style = "text-align:center;",
                shiny::actionButton(
                  session$ns("load_example_from_popup"),
                  "Load example dataset",
                  class = "btn btn-outline-info welcome-btn-sm"
                ),
                shiny::actionButton(
                  session$ns("upload_new_from_popup"),
                  "Upload new data",
                  class = "btn btn-outline-info welcome-btn-sm"
                ),
                shiny::actionButton(
                  session$ns("load_library_from_popup"),
                  "Load from library",
                  class = "btn btn-outline-primary welcome-btn-sm"
                )
              ),
              footer = shiny::modalButton("Cancel"),
              size = "s",
              easyClose = FALSE
            )
          )
        }
      })

      shiny::observeEvent(input$upload_new_from_popup, {
        shiny::removeModal()
        show_upload_modal()
      })

      shiny::observeEvent(input$load_library_from_popup, {
        shiny::removeModal()
        if (is.null(parent)) {
          warning("[qsee_server] !!! no parent session available for Library navigation")
          return(NULL)
        }
        bslib::nav_select("app-sidebar", "Library", session = parent)
      })

    } ## end-of-server
  )
}
