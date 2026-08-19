##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' Collapsible developer-options block for QSEE settings sidebars.
#'
#' Borderless and transparent so it inherits the settings pane (aliceblue)
#' instead of looking like a Bootstrap card.
qsee_developer_options <- function(id, ...) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::tags$style(HTML("
      .qsee-dev-options {
        --bs-accordion-bg: transparent;
        --bs-accordion-active-bg: transparent;
        --bs-accordion-btn-bg: transparent;
        --bs-accordion-border-width: 0;
        --bs-accordion-border-color: transparent;
        --bs-accordion-btn-focus-box-shadow: none;
        --bs-accordion-inner-border-radius: 0;
        --bs-accordion-border-radius: 0;
        --bs-accordion-btn-icon: url('data:image/svg+xml,%3csvg xmlns=%27http://www.w3.org/2000/svg%27 viewBox=%270 0 16 16%27 fill=%27%23dc3545%27%3e%3cpath fill-rule=%27evenodd%27 d=%27M1.646 4.646a.5.5 0 0 1 .708 0L8 10.293l5.646-5.647a.5.5 0 0 1 .708.708l-6 6a.5.5 0 0 1-.708 0l-6-6a.5.5 0 0 1 0-.708z%27/%3e%3c/svg%3e');
        --bs-accordion-btn-active-icon: url('data:image/svg+xml,%3csvg xmlns=%27http://www.w3.org/2000/svg%27 viewBox=%270 0 16 16%27 fill=%27%23dc3545%27%3e%3cpath fill-rule=%27evenodd%27 d=%27M1.646 4.646a.5.5 0 0 1 .708 0L8 10.293l5.646-5.647a.5.5 0 0 1 .708.708l-6 6a.5.5 0 0 1-.708 0l-6-6a.5.5 0 0 1 0-.708z%27/%3e%3c/svg%3e');
      }
      .qsee-dev-options .accordion-item,
      .qsee-dev-options .accordion-header,
      .qsee-dev-options .accordion-button,
      .qsee-dev-options .accordion-button:not(.collapsed),
      .qsee-dev-options .accordion-collapse,
      .qsee-dev-options .accordion-body {
        background-color: transparent !important;
        background: transparent !important;
        border: 0 !important;
        box-shadow: none !important;
      }
      .qsee-dev-options .accordion-button:focus {
        box-shadow: none !important;
        border-color: transparent !important;
      }
      .qsee-dev-options .accordion-button,
      .qsee-dev-options .accordion-body {
        padding-left: 0;
        padding-right: 0;
      }
      /* Bottom-aligned: closed points up (open toward content), open points down. */
      .qsee-dev-options .accordion-button::after {
        transform: rotate(-180deg);
      }
      .qsee-dev-options .accordion-button:not(.collapsed)::after {
        transform: rotate(0deg);
      }
    ")),
    bslib::accordion(
      id = ns("dev_options"),
      open = FALSE,
      class = "qsee-dev-options",
      bslib::accordion_panel(
        value = "developer-options",
        title = shiny::span("Developer options", class = "text-danger"),
        ...
      )
    )
  )
}

#' Standalone page wrapper for a QSEE (or other bigdash) applet.
#'
#' Fills the viewport (min-height 800px), loads Playground styles + Lato,
#' and attaches shinyjs / bigLoaders deps used by PlotModule.
qsee_app_shell <- function(..., theme = NULL, min_height = 800) {
  if (is.null(theme)) {
    theme <- bslib::bs_theme(base_font = bslib::font_google("Lato"))
  }
  mh <- as.integer(min_height)
  bslib::page_fillable(
    theme = theme,
    padding = 0,
    gap = 0,
    shiny::tags$head(
      shiny::tags$link(rel = "stylesheet", href = "custom/styles.min.css"),
      shiny::tags$link(
        rel = "stylesheet",
        href = "https://fonts.googleapis.com/css2?family=Lato:ital,wght@0,300;0,400;0,700;1,400&display=swap"
      ),
      shiny::tags$style(HTML(sprintf("
        :root {
          --bs-font-sans-serif: 'Lato', sans-serif;
          --bs-body-font-family: 'Lato', sans-serif;
        }
        html, body {
          font-family: 'Lato', sans-serif;
          height: 100%%;
          min-height: %dpx;
          margin: 0;
          overflow: hidden;
        }
        @media (max-height: %dpx) {
          html, body { overflow-y: auto; }
        }
        .bslib-page-fill,
        .bigdash-app,
        .big-full-page,
        .bigdash-sidebar-shell,
        .bigdash-settings-shell {
          height: 100%% !important;
          max-height: 100%%;
          min-height: %dpx;
        }
        .bigdash-app > .flex-grow-1 {
          min-height: 0;
          overflow: hidden;
        }
        html .content { min-height: 0; }
      ", mh, mh - 1L, mh)))
    ),
    shinyjs::useShinyjs(),
    bigLoaders::addBigLoaderDeps(),
    ...
  )
}

qsee_ui <- function(id, height = "100%", lazy = TRUE) {
  ns <- shiny::NS(id)

  ## This probe reports visibility of the whole Qsee module.  The individual
  ## board probes below only report their respective inner tabs. No plots
  ## are rendered directly under this module, so nothing is ever purged for
  ## it -- see the purge = FALSE passed to bd_is_visible() in qsee_server.R.
  ui <- shiny::tagList(
    shiny::tags$style(HTML("
      .bigdash-app,
      .big-full-page,
      .bigdash-sidebar-shell,
      .bigdash-settings-shell {
        height: 100% !important;
        max-height: 100%;
        min-height: 800px;
      }
      .bigdash-app > .flex-grow-1 {
        min-height: 0;
        overflow: hidden;
      }
      /* Beat Bootstrap d-md-block (display:block !important) so the
         settings pane is a column and developer options can sit at the bottom. */
      .bigdash-settings-shell.settings-expanded {
        display: flex !important;
        flex-direction: column !important;
        overflow: hidden;
      }
      .bigdash-settings-shell > .settings {
        flex: 1 1 auto;
        min-height: 0;
        display: flex;
        flex-direction: column;
        overflow: hidden;
      }
      .bigdash-settings-shell [id$='-settings-content'] {
        flex: 1 1 auto;
        min-height: 0;
        height: 100%;
        display: flex;
        flex-direction: column;
      }
      .bigdash-settings-shell .tab-settings:not(.d-none) {
        flex: 1 1 auto;
        min-height: 0;
        height: 100%;
        display: flex !important;
        flex-direction: column;
      }
      .qsee-dev-options {
        margin-top: auto;
      }
      html.bslib-has-full-screen .bslib-full-screen .html-widget.iheatmapr,
      html.bslib-has-full-screen .bslib-full-screen .iheatmapr.html-widget {
        height: calc(100vh - 90px) !important;
        width: 100% !important;
      }
    ")),
    bigdash::bd_visibility_probe(ns),
    bigdash::bigPage(
      id = id,
      #navbar = bigdash::navbar("Qsee/Bsee"),
      navbar = div(
        style = "visibility: hidden; display: none;"
      ),
      sidebar = bigdash::sidebar(
        "Qsee/Bsee",
        id = id,
        bigdash::sidebarItem("Normalization", ns("normalize-tab")),
        bigdash::sidebarItem("Missing values", ns("impute-tab")),
        bigdash::sidebarItem("PCA explorer", ns("pcaexplorer-tab")),
        bigdash::sidebarItem("Outlier analysis", ns("outlier-tab")),
        bigdash::sidebarItem("Batch-effects", ns("bsee-tab")),
        bigdash::sidebarItem("SD Filtering", ns("filtering-tab")),
        br(),
        div(
          style = "display: flex; flex-direction: column; gap: 6px;",
          shiny::actionButton(ns("upload"), "upload CSV", width = "100%"),
          shiny::actionButton(ns("load_example"), "load example", width = "100%"),
          shiny::actionButton(ns("load_pgx"), "use current pgx", width = "100%")
        )
      ),
      settings = bigdash::settings("Settings", id = id),
      ## bigdash::sidebarHelp(
      ##   id = id,
      ##   bigdash::sidebarTabHelp(
      ##     ns("normalize-tab"),
      ##     "Normalization",
      ##     "Compare different normalization methods on boxplots, densities and PCA. Includes noise injection for testing robustness."
      ##   ),
      ##   bigdash::sidebarTabHelp(
      ##     ns("impute-tab"),
      ##     "Missing values",
      ##     "Analyze imputation methods, missingness patterns (MAR/MNAR) and imputation accuracy using PCA, densities, heatmaps and validation plots."
      ##   ),
      ##   bigdash::sidebarTabHelp(
      ##     ns("outlier-tab"),
      ##     "Outlier analysis",
      ##     "Detect and visualize sample and feature outliers using z-scores, clustered heatmaps of extreme features, and PCA."
      ##   ),
      ##   bigdash::sidebarTabHelp(
      ##     ns("bsee-tab"),
      ##     "Batch-effects",
      ##     "Batch correction analysis. Clustering before/after correction, covariate correlation with phenotype vs PC, PVCA, and scores."
      ##   ),
      ##   bigdash::sidebarTabHelp(
      ##     ns("filtering-tab"),
      ##     "Filtering",
      ##     "SD filtering to remove low-variance features. Includes PCA of top-SD features, cumulative variance explained, and SD histogram."
      ##   )
      ## ),
      ## lazy = TRUE (default): each bigTabItem() holds only its sidebar
      ## inputs up front, and the board body (OmicsBoardUI + its visibility
      ## probe) is inserted lazily by bigdash::bigTabsLazy() in
      ## qsee_server.R, on that tab's first visit. lazy = FALSE builds every
      ## board's body up front instead, matching the pre-bigTabsLazy
      ## behaviour -- kept switchable for an A/B comparison of the two.
      ## Whatever is passed here must match the `lazy` given to
      ## qsee_server() for the same `id`.
      bigdash::bigTabs(
        id = id,
        bigdash::bigTabItem(
          ns("normalize-tab"),
          qsee_normalization_inputs(ns("normalize")),
          if (!lazy) qsee_normalization_ui(ns("normalize"))
        ),
        bigdash::bigTabItem(
          ns("impute-tab"),
          qsee_imputation_inputs(ns("impute")),
          if (!lazy) qsee_imputation_ui(ns("impute"))
        ),
        bigdash::bigTabItem(
          ns("pcaexplorer-tab"),
          qsee_pcaexplorer_inputs(ns("pcaexplorer")),
          if (!lazy) qsee_pcaexplorer_ui(ns("pcaexplorer"))
        ),
        bigdash::bigTabItem(
          ns("outlier-tab"),
          qsee_outlier_inputs(ns("outlier")),
          if (!lazy) qsee_outlier_ui(ns("outlier"))
        ),
        bigdash::bigTabItem(
          ns("bsee-tab"),
          qsee_bsee_inputs(ns("bsee")),
          if (!lazy) qsee_bsee_ui(ns("bsee"))
        ),
        bigdash::bigTabItem(
          ns("filtering-tab"),
          qsee_filtering_inputs(ns("filtering")),
          if (!lazy) qsee_filtering_ui(ns("filtering"))
        )
      )
    )
  )

  return(ui)
}
