##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


qsee_ui <- function(id, height = "100%", lazy = TRUE) {
  ns <- shiny::NS(id)

  ## This probe reports visibility of the whole Qsee module.  The individual
  ## board probes below only report their respective inner tabs. No plots
  ## are rendered directly under this module, so nothing is ever purged for
  ## it -- see the purge = FALSE passed to bd_is_visible() in qsee_server.R.
  ui <- shiny::tagList(
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
