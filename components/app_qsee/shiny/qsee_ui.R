##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


qsee_ui <- function(id, height = "100%") {
  ns <- shiny::NS(id)

  ## This probe reports visibility of the whole Qsee module.  The individual
  ## board probes below only report their respective inner tabs.
  ui <- shiny::tagList(
    qsee_visibility_probe(ns),
    qsee_plot_spinner_js(ns),
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
        shiny::actionButton(ns("upload"), "upload CSV"),
        shiny::actionButton(ns("load_example"), "load example"),
        shiny::actionButton(ns("load_pgx"), "use current dataset")
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
      bigdash::bigTabs(
        id = id,
        bigdash::bigTabItem(
          ns("normalize-tab"),
          qsee_normalization_inputs(ns("normalize")),
          qsee_normalization_ui(ns("normalize"))
        ),
        bigdash::bigTabItem(
          ns("impute-tab"),
          qsee_imputation_inputs(ns("impute")),
          qsee_imputation_ui(ns("impute"))
        ),
        bigdash::bigTabItem(
          ns("pcaexplorer-tab"),
          qsee_pcaexplorer_inputs(ns("pcaexplorer")),
          qsee_pcaexplorer_ui(ns("pcaexplorer"))
        ),
        bigdash::bigTabItem(
          ns("outlier-tab"),
          qsee_outlier_inputs(ns("outlier")),
          qsee_outlier_ui(ns("outlier"))
        ),
        bigdash::bigTabItem(
          ns("bsee-tab"),
          qsee_bsee_inputs(ns("bsee")),
          qsee_bsee_ui(ns("bsee"))
        ),
        bigdash::bigTabItem(
          ns("filtering-tab"),
          qsee_filtering_inputs(ns("filtering")),
          qsee_filtering_ui(ns("filtering"))
        )
      )
    )
  )

  return(ui)
}
