## make as R6 class?? e.g. add documentation, initialize object,
## object id.

## Registry for the Clustering menu group: which tabs it contains, their
## sidebar labels (module_menu), the shell each tab renders at boot
## (module_ui) and the UI + server bigTabsLazy() fills it with on first
## visit (module_lazy). Read by opg_ui() and opg_server().
##
## The boards themselves live in components/board.*; a group draws on
## whichever of them it needs.

MODULE.clustering <- list(
  ## MODULES[["multiomics"]] <- list(
  module_menu = function() {
    c(
      clustersamples = "Samples",
      clusterfeatures = "Features"
    )
  },
  module_ui = function(ns = identity) {
    list(
      bigdash::bigTabItem(
        ns("clustersamples-tab"),
        ClusteringInputs(ns("clustersamples")),
        create_loader(ns("clustersamples-loader"))
      ),
      bigdash::bigTabItem(
        ns("clusterfeatures-tab"),
        FeatureMapInputs(ns("clusterfeatures")),
        create_loader(ns("clusterfeatures-loader"))
      )
    )
  },
  module_lazy = function(PGX, labeltype = NULL, ns = identity, ...) {
    list(
      "clustersamples-tab" = list(
        ui     = function() ClusteringUI(ns("clustersamples")),
        server = function() ClusteringBoard("clustersamples", pgx = PGX, labeltype = labeltype)
      ),
      "clusterfeatures-tab" = list(
        ui     = function() FeatureMapUI(ns("clusterfeatures")),
        server = function() FeatureMapBoard("clusterfeatures", pgx = PGX, labeltype = labeltype)
      )
    )
  }
)
