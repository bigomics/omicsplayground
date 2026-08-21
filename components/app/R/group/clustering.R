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
  module_ui = function() {
    list(
      bigdash::bigTabItem(
        "clustersamples-tab",
        ClusteringInputs("clustersamples"),
        create_loader("clustersamples-loader")
      ),
      bigdash::bigTabItem(
        "clusterfeatures-tab",
        FeatureMapInputs("clusterfeatures"),
        create_loader("clusterfeatures-loader")
      )
    )
  },
  module_lazy = function(PGX, labeltype = NULL, ...) {
    list(
      "clustersamples-tab" = list(
        ui     = function() ClusteringUI("clustersamples"),
        server = function() ClusteringBoard("clustersamples", pgx = PGX, labeltype = labeltype)
      ),
      "clusterfeatures-tab" = list(
        ui     = function() FeatureMapUI("clusterfeatures"),
        server = function() FeatureMapBoard("clusterfeatures", pgx = PGX, labeltype = labeltype)
      )
    )
  }
)
