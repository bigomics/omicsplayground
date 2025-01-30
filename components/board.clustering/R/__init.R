## make as R6 class?? e.g. add documentation, initialize object,
## object id.

MODULE.clustering <- list(
  ## MODULES[["multiomics"]] <- list(
  module_menu = function() {
    c(
      clustersamples = "Samples",
      clusterfeatures = "Features"
    )
  },
  module_server = function(PGX, board_observers = NULL, labeltype = NULL) {
    info("[SERVER] calling ClusteringBoard module")
    ClusteringBoard("clustersamples",
      pgx = PGX, labeltype = labeltype,
      board_observers = board_observers
    )

    info("[SERVER] calling FeatureMapBoard module")
    FeatureMapBoard("clusterfeatures",
      pgx = PGX, labeltype = labeltype,
      board_observers = board_observers
    )
  },
  module_ui = function() {
    list(
      bigdash::bigTabItem(
        "clustersamples-tab",
        ClusteringInputs("clustersamples")
      ),
      bigdash::bigTabItem(
        "clusterfeatures-tab",
        FeatureMapInputs("clusterfeatures")
      )
    )
  },
  module_ui2 = function() {
    list(
      list(
        "clustersamples-tab",
        ClusteringUI("clustersamples")
      ),
      list(
        "clusterfeatures-tab",
        FeatureMapUI("clusterfeatures")
      )
    )
  }
)
