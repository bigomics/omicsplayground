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
  ## Per-board servers, keyed by tab, so opg_server can start just the board
  ## that was opened. module_server() above stays for the eager path.
  module_servers = list(
    "clustersamples-tab" = function(PGX, labeltype = NULL)
      ClusteringBoard("clustersamples", pgx = PGX, labeltype = labeltype),
    "clusterfeatures-tab" = function(PGX, labeltype = NULL)
      FeatureMapBoard("clusterfeatures", pgx = PGX, labeltype = labeltype)
  ),
  module_server = function(PGX, labeltype = NULL) {
    info("[SERVER] calling ClusteringBoard module")
    ClusteringBoard("clustersamples",
      pgx = PGX, labeltype = labeltype
    )

    info("[SERVER] calling FeatureMapBoard module")
    FeatureMapBoard("clusterfeatures",
      pgx = PGX, labeltype = labeltype
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
