## make as R6 class?? e.g. add documentation, initialize object,
## object id.

MODULE.clustering <- list(
  ## MODULES[["multiomics"]] <- list(
  module_menu = function() {
    c(
      clustersamples = "Samples",
      clusterfeatures = "Features",
      timeseries = "TimeSeries"
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

    info("[SERVER] calling TimeSeries module")
    TimeSeriesBoard("timeseries",
      pgx = PGX, labeltype = labeltype,
      board_observers = board_observers
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
      ),
      bigdash::bigTabItem(
        "timeseries-tab",
        TimeSeriesInputs("timeseries"),
        create_loader("timeseries-loader")
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
      ),
      list(
        "timeseries-tab",
        TimeSeriesUI("timeseries")
      )
    )
  }
)
