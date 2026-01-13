## make as R6 class?? e.g. add documentation, initialize object,
## object id.

MODULE.expression <- list(
  ## MODULES[["multiomics"]] <- list(
  module_menu = function() {
    c(
      diffexpr = "Differential expression",
      timeseries = "TimeSeries",
      corr = "Correlation analysis",
      bio = "Find biomarkers"
    )
  },
  module_server = function(PGX, labeltype = NULL) {
    info("[SERVER] calling CorrelationBoard module")
    CorrelationBoard("corr",
      pgx = PGX, labeltype = labeltype
    )

    info("[SERVER] calling BiomarkerBoard module")
    BiomarkerBoard(
      "bio",
      pgx = PGX
    )

    info("[SERVER] calling TimeSeries module")
    TimeSeriesBoard("timeseries",
      pgx = PGX, labeltype = labeltype
    )
  },
  module_ui = function() {
    list(
      bigdash::bigTabItem(
        "diffexpr-tab",
        ExpressionInputs("diffexpr"),
        create_loader("diffexpr-loader")
      ),
      bigdash::bigTabItem(
        "timeseries-tab",
        TimeSeriesInputs("timeseries"),
        create_loader("timeseries-loader")
      ),
      bigdash::bigTabItem(
        "corr-tab",
        CorrelationInputs("corr"),
        create_loader("corr-loader")
      ),
      bigdash::bigTabItem(
        "bio-tab",
        BiomarkerInputs("bio"),
        create_loader("bio-loader")
      )
    )
  },
  module_ui2 = function() {
    list(
      list(
        "corr-tab",
        CorrelationUI("corr")
      ),
      list(
        "timeseries-tab",
        TimeSeriesUI("timeseries")
      ),
      list(
        "diffexpr-tab",
        ExpressionUI("diffexpr")
      ),
      list(
        "bio-tab",
        BiomarkerUI("bio")
      )
    )
  }
)
