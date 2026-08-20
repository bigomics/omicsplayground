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
  module_lazy = function(PGX, labeltype = NULL, ...) {
    list(
      "corr-tab" = list(
        ui     = function() CorrelationUI("corr"),
        server = function() CorrelationBoard("corr", pgx = PGX, labeltype = labeltype)
      ),
      "timeseries-tab" = list(
        ui     = function() TimeSeriesUI("timeseries"),
        server = function() TimeSeriesBoard("timeseries", pgx = PGX, labeltype = labeltype)
      ),
      "bio-tab" = list(
        ui     = function() BiomarkerUI("bio"),
        server = function() BiomarkerBoard("bio", pgx = PGX)
      )
    )
  }
)
