## make as R6 class?? e.g. add documentation, initialize object,
## object id.

## Registry for the Expression menu group: which tabs it contains, their
## sidebar labels (module_menu), the shell each tab renders at boot
## (module_ui) and the UI + server bigTabsLazy() fills it with on first
## visit (module_lazy). Read by opg_ui() and opg_server().
##
## The boards themselves live in components/board.*; a group draws on
## whichever of them it needs.

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
  module_ui = function(ns = identity) {
    list(
      bigdash::bigTabItem(
        ns("diffexpr-tab"),
        ExpressionInputs(ns("diffexpr")),
        create_loader(ns("diffexpr-loader"))
      ),
      bigdash::bigTabItem(
        ns("timeseries-tab"),
        TimeSeriesInputs(ns("timeseries")),
        create_loader(ns("timeseries-loader"))
      ),
      bigdash::bigTabItem(
        ns("corr-tab"),
        CorrelationInputs(ns("corr")),
        create_loader(ns("corr-loader"))
      ),
      bigdash::bigTabItem(
        ns("bio-tab"),
        BiomarkerInputs(ns("bio")),
        create_loader(ns("bio-loader"))
      )
    )
  },
  module_lazy = function(PGX, labeltype = NULL, ns = identity, ...) {
    list(
      "corr-tab" = list(
        ui     = function() CorrelationUI(ns("corr")),
        server = function() CorrelationBoard("corr", pgx = PGX, labeltype = labeltype)
      ),
      "timeseries-tab" = list(
        ui     = function() TimeSeriesUI(ns("timeseries")),
        server = function() TimeSeriesBoard("timeseries", pgx = PGX, labeltype = labeltype)
      ),
      "bio-tab" = list(
        ui     = function() BiomarkerUI(ns("bio")),
        server = function() BiomarkerBoard("bio", pgx = PGX)
      )
    )
  }
)
