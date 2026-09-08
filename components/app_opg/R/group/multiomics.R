## make as R6 class?? e.g. add documentation, initialize object,
## object id.

## Registry for the MultiOmics menu group: which tabs it contains, their
## sidebar labels (module_menu), the shell each tab renders at boot
## (module_ui) and the UI + server bigTabsLazy() fills it with on first
## visit (module_lazy). Read by opg_ui() and opg_server().
##
## The boards themselves live in components/board.*; a group draws on
## whichever of them it needs.

MODULE.multiomics <- list(
  module_menu = function() {
    c(
      snf = "SNF",
      lasagna = "LASAGNA",
      mgsea = "Multiomics GSEA",
      mofa = "MOFA",
      deepnet = "DeepLearning"
    )
  },
  module_ui = function(ns = identity) {
    list(
      bigdash::bigTabItem(ns("mofa-tab"), MofaInputs(ns("mofa")), create_loader(ns("mofa-loader"))),
      bigdash::bigTabItem(ns("mgsea-tab"), MGseaInputs(ns("mgsea")), create_loader(ns("mgsea-loader"))),
      bigdash::bigTabItem(ns("snf-tab"), SNFInputs(ns("snf")), create_loader(ns("snf-loader"))),
      bigdash::bigTabItem(ns("lasagna-tab"), LasagnaInputs(ns("lasagna")), create_loader(ns("lasagna-loader"))),
      bigdash::bigTabItem(ns("deepnet-tab"), DeepNetInputs(ns("deepnet")), create_loader(ns("deepnet-loader")))
    )
  },
  module_lazy = function(PGX, ns = identity, ...) {
    list(
      "mofa-tab" = list(
        ui     = function() MofaUI(ns("mofa")),
        server = function() MofaBoard("mofa", pgx = PGX)
      ),
      "mgsea-tab" = list(
        ui     = function() MGseaUI(ns("mgsea")),
        server = function() MGseaBoard("mgsea", pgx = PGX)
      ),
      "snf-tab" = list(
        ui     = function() SNFUI(ns("snf")),
        server = function() SNFBoard("snf", pgx = PGX)
      ),
      "lasagna-tab" = list(
        ui     = function() LasagnaUI(ns("lasagna")),
        server = function() LasagnaBoard("lasagna", pgx = PGX)
      ),
      "deepnet-tab" = list(
        ui     = function() DeepNetUI(ns("deepnet")),
        server = function() DeepNetBoard("deepnet", pgx = PGX)
      )
    )
  },
  module_help = function(ns = identity) {
    list(
      bigdash::sidebarTabHelp(
        ns("mofa-tab"), "MOFA",
        tspan("Multi-omics Factor Analysis (MOFA) is a multi-omics
                  integration method based on matrix factorization.")
      ),
      bigdash::sidebarTabHelp(
        ns("mgsea-tab"), "multiGSEA",
        tspan("multiGSEA performs multi-omics integration on gene set level.")
      ),
      bigdash::sidebarTabHelp(
        ns("snf-tab"), "SNF",
        tspan("SNF clustering")
      ),
      bigdash::sidebarTabHelp(
        ns("lasagna-tab"), "Lasagna",
        tspan("LASAGNA is a stacked layer model for multi-omics integration where each layer corresponds to a datatype.")
      ),
      bigdash::sidebarTabHelp(
        ns("deepnet-tab"), "DeepLearning",
        tspan("Integration using DeepLearning")
      )
    )
  }
)
