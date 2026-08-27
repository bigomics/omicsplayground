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
  module_ui = function() {
    list(
      bigdash::bigTabItem("mofa-tab", MofaInputs("mofa"), create_loader("mofa-loader")),
      bigdash::bigTabItem("mgsea-tab", MGseaInputs("mgsea"), create_loader("mgsea-loader")),
      bigdash::bigTabItem("snf-tab", SNFInputs("snf"), create_loader("snf-loader")),
      bigdash::bigTabItem("lasagna-tab", LasagnaInputs("lasagna"), create_loader("lasagna-loader")),
      bigdash::bigTabItem("deepnet-tab", DeepNetInputs("deepnet"), create_loader("deepnet-loader"))
    )
  },
  module_lazy = function(PGX, ...) {
    list(
      "mofa-tab" = list(
        ui     = function() MofaUI("mofa"),
        server = function() MofaBoard("mofa", pgx = PGX)
      ),
      "mgsea-tab" = list(
        ui     = function() MGseaUI("mgsea"),
        server = function() MGseaBoard("mgsea", pgx = PGX)
      ),
      "snf-tab" = list(
        ui     = function() SNFUI("snf"),
        server = function() SNFBoard("snf", pgx = PGX)
      ),
      "lasagna-tab" = list(
        ui     = function() LasagnaUI("lasagna"),
        server = function() LasagnaBoard("lasagna", pgx = PGX)
      ),
      "deepnet-tab" = list(
        ui     = function() DeepNetUI("deepnet"),
        server = function() DeepNetBoard("deepnet", pgx = PGX)
      )
    )
  },
  module_help = function() {
    list(
      bigdash::sidebarTabHelp(
        "mofa-tab", "MOFA",
        tspan("Multi-omics Factor Analysis (MOFA) is a multi-omics
                  integration method based on matrix factorization.")
      ),
      bigdash::sidebarTabHelp(
        "mgsea-tab", "multiGSEA",
        tspan("multiGSEA performs multi-omics integration on gene set level.")
      ),
      bigdash::sidebarTabHelp(
        "snf-tab", "SNF",
        tspan("SNF clustering")
      ),
      bigdash::sidebarTabHelp(
        "lasagna-tab", "Lasagna",
        tspan("LASAGNA is a stacked layer model for multi-omics integration where each layer corresponds to a datatype.")
      ),
      bigdash::sidebarTabHelp(
        "deepnet-tab", "DeepLearning",
        tspan("Integration using DeepLearning")
      )
    )
  }
)
