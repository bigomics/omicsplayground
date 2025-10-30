## make as R6 class?? e.g. add documentation, initialize object,
## object id.

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
  module_ui2 = function() {
    list(
      list("mofa-tab", MofaUI("mofa")),
      list("mgsea-tab", MGseaUI("mgsea")),
      list("snf-tab", SNFUI("snf")),
      list("lasagna-tab", LasagnaUI("lasagna")),
      list("deepnet-tab", DeepNetUI("deepnet"))
    )
  },
  module_server = function(PGX) {
    info("[SERVER] calling MofaBoard module")
    MofaBoard("mofa", pgx = PGX)

    info("[SERVER] calling MGseaBoard module")
    MGseaBoard("mgsea", pgx = PGX)

    info("[SERVER] calling SNFBoard module")
    SNFBoard("snf", pgx = PGX)

    info("[SERVER] calling LasagnaBoard module")
    LasagnaBoard("lasagna", pgx = PGX)

    info("[SERVER] calling DeepNetBoard module")
    DeepNetBoard("deepnet", pgx = PGX)

  },
  module_help = function() {
    list(
      bigdash::sidebarTabHelp(
        "mofa-tab", "MOFA",
        tspan("Multi-omics Factor Analysis (MOFA) is a multi-omics
                  integration method based on matrix factorization.")),

      bigdash::sidebarTabHelp("mgsea-tab", "multiGSEA",
        tspan("multiGSEA performs multi-omics integration on gene set level.")),

      bigdash::sidebarTabHelp("snf-tab", "SNF",
        tspan("SNF clustering")),

      bigdash::sidebarTabHelp("lasagna-tab", "Lasagna",
        tspan("LASAGNA is a stacked layer model for multi-omics integration where each layer corresponds to a datatype.")),
      
      bigdash::sidebarTabHelp("deepnet-tab", "DeepLearning",
        tspan("Integration using DeepLearning"))
            
    )
  }
)
