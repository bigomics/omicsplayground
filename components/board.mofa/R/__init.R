## make as R6 class?? e.g. add documentation, initialize object,
## object id.

MODULE.multiomics <- list(
  module_menu = function() {
    c(
      snf = "SNF",
      lasagna = "LASAGNA",
      mgsea = "multiGSEA",
      mofa = "MOFA",
      deepnet = "DeepLearning"
    )
  },

  module_ui = function() {
    list(
      bigdash::bigTabItem(   ## call this in app.R??
        "mofa-tab",
        MofaInputs("mofa"),
        create_loader("mofa-loader")
      ),
      bigdash::bigTabItem(
        "mgsea-tab",
        MGseaInputs("mgsea"),
        create_loader("mgsea-loader")
      ),
      bigdash::bigTabItem(
        "snf-tab",
        SNF_Inputs("snf"),
        create_loader("snf-loader")
      ),
      bigdash::bigTabItem(
        "lasagna-tab",
        LasagnaInputs("lasagna"),
        create_loader("lasagna-loader")
      ),
      bigdash::bigTabItem(
        "deepnet-tab",
        DeepNetInputs("deepnet"),
        create_loader("deepnet-loader")
      )
    )
  },
  module_ui2 = function() {
    list(
      list(   ## call this in app.R??
        "mofa-tab",
        MofaUI("mofa")
      ),
      list(
        "mgsea-tab",
        MGseaUI("mgsea")
      ),
      list(
        "snf-tab",
        SNF_UI("snf")
      ),
      list(
        "lasagna-tab",
        LasagnaUI("lasagna")
      ),
      list(
        "deepnet-tab",
        DeepNetUI("deepnet")
      )
    )
  },
  module_server = function(PGX, board_observers = NULL) {
    info("[SERVER] calling MofaBoard module")
    MofaBoard("mofa", pgx = PGX,
      board_observers = board_observers
    )

    info("[SERVER] calling MGseaBoard module")
    MGseaBoard("mgsea", pgx = PGX,
      board_observers = board_observers
    )
          
    info("[SERVER] calling SNF_Board module")
    SNF_Board("snf", pgx = PGX,
      board_observers = board_observers
    )

    info("[SERVER] calling LasagnaBoard module")
    LasagnaBoard("lasagna", pgx = PGX,
      board_observers = board_observers
    )

    info("[SERVER] calling DeepNetBoard module")
    DeepNetBoard("deepnet", pgx = PGX,
      board_observers = board_observers
    )
  },

  module_help = function() {
    list(
      bigdash::sidebarTabHelp(  ## call this function in app.R?
        "mofa-tab",
        "MOFA",
        tspan("Multi-omics Factor Analysis (MOFA) is a multi-omics
                  integration method based on matrix factorization.")
      ),
      bigdash::sidebarTabHelp(
        "mgsea-tab",
        "multiGSEA",
        tspan("multiGSEA performs multi-omics integration on gene set level.")
      ),
      bigdash::sidebarTabHelp(
        "snf-tab",
        "SNF",
        tspan("SNF clustering")
      ),
      bigdash::sidebarTabHelp(
        "lasagna-tab",
        "Lasagna",
        tspan("LASAGNA is a stacked layer model for multi-omics integration where each layer corresponds to a datatype.")
      ),
      bigdash::sidebarTabHelp(
        "deepnet-tab",
        "DeepLearning",
        tspan("Integration using DeepLearning")
      )
      ## bigdash::sidebarTabHelp(
      ##   "mpcsf-tab",
      ##   "multiPCSF",
      ##   tspan("Perfom multiPCSF analysis")
      ## )
    )
  }
  
  
)
