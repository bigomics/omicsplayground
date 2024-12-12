## make as R6 class?? e.g. add documentation, initialize object,
## object id.

MODULE.multiomics <- list( 
##MODULES[["multiomics"]] <- list(

  module_menu = function() {
    c(
      mofa = "MOFA",
      mgsea = "multiGSEA",
      snf = "SNF",
      lasagna = "Lasagna"          
    )
  },
  
  module_server = function(PGX) {
    info("[SERVER] calling MofaBoard module")
    MofaBoard("mofa", pgx = PGX)

    info("[SERVER] calling MGseaBoard module")
    MGseaBoard("mgsea", pgx = PGX)
          
    info("[SERVER] calling SNF_Board module")
    SNF_Board("snf", pgx = PGX)

    info("[SERVER] calling LasagnaBoard module")
    LasagnaBoard("lasagna", pgx = PGX)
  },

  module_ui = function() {
    list(
      bigdash::bigTabItem(   ## call this in app.R??
        "mofa-tab",
        MofaInputs("mofa"),
        MofaUI("mofa")
      ),        
      bigdash::bigTabItem(
        "mgsea-tab",
        MGseaInputs("mgsea"),
        MGseaUI("mgsea")
      ),        
      bigdash::bigTabItem(
        "snf-tab",
        SNF_Inputs("snf"),
        SNF_UI("snf")
      ),        
      bigdash::bigTabItem(
        "lasagna-tab",
        LasagnaInputs("lasagna"),
        LasagnaUI("lasagna")
      )        
    )
  },    

  module_help = function() {
    list(
      bigdash::sidebarTabHelp(  ## call this function in app.R?
        "mofa-tab",
        "MOFA",
        tspan("Multi-omics Factor Analysis (MOFA) is a multi-omics
                  integration method based on multi-omcis factor analysis.")
      ),
      bigdash::sidebarTabHelp(
        "mgsea-tab",
        "multiGSEA",
        tspan("multiGSEA perform multi-omics integration on gene set level.")
      ),
      bigdash::sidebarTabHelp(
        "snf-tab",
        "SNF",
        tspan("SNF clustering")
      )
      ## bigdash::sidebarTabHelp(
      ##   "mpcsf-tab",
      ##   "multiPCSF",
      ##   tspan("Perfom multiPCSF analysis")
      ## )
    )
  }
  
  
)
