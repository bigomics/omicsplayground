#
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

MODULE.wgcna <- list(
  module_menu = function() {
    c(
      wgcna = HTML("Standard WGCNA"),
      consensus = "Consensus WGCNA",
      preservation = "Preservation WGCNA",
      mwgcna = "Multiomics WGCNA"
    )
  },
  
  module_ui = function() {
    list(
      bigdash::bigTabItem(
        "wgcna-tab",
        WgcnaInputs("wgcna"),
        create_loader("wgcna-loader")
      ),
      bigdash::bigTabItem(
        "consensus-tab",
        ConsensusWGCNA_Inputs("consensus"),
        create_loader("consensus-loader")
      ),
      bigdash::bigTabItem(
        "preservation-tab",
        PreservationWGCNA_Inputs("preservation"),
        create_loader("preservation-loader")
      ),
      bigdash::bigTabItem(
        "mwgcna-tab",
        MultiWGCNA_Inputs("mwgcna"),
        create_loader("mwgcna-loader")
      )
    )
  },
  
  module_ui2 = function() {
    list(
      list(
        "wgcna-tab",
        WgcnaUI("wgcna")
      ),
      list(
        "consensus-tab",
        ConsensusWGCNA_UI("consensus")
      ),
      list(
        "preservation-tab",
        PreservationWGCNA_UI("preservation")
      ),
      list(
        "mwgcna-tab",
        MultiWGCNA_UI("mwgcna")
      )
    )
  },
  
  module_server = function(PGX) {

    WgcnaBoard("wgcna",
      pgx = PGX
    )

    ConsensusWGCNA_Board(
      id = "consensus",
      pgx = PGX
    )     

    PreservationWGCNA_Board(
      id = "preservation",
      pgx = PGX
    )     

    MultiWGCNA_Board(
      "mwgcna",
      pgx = PGX
    )
  },

  module_help = function() {
    list(
      bigdash::sidebarTabHelp("mwgcna-tab", "MultiOmics WGCNA",
        tspan("WGCNA for multi-omics"))
    )
  }
  
)
