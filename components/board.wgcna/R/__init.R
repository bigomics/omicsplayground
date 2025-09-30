#
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

MODULE.systems <- list(
  module_menu = function() {
    c(
      drug = "Drug connectivity",
      wgcna = "Network analysis",
      tcga = "TCGA analysis",
      cell = "Single cell",
      pcsf = "Protein networks",
      consensus = "Consensus WGCNA",
      preservation = "Preservation WGCNA"
    )
  },
  module_ui = function() {
    list(
      bigdash::bigTabItem(
        "drug-tab",
        DrugConnectivityInputs("drug"),
        create_loader("drug-loader")
      ),
      bigdash::bigTabItem(
        "wgcna-tab",
        WgcnaInputs("wgcna"),
        create_loader("wgcna-loader")
      ),
      bigdash::bigTabItem(
        "tcga-tab",
        TcgaInputs("tcga"),
        create_loader("tcga-loader")
      ),
      bigdash::bigTabItem(
        "cell-tab",
        SingleCellInputs("cell"),
        create_loader("cell-loader")
      ),
      bigdash::bigTabItem(
        "pcsf-tab",
        PcsfInputs("pcsf"),
        create_loader("pcsf-loader")
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
      )
    )
  },
  module_ui2 = function() {
    list(
      list(
        "drug-tab",
        DrugConnectivityUI("drug")
      ),
      list(
        "wgcna-tab",
        WgcnaUI("wgcna")
      ),
      list(
        "tcga-tab",
        TcgaUI("tcga")
      ),
      list(
        "cell-tab",
        SingleCellUI("cell")
      ),
      list(
        "pcsf-tab",
        PcsfUI("pcsf")
      ),
      list(
        "consensus-tab",
        ConsensusWGCNA_UI("consensus")
      ),
      list(
        "preservation-tab",
        PreservationWGCNA_UI("preservation")
      )
    )
  },
  module_server = function(PGX) {
    DrugConnectivityBoard("drug",
      pgx = PGX
    )

    WgcnaBoard("wgcna",
      pgx = PGX
    )

    TcgaBoard("tcga",
      pgx = PGX
    )

    SingleCellBoard("cell",
      pgx = PGX
    )

    PcsfBoard("pcsf",
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
    
  }
)
