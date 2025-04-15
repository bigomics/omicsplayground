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
      pcsf = "Protein networks"
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
      )
    )
  },
  module_server = function(PGX, board_observers = NULL) {
    DrugConnectivityBoard("drug",
      pgx = PGX,
      board_observers = board_observers
    )

    WgcnaBoard("wgcna",
      pgx = PGX#,
      #board_observers = board_observers
    )

    TcgaBoard("tcga",
      pgx = PGX,
      board_observers = board_observers
    )

    SingleCellBoard("cell",
      pgx = PGX,
      board_observers = board_observers
    )

    PcsfBoard("pcsf",
      pgx = PGX,
      board_observers = board_observers
    )
  }
)

