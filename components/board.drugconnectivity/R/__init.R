#
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

MODULE.systems <- list(
  module_menu = function() {
    c(
      drug = "Drug connectivity",
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
  ## Per-board servers, keyed by tab, so opg_server can start just the board
  ## that was opened. module_server() above stays for the eager path.
  module_servers = list(
    "drug-tab" = function(PGX)
      DrugConnectivityBoard("drug", pgx = PGX),
    "tcga-tab" = function(PGX)
      TcgaBoard("tcga", pgx = PGX),
    "cell-tab" = function(PGX)
      SingleCellBoard("cell", pgx = PGX),
    "pcsf-tab" = function(PGX)
      PcsfBoard("pcsf", pgx = PGX)
  ),
  module_server = function(PGX) {
    DrugConnectivityBoard("drug",
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
  }
)
