#
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## Registry for the SystemsBio menu group: which tabs it contains, their
## sidebar labels (module_menu), the shell each tab renders at boot
## (module_ui) and the UI + server bigTabsLazy() fills it with on first
## visit (module_lazy). Read by opg_ui() and opg_server().
##
## The boards themselves live in components/board.*; a group draws on
## whichever of them it needs.

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
  module_lazy = function(PGX, ...) {
    list(
      "drug-tab" = list(
        ui     = function() DrugConnectivityUI("drug"),
        server = function() DrugConnectivityBoard("drug", pgx = PGX)
      ),
      "tcga-tab" = list(
        ui     = function() TcgaUI("tcga"),
        server = function() TcgaBoard("tcga", pgx = PGX)
      ),
      "cell-tab" = list(
        ui     = function() SingleCellUI("cell"),
        server = function() SingleCellBoard("cell", pgx = PGX)
      ),
      "pcsf-tab" = list(
        ui     = function() PcsfUI("pcsf"),
        server = function() PcsfBoard("pcsf", pgx = PGX)
      )
    )
  }
)
