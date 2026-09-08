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
  module_ui = function(ns = identity) {
    list(
      bigdash::bigTabItem(
        ns("drug-tab"),
        DrugConnectivityInputs(ns("drug")),
        create_loader(ns("drug-loader"))
      ),
      bigdash::bigTabItem(
        ns("tcga-tab"),
        TcgaInputs(ns("tcga")),
        create_loader(ns("tcga-loader"))
      ),
      bigdash::bigTabItem(
        ns("cell-tab"),
        SingleCellInputs(ns("cell")),
        create_loader(ns("cell-loader"))
      ),
      bigdash::bigTabItem(
        ns("pcsf-tab"),
        PcsfInputs(ns("pcsf")),
        create_loader(ns("pcsf-loader"))
      )
    )
  },
  module_lazy = function(PGX, ns = identity, ...) {
    list(
      "drug-tab" = list(
        ui     = function() DrugConnectivityUI(ns("drug")),
        server = function() DrugConnectivityBoard("drug", pgx = PGX)
      ),
      "tcga-tab" = list(
        ui     = function() TcgaUI(ns("tcga")),
        server = function() TcgaBoard("tcga", pgx = PGX)
      ),
      "cell-tab" = list(
        ui     = function() SingleCellUI(ns("cell")),
        server = function() SingleCellBoard("cell", pgx = PGX)
      ),
      "pcsf-tab" = list(
        ui     = function() PcsfUI(ns("pcsf")),
        server = function() PcsfBoard("pcsf", pgx = PGX)
      )
    )
  }
)
