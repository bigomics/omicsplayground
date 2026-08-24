#
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## Registry for the WGCNA menu group: which tabs it contains, their
## sidebar labels (module_menu), the shell each tab renders at boot
## (module_ui) and the UI + server bigTabsLazy() fills it with on first
## visit (module_lazy). Read by opg_ui() and opg_server().
##
## The boards themselves live in components/board.*; a group draws on
## whichever of them it needs.

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
  module_lazy = function(PGX, save_pgx = NULL, ...) {
    list(
      "wgcna-tab" = list(
        ui     = function() WgcnaUI("wgcna"),
        server = function() WgcnaBoard("wgcna", pgx = PGX, save_pgx = save_pgx)
      ),
      "consensus-tab" = list(
        ui     = function() ConsensusWGCNA_UI("consensus"),
        server = function() ConsensusWGCNA_Board(id = "consensus", pgx = PGX)
      ),
      "preservation-tab" = list(
        ui     = function() PreservationWGCNA_UI("preservation"),
        server = function() PreservationWGCNA_Board(id = "preservation", pgx = PGX)
      ),
      "mwgcna-tab" = list(
        ui     = function() MultiWGCNA_UI("mwgcna"),
        server = function() MultiWGCNA_Board("mwgcna", pgx = PGX, save_pgx = save_pgx)
      )
    )
  },
  module_help = function() {
    list(
      bigdash::sidebarTabHelp(
        "mwgcna-tab", "MultiOmics WGCNA",
        tspan("WGCNA for multi-omics")
      )
    )
  }
)
