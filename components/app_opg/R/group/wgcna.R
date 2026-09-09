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
  module_ui = function(ns = identity) {
    list(
      bigdash::bigTabItem(
        ns("wgcna-tab"),
        WgcnaInputs(ns("wgcna")),
        create_loader(ns("wgcna-loader"))
      ),
      bigdash::bigTabItem(
        ns("consensus-tab"),
        ConsensusWGCNA_Inputs(ns("consensus")),
        create_loader(ns("consensus-loader"))
      ),
      bigdash::bigTabItem(
        ns("preservation-tab"),
        PreservationWGCNA_Inputs(ns("preservation")),
        create_loader(ns("preservation-loader"))
      ),
      bigdash::bigTabItem(
        ns("mwgcna-tab"),
        MultiWGCNA_Inputs(ns("mwgcna")),
        create_loader(ns("mwgcna-loader"))
      )
    )
  },
  module_lazy = function(PGX, save_pgx = NULL, ns = identity, ...) {
    list(
      "wgcna-tab" = list(
        ui     = function() WgcnaUI(ns("wgcna")),
        server = function() WgcnaBoard("wgcna", pgx = PGX, save_pgx = save_pgx)
      ),
      "consensus-tab" = list(
        ui     = function() ConsensusWGCNA_UI(ns("consensus")),
        server = function() ConsensusWGCNA_Board(id = "consensus", pgx = PGX)
      ),
      "preservation-tab" = list(
        ui     = function() PreservationWGCNA_UI(ns("preservation")),
        server = function() PreservationWGCNA_Board(id = "preservation", pgx = PGX)
      ),
      "mwgcna-tab" = list(
        ui     = function() MultiWGCNA_UI(ns("mwgcna")),
        server = function() MultiWGCNA_Board("mwgcna", pgx = PGX, save_pgx = save_pgx)
      )
    )
  },
  module_help = function(ns = identity) {
    list(
      bigdash::sidebarTabHelp(
        ns("mwgcna-tab"), "MultiOmics WGCNA",
        tspan("WGCNA for multi-omics")
      )
    )
  }
)
