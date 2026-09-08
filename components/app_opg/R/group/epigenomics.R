## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

## Registry for the Epigenomics menu group: which tabs it contains, their
## sidebar labels (module_menu), the shell each tab renders at boot
## (module_ui) and the UI + server bigTabsLazy() fills it with on first
## visit (module_lazy). Read by opg_ui() and opg_server().
##
## The boards themselves live in components/board.*; a group draws on
## whichever of them it needs.

MODULE.epigenomics <- list(
  module_menu = function() {
    c(
      ideograms = "Beta Ideograms"
    )
  },
  module_ui = function(ns = identity) {
    list(
      bigdash::bigTabItem(
        ns("ideograms-tab"),
        EpigenomicsInputs(ns("epigenomics")),
        create_loader(ns("ideograms-loader"))
      )
    )
  },
  module_lazy = function(PGX, ns = identity, ...) {
    list(
      "ideograms-tab" = list(
        ui     = function() EpigenomicsUI(ns("epigenomics")),
        server = function() EpigenomicsBoard("epigenomics", pgx = PGX)
      )
    )
  },
  module_help = function(ns = identity) {
    list(
      bigdash::sidebarTabHelp(
        ns("ideograms-tab"),
        "Beta Ideograms",
        tspan("Epigenomics visualizations and analyses for methylomics data.")
      )
    )
  }
)
