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
  module_ui = function() {
    list(
      bigdash::bigTabItem(
        "ideograms-tab",
        EpigenomicsInputs("epigenomics"),
        create_loader("ideograms-loader")
      )
    )
  },
  module_lazy = function(PGX, ...) {
    list(
      "ideograms-tab" = list(
        ui     = function() EpigenomicsUI("epigenomics"),
        server = function() EpigenomicsBoard("epigenomics", pgx = PGX)
      )
    )
  },
  module_help = function() {
    list(
      bigdash::sidebarTabHelp(
        "ideograms-tab",
        "Beta Ideograms",
        tspan("Epigenomics visualizations and analyses for methylomics data.")
      )
    )
  }
)
