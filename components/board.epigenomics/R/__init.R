## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

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
  module_ui2 = function() {
    list(
      list(
        "ideograms-tab",
        EpigenomicsUI("epigenomics")
      )
    )
  },
  ## Per-board servers, keyed by tab, so opg_server can start just the board
  ## that was opened. module_server() above stays for the eager path.
  module_servers = list(
    "ideograms-tab" = function(PGX)
      EpigenomicsBoard("epigenomics", pgx = PGX)
  ),
  module_server = function(PGX) {
    EpigenomicsBoard("epigenomics", pgx = PGX)
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
