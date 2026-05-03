## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.

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
