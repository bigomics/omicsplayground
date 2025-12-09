##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

MODULE.across <- list(
  module_menu = function() {
    c(
      across = "Across datasets"
    )
  },
  module_ui = function() {
    list(
      bigdash::bigTabItem(
        "across-tab",
        AcrossInputs("across"),
        create_loader("across-loader")
      )
    )
  },
  module_ui2 = function() {
    list(
      list(
        "across-tab",
        AcrossUI("across")
      )
    )
  },
  module_server = function(PGX, labeltype = NULL, auth = NULL, env = NULL, reload_pgxdir = NULL) {
    AcrossBoard("across",
      pgx = PGX,
      pgx_dir = reactive(auth$user_dir),
      labeltype = labeltype
    )
  }
)
