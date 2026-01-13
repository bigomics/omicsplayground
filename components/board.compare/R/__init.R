##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

MODULE.compare <- list(
  module_menu = function() {
    c(
      isect = "Compare signatures",
      comp = "Compare datasets",
      cmap = "Similar experiments"
    )
  },
  module_ui = function() {
    list(
      bigdash::bigTabItem(
        "isect-tab",
        IntersectionInputs("isect"),
        create_loader("isect-loader")
      ),
      bigdash::bigTabItem(
        "comp-tab",
        CompareInputs("comp"),
        create_loader("comp-loader")
      ),
      bigdash::bigTabItem(
        "cmap-tab",
        ConnectivityInputs("cmap"),
        create_loader("cmap-loader")
      )
    )
  },
  module_ui2 = function() {
    list(
      list(
        "isect-tab",
        IntersectionUI("isect")
      ),
      list(
        "comp-tab",
        CompareUI("comp")
      ),
      list(
        "cmap-tab",
        ConnectivityUI("cmap")
      )
    )
  },
  module_server = function(PGX, labeltype = NULL, auth = NULL, env = NULL, reload_pgxdir = NULL) {
    IntersectionBoard("isect",
      pgx = PGX,
      selected_gxmethods = env$diffexpr$selected_gxmethods,
      selected_gsetmethods = env$enrich$selected_gsetmethods
    )

    CompareBoard(
      "comp",
      pgx = PGX,
      pgx_dir = reactive(auth$user_dir),
      labeltype = labeltype
    )

    ConnectivityBoard("cmap",
      pgx = PGX,
      auth = auth,
      reload_pgxdir = reload_pgxdir
    )
  }
)
