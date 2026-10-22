##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## Registry for the Compare menu group: which tabs it contains, their
## sidebar labels (module_menu), the shell each tab renders at boot
## (module_ui) and the UI + server bigTabsLazy() fills it with on first
## visit (module_lazy). Read by opg_ui() and opg_server().
##
## The boards themselves live in components/board.*; a group draws on
## whichever of them it needs.

MODULE.compare <- list(
  module_menu = function() {
    c(
      isect = "Compare signatures",
      comp = "Compare datasets",
      cmap = "Similar experiments"
    )
  },
  module_ui = function(ns = identity) {
    list(
      bigdash::bigTabItem(
        ns("isect-tab"),
        IntersectionInputs(ns("isect")),
        create_loader(ns("isect-loader"))
      ),
      bigdash::bigTabItem(
        ns("comp-tab"),
        CompareInputs(ns("comp")),
        create_loader(ns("comp-loader"))
      ),
      bigdash::bigTabItem(
        ns("cmap-tab"),
        ConnectivityInputs(ns("cmap")),
        create_loader(ns("cmap-loader"))
      )
    )
  },
  module_lazy = function(PGX, labeltype = NULL, auth = NULL, env = NULL, reload_pgxdir = NULL, ns = identity, ...) {
    list(
      "isect-tab" = list(
        ui     = function() IntersectionUI(ns("isect")),
        server = function() IntersectionBoard("isect",
          pgx = PGX,
          selected_gxmethods = env$diffexpr$selected_gxmethods,
          selected_gsetmethods = env$enrich$selected_gsetmethods
        )
      ),
      "comp-tab" = list(
        ui     = function() CompareUI(ns("comp")),
        server = function() CompareBoard("comp",
          pgx = PGX,
          pgx_dir = reactive(auth$user_dir),
          labeltype = labeltype
        )
      ),
      "cmap-tab" = list(
        ui     = function() ConnectivityUI(ns("cmap")),
        server = function() ConnectivityBoard("cmap",
          pgx = PGX,
          auth = auth,
          reload_pgxdir = reload_pgxdir
        )
      )
    )
  }
)
