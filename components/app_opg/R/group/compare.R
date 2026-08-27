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
  module_lazy = function(PGX, labeltype = NULL, auth = NULL, env = NULL, reload_pgxdir = NULL, ...) {
    list(
      "isect-tab" = list(
        ui     = function() IntersectionUI("isect"),
        server = function() IntersectionBoard("isect",
          pgx = PGX,
          selected_gxmethods = env$diffexpr$selected_gxmethods,
          selected_gsetmethods = env$enrich$selected_gsetmethods
        )
      ),
      "comp-tab" = list(
        ui     = function() CompareUI("comp"),
        server = function() CompareBoard("comp",
          pgx = PGX,
          pgx_dir = reactive(auth$user_dir),
          labeltype = labeltype
        )
      ),
      "cmap-tab" = list(
        ui     = function() ConnectivityUI("cmap"),
        server = function() ConnectivityBoard("cmap",
          pgx = PGX,
          auth = auth,
          reload_pgxdir = reload_pgxdir
        )
      )
    )
  }
)
