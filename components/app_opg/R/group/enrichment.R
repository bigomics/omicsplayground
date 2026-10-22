## make as R6 class?? e.g. add documentation, initialize object,
## object id.

## Registry for the GeneSets menu group: which tabs it contains, their
## sidebar labels (module_menu), the shell each tab renders at boot
## (module_ui) and the UI + server bigTabsLazy() fills it with on first
## visit (module_lazy). Read by opg_ui() and opg_server().
##
## The boards themselves live in components/board.*; a group draws on
## whichever of them it needs.

MODULE.enrichment <- list(
  module_menu = function() {
    c(
      enrich = "Geneset enrichment",
      sig = "Test geneset",
      pathway = "Pathway analysis",
      wordcloud = "Word cloud"
    )
  },
  module_ui = function(ns = identity) {
    list(
      bigdash::bigTabItem(
        ns("enrich-tab"),
        tagList(
          EnrichmentInputs(ns("enrich")),
          create_loader(ns("enrich-loader"))
        )
      ),
      bigdash::bigTabItem(
        ns("sig-tab"),
        tagList(
          SignatureInputs(ns("sig")),
          create_loader(ns("sig-loader"))
        )
      ),
      bigdash::bigTabItem(
        ns("pathway-tab"),
        tagList(
          PathwayInputs(ns("pathway")),
          create_loader(ns("pathway-loader"))
        )
      ),
      bigdash::bigTabItem(
        ns("wordcloud-tab"),
        tagList(
          WordCloudInputs(ns("wordcloud")),
          create_loader(ns("wordcloud-loader"))
        )
      )
    )
  },
  module_lazy = function(PGX, labeltype = NULL, env = NULL, ns = identity, ...) {
    list(
      ## "enrich-tab" = list(
      ##   ui     = function() EnrichmentUI(ns("enrich")),
      ##   server = function() EnrichmentBoard("enrich", pgx = PGX, selected_gxmethods = env$diffexpr$selected_gxmethods)
      ## ),
      "sig-tab" = list(
        ui     = function() SignatureUI(ns("sig")),
        server = function() SignatureBoard("sig", pgx = PGX, selected_gxmethods = env$diffexpr$selected_gxmethods)
      ),
      "pathway-tab" = list(
        ui     = function() PathwayUI(ns("pathway")),
        server = function() PathwayBoard("pathway", pgx = PGX, selected_gsetmethods = env$enrich$selected_gsetmethods)
      ),
      "wordcloud-tab" = list(
        ui     = function() WordCloudUI(ns("wordcloud")),
        server = function() WordCloudBoard("wordcloud", pgx = PGX)
      )
    )
  }
)
