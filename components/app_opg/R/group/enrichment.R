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
  module_ui = function() {
    list(
      bigdash::bigTabItem(
        "enrich-tab",
        tagList(
          EnrichmentInputs("enrich"),
          create_loader("enrich-loader")
        )
      ),
      bigdash::bigTabItem(
        "sig-tab",
        tagList(
          SignatureInputs("sig"),
          create_loader("sig-loader")
        )
      ),
      bigdash::bigTabItem(
        "pathway-tab",
        tagList(
          PathwayInputs("pathway"),
          create_loader("pathway-loader")
        )
      ),
      bigdash::bigTabItem(
        "wordcloud-tab",
        tagList(
          WordCloudInputs("wordcloud"),
          create_loader("wordcloud-loader")
        )
      )
    )
  },
  module_lazy = function(PGX, labeltype = NULL, env = NULL, ...) {
    list(
      ## "enrich-tab" = list(
      ##   ui     = function() EnrichmentUI("enrich"),
      ##   server = function() EnrichmentBoard("enrich", pgx = PGX, selected_gxmethods = env$diffexpr$selected_gxmethods)
      ## ),
      "sig-tab" = list(
        ui     = function() SignatureUI("sig"),
        server = function() SignatureBoard("sig", pgx = PGX, selected_gxmethods = env$diffexpr$selected_gxmethods)
      ),
      "pathway-tab" = list(
        ui     = function() PathwayUI("pathway"),
        server = function() PathwayBoard("pathway", pgx = PGX, selected_gsetmethods = env$enrich$selected_gsetmethods)
      ),
      "wordcloud-tab" = list(
        ui     = function() WordCloudUI("wordcloud"),
        server = function() WordCloudBoard("wordcloud", pgx = PGX)
      )
    )
  }
)
