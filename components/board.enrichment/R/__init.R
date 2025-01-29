## make as R6 class?? e.g. add documentation, initialize object,
## object id.

MODULE.enrichment <- list(
  module_menu = function() {
    c(
      enrich = "Geneset enrichment",
      sig = "Test geneset",
      pathway = "Pathway analysis",
      wordcloud = "Word cloud"
    )
  },
  module_server = function(PGX, board_observers = NULL, labeltype = NULL, env = NULL) {
    info("[SERVER] calling SignatureBoard module")
    SignatureBoard("sig",
      pgx = PGX,
      selected_gxmethods = env$diffexpr$selected_gxmethods,
      board_observers = board_observers
    )

    info("[SERVER] calling PathwayBoard module")
    PathwayBoard("pathway",
      pgx = PGX,
      selected_gsetmethods = env$enrich$selected_gsetmethods,
      board_observers = board_observers
    )

    info("[SERVER] calling WordCloudBoard module")
    WordCloudBoard("wordcloud", pgx = PGX)
  },
  module_ui = function() {
    list(
      bigdash::bigTabItem(
        "enrich-tab",
        tagList(
          EnrichmentInputs("enrich"),
          EnrichmentUI("enrich")
        )
      ),
      bigdash::bigTabItem(
        "sig-tab",
        tagList(
          SignatureInputs("sig"),
          SignatureUI("sig")
        )
      ),
      bigdash::bigTabItem(
        "pathway-tab",
        tagList(
          PathwayInputs("pathway"),
          PathwayUI("pathway")
        )
      ),
      bigdash::bigTabItem(
        "wordcloud-tab",
        tagList(
          WordCloudInputs("wordcloud"),
          WordCloudUI("wordcloud")
        )
      )
    )
  }#,
  # module_ui2 = function() {
  #   list(
  #     list(
  #       "enrich-tab",
  #       EnrichmentUI("enrich")
  #     ),
  #     list(
  #       "sig-tab",
  #       SignatureUI("sig")
  #     ),
  #     list(
  #       "pathway-tab",
  #       PathwayUI("pathway")
  #     ),
  #     list(
  #       "wordcloud-tab",
  #       WordCloudUI("wordcloud")
  #     )
  #   )
  # }
)
