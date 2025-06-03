##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

WordCloudBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    wc_infotext <- tspan(paste("This module performs WordCloud analysis or 'keyword enrichment', i.e. it computes the enrichment of keywords for the contrasts. Frequently appearing words in the top ranked gene sets form an unbiased description of the contrast.
<br><br><br><br>
<center><iframe width='560' height='315' src='https://www.youtube.com/embed/BmPTfanUnR0?si=2irSbCjCBRgQf5Wd&amp;start=190' title='YouTube video player' frameborder='0' allow='accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share' referrerpolicy='strict-origin-when-cross-origin' allowfullscreen></iframe></center>
"), js = FALSE)

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$wc_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>WordCloud Analysis Board</strong>"),
        shiny::HTML(wc_infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    shiny::observe({
      #shiny::req(pgx$gset.meta)
      #ct <- names(pgx$gset.meta$meta)
      shiny::req(pgx)
      ct <- playbase::pgx.getContrasts(pgx)
      ct <- sort(ct[!grepl("^IA:", ct)])
      shiny::updateSelectInput(session, "wc_contrast", choices = ct)
    })

    ## ---------------------------------------------------------------
    ## ------------- Functions for WordCloud -------------------------
    ## ---------------------------------------------------------------

    getWordFreqResults <- shiny::reactive({
      shiny::req(pgx$gset.meta)
      if ("wordcloud" %in% names(pgx)) {
        res <- pgx$wordcloud
      } else {
        dbg("**** CALCULATING WORDCLOUD ****\n")
        progress <- shiny::Progress$new()
        res <- playbase::pgx.calculateWordCloud(pgx, progress = progress, pg.unit = 1)
        on.exit(progress$close())

        ## save in object??
        pgx[["wordcloud"]] <- res
      }
      return(res)
    })

    getCurrentWordEnrichment <- shiny::reactive({
      res <- getWordFreqResults()
      shiny::req(res, input$wc_contrast)

      contr <- input$wc_contrast
      gsea1 <- res$gsea[[contr]]

      # sometimes we have words that NA is tsne, make sure we remove them (likely special characters) in windows or wsl
      res$tsne <- res$tsne[!is.na(rownames(res$tsne)), ]
      res$umap <- res$umap[!is.na(rownames(res$umap)), ]

      ordered_words <- gsea1$word
      res$tsne <- res$tsne[ordered_words, ]
      res$umap <- res$umap[ordered_words, ]

      # end sometimes we have words that NA is tsne, make sure we remove them (likely special characters) in windows or wsl

      topFreq <- data.frame(gsea1, tsne = res$tsne, umap = res$umap)

      topFreq <- topFreq[order(-topFreq$NES), ]

      return(topFreq)
    })

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    # Enrichment plots

    wordcloud_plot_enrichment_server(
      "gseaplots",
      getWordFreqResults = getWordFreqResults,
      getCurrentWordEnrichment = getCurrentWordEnrichment,
      wordcloud_enrichmentTable = wordcloud_enrichmentTable,
      watermark = WATERMARK
    )

    # Word cloud

    wordcloud_plot_wordcloud_server(
      "wordcloud",
      getCurrentWordEnrichment = getCurrentWordEnrichment,
      watermark = WATERMARK
    )

    # Word t-SNE

    wordcloud_plot_wordtsne_server(
      "wordtsne",
      getCurrentWordEnrichment = getCurrentWordEnrichment,
      watermark = WATERMARK
    )

    # Enrichment table

    wordcloud_enrichmentTable <- wordcloud_table_enrichment_server(
      "wordcloud_enrichmentTable",
      getCurrentWordEnrichment = getCurrentWordEnrichment
    )

    # Leading-edge table

    wordcloud_table_leading_edge_server(
      "wordcloud_leadingEdgeTable",
      pgx = pgx,
      wc_contrast = shiny::reactive(input$wc_contrast),
      wordcloud_enrichmentTable = wordcloud_enrichmentTable,
      getCurrentWordEnrichment = getCurrentWordEnrichment
    )
  })
}
