##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wordcloud_plot_enrichment_ui <- function(
  id,
  title,
  info.text,
  caption,
  height) {
  ns <- shiny::NS(id)

  gseaplots_opts <- shiny::tagList(
    withTooltip(shiny::textInput(ns("gseaplots_keywords"), "Keyword:", "cell cycle"),
      "Paste a keyword such as 'apoptosis', 'replication' or 'cell cycle'.",
      placement = "top", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = "a",
    info.text = info.text,
    caption = caption,
    options = gseaplots_opts,
    height = height,
    download.fmt = c("png", "pdf")
  )
}

wordcloud_plot_enrichment_server <- function(id,
                                             getWordFreqResults,
                                             getCurrentWordEnrichment,
                                             wordcloud_enrichmentTable,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    gseaplots.RENDER <- shiny::reactive({
      res <- getWordFreqResults()
      topFreq <- getCurrentWordEnrichment()

      shiny::req(res, topFreq)

      S <- res$S ## geneset expressions
      keyword <- "ribosomal"
      keyword <- "lipid"
      keyword <- "apoptosis"
      keyword <- "cell.cycle"

      sel.row <- wordcloud_enrichmentTable$rows_selected()
      shiny::req(sel.row)
      keyword <- topFreq$word[sel.row]

      if (length(keyword) == 0 || keyword[1] %in% c(NA, "")) keyword <- "cell.cycle"
      targets <- names(which(res$W[, keyword] == 1))

      nes <- unlist(sapply(res[["gsea"]], function(G) G[match(keyword, G$word), "NES"]))
      pv <- unlist(sapply(res[["gsea"]], function(G) G[match(keyword, G$word), "pval"]))
      qv <- unlist(sapply(res[["gsea"]], function(G) G[match(keyword, G$word), "padj"]))
      names(qv) <- names(pv) <- names(nes) <- sub("[.]NES", "", names(nes))

      top <- names(pv)[order(-abs(nes), pv)]
      top <- intersect(top, colnames(S))

      par(mfrow = c(3, 3), mar = c(0.2, 3.2, 3.2, 0.2), mgp = c(1.8, 0.7, 0))
      i <- 1
      for (i in 1:9) {
        if (i > length(top)) {
          frame()
        } else {
          a <- top[i]
          playbase::gsea.enplot(S[, a], targets,
            names = NULL, ## main=gs,
            main = paste0("#", toupper(keyword), "\n@", a),
            cex.main = 0.9, len.main = 80, xlab = ""
          )
          qv1 <- formatC(qv[a], format = "e", digits = 3)
          nes1 <- formatC(nes[a], format = "f", digits = 3)
          tt <- c(paste("NES=", nes1), paste("q=", qv1))
          legend("topright", tt, bty = "n", cex = 0.85)
        }
      }
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = gseaplots.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = 90,
      add.watermark = watermark
    )
  })
}
