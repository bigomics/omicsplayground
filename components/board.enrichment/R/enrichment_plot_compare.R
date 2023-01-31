##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

enrichment_plot_compare_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "Under the <strong>Compare</strong> tab, enrichment profiles of the selected geneset in enrichment Table <code>I</code> can be visualised against all available contrasts."

  PlotModuleUI(
    ns("plot"),
    title = "Enrichment of geneset across multiple contrasts",
    info.text = info_text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

enrichment_plot_compare_server <- function(id,
                                           inputData,
                                           gs_contrast,
                                           gset_selected,
                                           selected_gsetmethods,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    compare.RENDER <- shiny::reactive({
      ngs <- inputData()
      shiny::req(ngs, gs_contrast())

      comp <- 1
      comp <- gs_contrast()
      if (is.null(comp)) {
        return(NULL)
      }

      gset <- rownames(ngs$gsetX)[1]
      gset <- gset_selected()
      if (is.null(gset)) {
        return(NULL)
      }
      gset <- gset[1]

      score <- sapply(ngs$gset.meta$meta, function(x) x[gset, "meta.fx"])

      top.up <- names(sort(score[which(score > 0)], decreasing = TRUE))
      top.dn <- names(sort(score[which(score < 0)]))
      genes <- names(which(ngs$GMT[, gset] != 0))
      genes <- toupper(sub(".*:", "", genes))
      gx.meta <- ngs$gx.meta$meta

      gsmethods <- selected_gsetmethods()

      par(mfrow = c(2, 5), mar = c(0.5, 3.2, 2.6, 0.5), mgp = c(2, 0.8, 0))
      i <- 1
      for (i in 1:5) {
        if (i > length(top.up)) {
          frame()
        } else {
          cmp <- top.up[i]
          rnk0 <- gx.meta[[cmp]]$meta.fx
          names(rnk0) <- rownames(gx.meta[[1]])
          names(rnk0) <- toupper(sub(".*:", "", names(rnk0)))

          gs.meta <- ngs$gset.meta$meta[[cmp]]
          qv0 <- max(gs.meta[gset, "q"][, gsmethods], na.rm = TRUE)

          gs1 <- breakstring(gset, 28, 50, force = FALSE)
          cmp <- paste0(gset, "\n@", cmp)
          gsea.enplot(rnk0, genes,
            names = NULL, ## main=gs,
            main = cmp, xlab = "",
            cex.main = 0.80, len.main = 72
          )
          qv1 <- formatC(qv0, format = "e", digits = 2)
          legend("topright", paste("q=", qv1), bty = "n", cex = 0.85)
        }
      }
      for (i in 1:5) {
        if (i > length(top.dn)) {
          frame()
        } else {
          cmp <- top.dn[i]
          rnk0 <- gx.meta[[cmp]]$meta.fx
          names(rnk0) <- rownames(gx.meta[[1]])
          names(rnk0) <- toupper(sub(".*:", "", names(rnk0)))

          gs.meta <- ngs$gset.meta$meta[[cmp]]
          qv0 <- max(gs.meta[gset, "q"][, gsmethods], na.rm = TRUE)

          gs1 <- breakstring(gset, 28, 50, force = FALSE)
          cmp <- paste0(gset, "\n@", cmp)
          gsea.enplot(rnk0, genes,
            names = NULL, ## main=gs,
            main = cmp, xlab = "",
            cex.main = 0.80, len.main = 72
          )
          qv1 <- formatC(qv0, format = "e", digits = 2)
          legend("topright", paste("q=", qv1), bty = "n", cex = 0.85)
        }
      }
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = compare.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(95, 100),
      add.watermark = watermark
    )
  })
}
