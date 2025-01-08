##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_compare_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    caption = caption,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

enrichment_plot_compare_server <- function(id,
                                           pgx,
                                           gs_contrast,
                                           gset_selected,
                                           selected_gsetmethods,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    render_compare <- function() {
      shiny::req(pgx$X, gs_contrast())

      comp <- 1
      comp <- gs_contrast()
      if (is.null(comp)) {
        return(NULL)
      }

      gset <- rownames(pgx$gsetX)[1]
      gset <- gset_selected()
      if (is.null(gset)) {
        frame()
        text(0.5, 0.5, tspan("Please select a geneset", js = FALSE), col = "grey50")
        return()
      }
      gset <- gset[1]

      score <- sapply(pgx$gset.meta$meta, function(x) x[gset, "meta.fx"])

      top.up <- names(sort(score[which(score > 0)], decreasing = TRUE))
      top.dn <- names(sort(score[which(score < 0)]))
      genes <- names(which(pgx$GMT[, gset] != 0))
      genes <- toupper(sub(".*:", "", genes))
      gx.meta <- pgx$gx.meta$meta

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

          gs.meta <- pgx$gset.meta$meta[[cmp]]
          qv0 <- max(gs.meta[gset, "q"][, gsmethods], na.rm = TRUE)

          gs1 <- playbase::breakstring(gset, 28, 50, force = FALSE)
          cmp <- paste0(gset, "\n@", cmp)
          
          playbase::gsea.enplot(
            rnk0,
            genes,
            names = NULL,
            main = cmp,
            xlab = "",
            cex.main = 0.80,
            len.main = 72
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

          gs.meta <- pgx$gset.meta$meta[[cmp]]
          qv0 <- max(gs.meta[gset, "q"][, gsmethods], na.rm = TRUE)

          gs1 <- playbase::breakstring(gset, 28, 50, force = FALSE)
          cmp <- paste0(gset, "\n@", cmp)

          playbase::gsea.enplot(
            rnk0,
            genes,
            names = NULL,
            main = cmp,
            xlab = "",
            cex.main = 0.80,
            len.main = 72
          )
          qv1 <- formatC(qv0, format = "e", digits = 2)
          legend("topright", paste("q=", qv1), bty = "n", cex = 0.85)
        }
      }

    }

    compare.RENDER <- function() { render_compare() }

    compare.RENDER2 <- function() { render_compare() }

    PlotModuleServer(
      "plot",
      func = compare.RENDER,
      func2 = compare.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(95, 100),
      add.watermark = watermark
    )
  })
}
