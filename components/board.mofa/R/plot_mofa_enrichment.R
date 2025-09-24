##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_enrichment_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::selectInput(ns("labeltype"), "Label type",
      choices = c("feature", "symbol", "gene_title")
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

mofa_plot_enrichment_server <- function(id,
                                        mofa,
                                        pgx,
                                        input_k = reactive(1),
                                        req.selection = FALSE,
                                        select = reactive(NULL),
                                        ntop = 15,
                                        watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function(ntop = 16, par = TRUE) {
      mofa <- mofa()
      gsea <- mofa$gsea
      validate(need(!is.null(gsea), "missing MOFA data."))

      k <- input_k()
      shiny::req(k)
      shiny::req(k %in% names(gsea$table))
      sel <- select()

      nx <- grep("[:]", rownames(mofa$W), value = TRUE)
      dtypes <- unique(sub(":.*", "", nx))

      if (length(sel) == 1) {
        W <- playbase::rename_by2(mofa$W, pgx$genes, "symbol", keep.prefix = TRUE)
        rnk <- playbase::normalize_multirank(W[, k])

        gset <- sub(".*:", "", names(which(pgx$GMT[, sel] != 0)))
        gset <- as.vector(sapply(dtypes, paste0, ":", gset))
        gset <- intersect(gset, names(rnk))

        par(mfrow = c(1, 1), mar = c(3, 4, 2, 1))
        playbase::gsea.enplot(
          rnk, gset,
          cex = 1.3,
          main = sel, cex.main = 1.2, main.line = 0.6,
          ylab = "Normalized rank metric",
          cex.lab = 1, lab.line = c(1, 2.4)
        )
      } else {
        if (par) {
          par(mfrow = c(1, 2), mar = c(4, 7, 0, 2))
          plot.new()
        }
        playbase::mofa.plot_enrichment(
          gsea$table[[k]],
          type = "barplot",
          ntop = ntop,
          select = select(),
          strip.names = TRUE,
          sort = FALSE,
          par = FALSE,
          title = ""
        )
      }
    }

    plot.RENDER2 <- function() {
      par(mfrow = c(1, 2), mar = c(4, 14, 1, 2))
      plot.new()
      plot.RENDER(ntop = 30, par = FALSE)
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,
      pdf.width = 9, pdf.height = 5,
      res = c(72, 110),
      add.watermark = watermark
    )
  })
}
