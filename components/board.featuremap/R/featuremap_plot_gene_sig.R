##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

featuremap_plot_gene_sig_ui <- function(
    id,
    label = "",
    title,
    caption,
    info.text,
    info.methods,
    info.extra_link,
    height,
    width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("gene_sig"),
    title = title,
    label = "b",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

featuremap_plot_gene_sig_server <- function(id,
                                            pgx,
                                            sigvar,
                                            ref_group,
                                            plotFeaturesPanel,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      pheno <- sigvar()
      ref <- ref_group()

      pos <- pgx$cluster.genes$pos[["umap2d"]]
      if (any(pheno %in% colnames(pgx$samples))) {
        y <- pgx$samples[, pheno]
        if (ref == "<average>") {
          refX <- rowMeans(pgx$X, na.rm = TRUE)
        } else {
          kk <- which(y == ref)
          refX <- rowMeans(pgx$X[, kk], na.rm = TRUE)
        }
        X <- pgx$X - refX
        y <- pgx$samples[, pheno]
        F <- do.call(cbind, tapply(1:ncol(X), y, function(i) {
          rowMeans(X[, i, drop = FALSE], na.rm = TRUE)
        }))
      } else {
        F <- playbase::pgx.getMetaMatrix(pgx, level = "gene")$fc
        kk <- intersect(pheno, colnames(F))
        F <- F[, kk, drop = FALSE]
      }
      if (nrow(F) == 0 || NCOL(F) == 0) {
        return(NULL)
      }
      return(list(F, pos))
    })

    renderPlots <- function() {
      dt <- plot_data()

      F <- dt[[1]]
      pos <- dt[[2]]
      shiny::req(F, pos)

      nc <- ceiling(sqrt(1.33 * ncol(F)))
      nr <- ceiling(ncol(F) / nc)

      par(mfrow = c(nr, nc), mar = c(2, 1, 1, 0), mgp = c(1.6, 0.55, 0), las = 0)
      progress <- NULL
      if (!interactive()) {
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Computing feature plots...", value = 0)
      }
      plotFeaturesPanel(pos, F, ntop = ntop, nr, nc, sel = NULL, progress)
    }

    PlotModuleServer(
      "gene_sig",
      plotlib = "base",
      func = renderPlots,
      csvFunc = plot_data,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 90),
      add.watermark = watermark
    )
  })
}
