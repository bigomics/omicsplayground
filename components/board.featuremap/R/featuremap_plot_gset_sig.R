##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

featuremap_plot_gset_sig_ui <- function(
    id,
    label = "",
    title,
    info.text,
    info.methods,
    info.extra_link,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  info_text <- "<b>Geneset signature maps.</b> UMAP clustering of genes colored by relative log-expression of the phenotype group. The distance metric is covariance. Genes that are clustered nearby have high covariance."

  PlotModuleUI(
    ns("gset_sig"),
    title = title,
    label = "b",
    info.text = info_text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

featuremap_plot_gset_sig_server <- function(id,
                                            pgx,
                                            getGsetUMAP,
                                            sigvar,
                                            ref_group,
                                            plotFeaturesPanel,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)

      pos <- pgx$cluster.gsets$pos[["umap2d"]]
      hilight <- NULL
      pheno <- "tissue"
      pheno <- sigvar()
      if (any(pheno %in% colnames(pgx$samples))) {
        y <- pgx$samples[, pheno]
        ref <- ref_group()
        if (ref == "<average>") {
          refX <- rowMeans(pgx$gsetX, na.rm = TRUE)
        } else {
          kk <- which(y == ref)
          refX <- rowMeans(pgx$gsetX[, kk], na.rm = TRUE)
        }
        X <- pgx$gsetX - refX
        F <- do.call(cbind, tapply(1:ncol(X), y, function(i) {
          rowMeans(X[, i, drop = FALSE], na.rm = TRUE)
        }))
      } else {
        F <- playbase::pgx.getMetaMatrix(pgx, level = "geneset")$fc
        kk <- intersect(pheno, colnames(F))
        F <- F[, kk, drop = FALSE]
      }
      if (nrow(F) == 0) {
        return(NULL)
      }
      return(list(F, pos))
    })

    renderPlots <- function() {
      dt <- plot_data()
      F <- dt[[1]]
      pos <- dt[[2]]
      shiny::req(F, pos)

      kk <- intersect(rownames(pos), rownames(F))
      F <- F[kk, , drop = FALSE]
      pos <- pos[kk, , drop = FALSE]

      ntop <- 15
      nc <- ceiling(sqrt(1.33 * ncol(F)))
      nr <- ceiling(ncol(F) / nc)

      par(mfrow = c(nr, nc), mar = c(3, 1, 1, 0.5), mgp = c(1.6, 0.55, 0))
      progress <- NULL
      if (!interactive()) {
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Computing feature plots...", value = 0)
      }
      plotFeaturesPanel(pos, F, ntop, nr, nc, sel = NULL, progress)
    }

    PlotModuleServer(
      "gset_sig",
      func = renderPlots,
      csvFunc = plot_data,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 90),
      add.watermark = watermark
    )
  })
}
