##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

featuremap_plot_gset_sig_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)

  info_text <- "<b>Geneset signature maps.</b> UMAP clustering of genes colored by relative log-expression of the phenotype group. The distance metric is covariance. Genes that are clustered nearby have high covariance."

  PlotModuleUI(
    ns("gset_sig"),
    title = "Geneset Signatures",
    label = "b",
    info.text = info_text,
    height = c(600, 700),
    width = c("auto", "100%"),
    download.fmt = c("png", "pdf")
  )
}

featuremap_plot_gset_sig_server <- function(id,
                                            pgx,
                                            getGsetUMAP,
                                            sigvar,
                                            plotFeaturesPanel,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    gsetSigPlots.plot_data <- shiny::reactive({
      shiny::req(pgx)

      pos <- getGsetUMAP()
      hilight <- NULL

      pheno <- "tissue"
      pheno <- sigvar()
      if (pheno %in% colnames(pgx$samples)) {
        X <- pgx$gsetX - rowMeans(pgx$gsetX)
        y <- pgx$samples[, pheno]
        F <- do.call(cbind, tapply(1:ncol(X), y, function(i) {
          rowMeans(X[, i, drop = FALSE])
        }))
      } else {
        F <- pgx.getMetaMatrix(pgx, level = "geneset")$fc
      }
      if (nrow(F) == 0) {
        return(NULL)
      }
      return(list(F, pos))
    })

    gsetSigPlots.RENDER <- shiny::reactive({
      dt <- gsetSigPlots.plot_data()
      F <- dt[[1]]
      pos <- dt[[2]]

      ntop <- 15
      nc <- ceiling(sqrt(ncol(F)))
      nr <- ceiling(ncol(F) / nc)
      nr2 <- ifelse(nr <= 3, nc, nr)
      par(mfrow = c(nr2, nc), mar = c(3, 2, 1.5, 2), mgp = c(1.6, 0.6, 0))
      par(mfrow = c(nr2, nc), mar = c(3, 1, 1, 0.5), mgp = c(1.6, 0.55, 0))
      progress <- NULL
      if (!interactive()) {
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Computing feature plots...", value = 0)
      }
      plotFeaturesPanel(pos, F, ntop, nr, nc, sel = NULL, progress)
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "gset_sig",
      func = gsetSigPlots.RENDER,
      csvFunc = gsetSigPlots.plot_data,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 90),
      add.watermark = watermark
    )
  })
}
