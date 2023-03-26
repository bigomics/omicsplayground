##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
correlation_plot_correlation_UMAP_ui <- function(id,
                                                    height = c(600, 800)) {
  ns <- shiny::NS(id)
  info_text <- "<b>Correlation UMAP.</b> UMAP clustering of genes using covariance as distance metric and colored by correlation (or covariance). Genes that are correlated are generally positioned close to each other. Red corresponds to positive correlation/covariance, blue for negative."
  cor_umap.opts <- shiny::tagList(
    shiny::radioButtons(ns("umap_param"), "color by:", choices = c("cor", "cov"), inline = TRUE)
  )

  PlotModuleUI(ns("plot"),
    title = "Correlation UMAP",
    plotlib = "plotly",
    label = "b",
    info.text = info_text,
    options = cor_umap.opts,
    download.fmt = c("png", "pdf", "csv"),
    height = c(700, 750),
    width = c("auto", "100%")
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
correlation_plot_correlation_UMAP_server <- function(id,
                                                        pgx,
                                                        cor_gene,
                                                        getFullGeneCorr,
                                                        getGeneCorr,
                                                        watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    # reactive function listeninng for changes in input
    cor_umap.DATA <- shiny::reactive({
      shiny::req(pgx)
      shiny::req(cor_gene)

      if (!"cluster.genes" %in% names(pgx)) {
        par(mfrow = c(1, 1))
        frame()
        text(0.5, 0.6, "Error: gene cluster position in PGX object", col = "red3")
        return(NULL)
      }

      R0 <- getFullGeneCorr()
      R1 <- getGeneCorr()

      if (is.null(R1)) {
        return(NULL)
      }
      pos <- pgx$cluster.genes$pos[["umap2d"]]
      if (input$umap_param == "cov") {
        rho0 <- R0[, "cov"]
        rho1 <- R1[, "cov"]
      } else {
        rho0 <- R0[, "cor"]
        rho1 <- R1[, "cor"]
      }
      rho0 <- rho0[match(rownames(pos), names(rho0))]
      rho1 <- rho1[match(rownames(pos), names(rho1))]
      names(rho0) <- rownames(pos)
      names(rho1) <- rownames(pos)
      ## attenuate genes *not* in current selected geneset
      ii <- which(!names(rho0) %in% names(rho1))
      rho0[ii] <- 0.2 * rho0[ii]
      return(cbind(pos, rho0, rho1))
    })

    cor_umap.PLOTFUN <- shiny::reactive({
      dt <- cor_umap.DATA()

      pos <- dt[, 1:2]
      rho0 <- dt[, 3]
      rho1 <- dt[, 4]

      gene <- cor_gene
      higenes <- c(gene)
      higenes <- names(tail(sort(rho1**2), 20))
      higenes <- unique(names(c(head(sort(rho1), 10), tail(sort(rho1), 10))))
      cexlab <- ifelse(length(higenes) == 1, 2.2, 1.3)

      p <- pgx.plotGeneUMAP(
        pgx,
        pos = pos, ## contrast=ct,
        value = rho0, title = "",
        cex = 0.9, cex.lab = cexlab,
        hilight = higenes, ntop = 20,
        plotlib = "plotly"
      )

      if (!is.null(p)) {
        return(p)
      }
    })

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = cor_umap.PLOTFUN,
      csvFunc = cor_umap.DATA,
      res = c(72, 80), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
