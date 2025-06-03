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
correlation_plot_correlation_UMAP_ui <- function(
    id,
    title,
    caption,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    height,
    width) {
  ns <- shiny::NS(id)

  cor_umap.opts <- shiny::tagList(
    shiny::radioButtons(ns("umap_param"), "color by:", choices = c("cor", "cov"), inline = TRUE)
  )

  PlotModuleUI(ns("plot"),
    title = title,
    caption = caption,
    plotlib = "plotly",
    label = "b",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    options = cor_umap.opts,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width
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
                                                     gene,
                                                     getFullGeneCorr,
                                                     getGeneCorr,
                                                     watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    # reactive function listeninng for changes in input
    cor_umap.DATA <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(gene())

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
      R0 <- R0[intersect(rownames(pos), rownames(R0)), , drop = FALSE]
      R1 <- R1[intersect(rownames(pos), rownames(R1)), , drop = FALSE]
      pos <- pos[intersect(rownames(pos), rownames(R1)), , drop = FALSE]

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

      gene <- gene()
      higenes <- c(gene)
      higenes <- names(tail(sort(rho1**2), 20))
      higenes <- unique(names(c(head(sort(rho1), 10), tail(sort(rho1), 10))))
      cexlab <- ifelse(length(higenes) == 1, 2.2, 1.3)

      p <- playbase::pgx.plotGeneUMAP(
        pgx,
        pos = pos,
        value = rho0,
        title = "",
        cex = 0.9,
        cex.lab = cexlab,
        hilight = higenes,
        ntop = 20,
        labeltype = "gene_name",
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
