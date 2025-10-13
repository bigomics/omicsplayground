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
signature_plot_markers_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption,
    height) {
  ns <- shiny::NS(id)

  markers.opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("markers_sortby"), "Sort by:",
        choices = c("correlation", "probability", "name"), inline = TRUE
      ),
      "Sort by correlation, probability or name.",
      placement = "top",
      options = list(container = "body")
    ),
    withTooltip(
      shiny::radioButtons(ns("markers_layout"), "Layout:",
        choices = c("4x4", "6x6"),
        inline = TRUE
      ),
      "Choose layout.",
      placement = "top", options = list(container = "body")
    ),
  )

  PlotModuleUI(
    id = ns("plot"),
    plotlib = "plotly",
    title = title,
    caption = caption,
    options = markers.opts,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    download.fmt = c("png", "pdf", "svg"),
    height = height
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
signature_plot_markers_server <- function(id,
                                          pgx,
                                          getCurrentMarkers,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    calcSingleSampleValues <- function(X, y, method = c("rho", "gsva")) {
      ## Calculates single-sample enrichment values for given matrix and
      ## binarized signature vector.
      ## very fast rank difference

      if (is.null(names(y)) && length(y) != nrow(X)) {
        cat("<signature:calcSingleSampleValues> FATAL ERROR: y must be named if not matched\n")
        return(NULL)
      }

      if (!is.null(names(y)) && length(y) != nrow(X)) {
        y <- y[match(rownames(X), names(y))]
      }

      names(y) <- rownames(X)
      jj <- which(!is.na(y))
      X <- X[jj, ]
      y <- y[jj]

      if (sum(y != 0) == 0) {
        cat("<signature:calcSingleSampleValues> WARNING: y is all zero!\n")
        matzero <- matrix(0, nrow = ncol(X), ncol = length(method))
        colnames(matzero) <- method
        rownames(matzero) <- colnames(X)
        return(matzero)
      }

      ss.rank <- function(x) scale(sign(x) * rank(abs(x)), center = FALSE)[, 1]

      S <- list()

      if ("rho" %in% method) {
        S[["rho"]] <- cor(apply(X, 2, ss.rank), y, use = "pairwise")[, 1]
      }

      if ("gsva" %in% method) {
        gset <- names(y)[which(y != 0)]
        gmt <- list("gmt" = gset)

        new.gsva <- exists("gsvaParam", where = asNamespace("GSVA"), mode = "function")

        if (new.gsva) {
          res.gsva <- GSVA::gsva(GSVA::gsvaParam(X, gmt, maxDiff = TRUE)) ## new style :(
        } else {
          res.gsva <- GSVA::gsva(X, gmt, method = "gsva", parallel.sz = 1) ## old style...
        }

        res.colnames <- colnames(res.gsva)
        fc <- as.vector(res.gsva[1, ])
        names(fc) <- res.colnames
        S[["gsva"]] <- fc[colnames(X)]
      }

      if (length(S) > 1) {
        S1 <- do.call(cbind, S)
      } else {
        S1 <- S[[1]]
      }

      S1 <- as.matrix(S1)
      rownames(S1) <- colnames(X)
      colnames(S1) <- names(S)

      return(S1)
    }

    getSingleSampleEnrichment <- shiny::reactive({
      ## Calls calcSingleSampleValues() and calculates single-sample
      ## enrichment values for complete data matrix and reduced data by
      ## group (for currentmarkers)

      if (is.null(pgx$X)) {
        return(NULL)
      }

      ## get the signature
      markers <- getCurrentMarkers()
      this.gset <- markers$symbols
      if (is.null(this.gset)) {
        return(NULL)
      }

      X <- grp <- NULL
      is.sc <- !is.null(pgx$datatype) && pgx$datatype == "scRNAseq"
      if (is.sc) {
        message("[Test signatures: Markers] Computing supercells.")
        ct <- pgx$samples[, "celltype"]
        block.group <- paste0(ct, ":", apply(pgx$contrasts, 1, paste, collapse = "_"))
        if ("batch" %in% colnames(pgx$samples)) {
          block.group <- paste0(block.group, ":", pgx$samples[, "batch"])
        }
        nb <- round(ncol(pgx$counts) / 1000) # strong compression
        message("[pgx.wgcna] running SuperCell. nb = ", nb)
        sc <- playbase::pgx.supercell(pmax(2**pgx$X - 1, 0), pgx$samples, block.group, nb)
        message("[pgx.wgcna] SuperCell done: ", ncol(pgx$counts), " ->", ncol(sc$counts))
        message("[pgx.wgcna] Normalizing supercell matrix (logCPM)")
        X <- as.matrix(playbase::logCPM(sc$counts, total = 1e4, prior = 1))
        samples <- sc$meta
        grp <- samples[, "celltype"]
        remove(ct, block.group, nb, sc)
        gc()
      }

      X <- pgx$X
      if (any(is.na(X))) X <- playbase::imputeMissing(X, method = "SVD2")
      xsymbol <- pgx$genes[rownames(X), "symbol"]
      X <- playbase::rowmean(X, xsymbol)
      y <- 1 * (rownames(X) %in% this.gset)
      names(y) <- rownames(X)

      ## expression by group
      if (is.null(grp)) grp <- pgx$model.parameters$group
      groups <- unique(grp)
      gX <- sapply(groups, function(g) rowMeans(X[, which(grp == g), drop = FALSE], na.rm = TRUE))
      colnames(gX) <- groups

      ## for large datasets pre-grouping is faster
      ss.bygroup <- calcSingleSampleValues(gX, y, method = c("rho", "gsva"))

      ## by sample is slow... so no gsva
      ss1 <- calcSingleSampleValues(X, y, method = "rho")
      ss.bysample <- cbind(rho = ss1)

      res <- list(by.sample = ss.bysample, by.group = ss.bygroup)

      return(res)
    })

    get_plots <- function() {
      ## get markers
      markers <- getCurrentMarkers()
      shiny::req(markers)

      ## get GSVA values
      res <- getSingleSampleEnrichment()
      shiny::req(res)

      features <- markers$features
      # make sure features are in pgx$X
      features <- intersect(features, rownames(pgx$X))
      gx <- pgx$X[features, , drop = FALSE]
      if (nrow(gx) == 0) {
        cat("WARNING:: Markers:: markers do not match!!\n")
        return(NULL)
      }

      ## get t-SNE positions of samples
      pos <- pgx$tsne2d[colnames(gx), ]
      gx <- gx - min(gx, na.rm = TRUE) + 0.001 ## subtract background
      grp <- pgx$model.parameters$group
      zx <- t(apply(gx, 1, function(x) tapply(x, as.character(grp), mean)))
      gx <- gx[order(-apply(zx, 1, sd, na.rm = TRUE)), , drop = FALSE]
      rownames(gx) <- sub(".*:", "", rownames(gx))

      ## get GSVA values and make some non-linear value fc1
      S <- res$by.sample
      if (NCOL(S) == 1) {
        fc <- S[, 1]
      } else {
        fc <- colMeans(t(S) / (1e-8 + sqrt(colSums(S**2)))) ## scaled mean
      }
      fc <- scale(fc)[, 1] ## scale??
      names(fc) <- rownames(S)
      fc1 <- tanh(1.0 * fc / (1e-4 + sd(fc)))
      fc1 <- fc1[rownames(pos)]

      cex1 <- 1.2
      cex1 <- 0.7 * c(1.6, 1.2, 0.8, 0.5)[cut(nrow(pos), breaks = c(-1, 40, 200, 1000, 1e10))]

      nmax <- NULL
      if (input$markers_layout == "6x6") nmax <- 35
      if (input$markers_layout == "4x4") nmax <- 15

      top.gx <- head(gx, nmax)
      if (input$markers_sortby == "name") {
        top.gx <- top.gx[order(rownames(top.gx)), , drop = FALSE]
      }
      if (input$markers_sortby == "probability") {
        top.gx <- top.gx[order(-rowMeans(top.gx, na.rm = TRUE)), , drop = FALSE]
      }
      if (input$markers_sortby == "correlation") {
        rho <- cor(t(top.gx), fc1)[, 1]
        top.gx <- top.gx[order(-rho), , drop = FALSE]
      }

      plt <- list()
      i <- 0
      for (i in 0:min(nmax, nrow(top.gx))) {
        jj <- 1:ncol(top.gx)
        if (i == 0) {
          klrpal <- grDevices::colorRampPalette(c(omics_colors("brand_blue"), omics_colors("grey"), omics_colors("red")))(16)
          colvar <- fc1
          klr1 <- klrpal[8 + round(7 * fc1)]
          tt <- "INPUT SIGNATURE"
          jj <- order(abs(fc1))
        } else {
          klrpal <- colorRampPalette(c("grey90", "grey60", omics_colors("red")))(16)
          colvar <- pmax(top.gx[i, ], 0)
          colvar <- 1 + round(15 * (colvar / (0.7 * max(colvar, na.rm = TRUE) + 0.3 * max(top.gx, na.rm = TRUE))))
          klr1 <- klrpal[colvar]
          gene <- substring(sub(".*:", "", rownames(top.gx)[i]), 1, 80)
          tt <- playbase::breakstring(gene, n = 20, force = TRUE)
          jj <- order(abs(top.gx[i, ]))
        }
        klr1 <- paste0(gplots::col2hex(klr1), "99")

        ## ------- start plot ----------
        p <- playbase::pgx.scatterPlotXY.PLOTLY(
          pos[jj, ],
          var = colvar[jj],
          col = klrpal,
          cex = 1.0 * cex1,
          xlab = "",
          ylab = "",
          xlim = 1.2 * range(pos[, 1]),
          ylim = 1.2 * range(pos[, 2]),
          axis = FALSE,
          title = tt,
          cex.title = 0.85,
          title.y = 0.86,
          label.clusters = FALSE,
          legend = FALSE,
          gridcolor = "fff"
        ) %>% plotly::layout(
          plot_bgcolor = "#f8f8f8"
        )
        plt[[i + 1]] <- p
      }
      return(plt)
    }


    plotly.RENDER <- function() {
      plt <- get_plots()
      shiny::req(plt)
      nr <- ceiling(sqrt(length(plt)))

      fig <- plotly::subplot(
        plt,
        nrows = nr,
        margin = 0.01
      ) %>%
        plotly_default() %>%
        plotly::layout(
          title = list(text = tspan("genes in signature", js = FALSE), size = 12),
          margin = list(l = 0, r = 0, b = 0, t = 30) # lrbt
        )
      return(fig)
    }

    plotly.RENDER_MODAL <- function() {
      fig <- plotly.RENDER() %>%
        plotly_modal_default() %>%
        plotly::layout(
          margin = list(l = 0, r = 0, b = 0, t = 40), # lfbt
          title = list(size = 18)
        )
      return(fig)
    }

    PlotModuleServer(
      "plot",
      func = plotly.RENDER,
      func2 = plotly.RENDER_MODAL,
      plotlib = "plotly",
      res = c(100, 95), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
