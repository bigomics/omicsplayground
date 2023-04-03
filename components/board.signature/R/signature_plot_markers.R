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
signature_plot_markers_ui <- function(id, height) {
  ns <- shiny::NS(id)
  info_text <- "<b>Markers t-SNE plot</b>. T-SNE plot for each gene, where the dot (corresponding to samples) are colored depending on the upregulation (in red) or downregulation (in blue) of that particular gene. After uploading a gene list, the <strong>Markers</strong> section produces a t-SNE plot of samples for each gene, where the samples are colored with respect to the upregulation (in red) or downregulation (in blue) of that particular gene."

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

  PlotModuleUI(ns("plot"),
    title = "Markers plot",
    options = markers.opts,
    info.text = info_text,
    download.fmt = c("png", "pdf"),
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
                                          IMMCHECK.GENES,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    calcSingleSampleValues <- function(X, y, method = c("rho", "gsva")) {
      ##
      ## Calculates single-sample enrichment values for given matrix and
      ## binarized signature vector.
      ##
      ##
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

      ## calculate GSVA
      if ("gsva" %in% method) {
        gset <- names(y)[which(y != 0)]
        gmt <- list("gmt" = gset)
        res.gsva <- GSVA::gsva(X, gmt, method = "gsva", parallel.sz = 1) ## parallel=buggy
        res.colnames <- colnames(res.gsva)
        fc <- as.vector(res.gsva[1, ])
        names(fc) <- res.colnames
        S[["gsva"]] <- fc[colnames(X)]
      }
      s.names <- names(S)
      if (length(S) > 1) {
        S1 <- do.call(cbind, S)
      } else {
        S1 <- S[[1]]
      }

      S1 <- as.matrix(S1)
      rownames(S1) <- colnames(X)
      colnames(S1) <- s.names
      return(S1)
    }

    getSingleSampleEnrichment <- shiny::reactive({
      ##
      ## Calls calcSingleSampleValues() and calculates single-sample
      ## enrichment values for complete data matrix and reduced data by
      ## group (for currentmarkers)
      ##
      ##
      if (is.null(pgx)) {
        return(NULL)
      }

      ## select samples
      X <- pgx$X
      sel <- colnames(X)
      X <- X[, sel]

      ## get the signature
      gset <- strsplit(IMMCHECK.GENES, split = " ")[[1]]
      gset <- getCurrentMarkers()
      if (is.null(gset)) {
        return(NULL)
      }

      xgene <- pgx$genes[rownames(X), "gene_name"]
      y <- 1 * (toupper(xgene) %in% toupper(gset))
      names(y) <- rownames(X)

      ## expression by group
      ## grp = pgx$samples[colnames(X),"group"]
      grp <- pgx$model.parameters$group
      groups <- unique(grp)
      gX <- sapply(groups, function(g) rowMeans(X[, which(grp == g), drop = FALSE]))
      colnames(gX) <- groups

      ## for large datasets pre-grouping is faster
      ss.bygroup <- calcSingleSampleValues(gX, y, method = c("rho", "gsva"))
      do.rho <- TRUE
      ss1 <- calcSingleSampleValues(X[, ], y, method = c("rho"))
      ss.bysample <- cbind(rho = ss1)

      res <- list(by.sample = ss.bysample, by.group = ss.bygroup)
      return(res)
    })

    markers.RENDER <- shiny::reactive({
      if (is.null(pgx)) {
        return(NULL)
      }

      markers <- pgx$families[[2]]
      markers <- COLLECTIONS[[10]]
      markers <- getCurrentMarkers()
      if (is.null(markers)) {
        return(NULL)
      }

      level <- "gene"
      xgene <- pgx$genes[rownames(pgx$X), ]$gene_name
      jj <- match(toupper(markers), toupper(xgene))
      jj <- setdiff(jj, NA)
      gx <- pgx$X[jj, , drop = FALSE]

      if (nrow(gx) == 0) {
        cat("WARNING:: Markers:: markers do not match!!\n")
        return(NULL)
      }

      ## get t-SNE positions of samples
      pos <- pgx$tsne2d[colnames(gx), ]
      gx <- gx - min(gx, na.rm = TRUE) + 0.001 ## subtract background
      grp <- pgx$model.parameters$group
      zx <- t(apply(gx, 1, function(x) tapply(x, as.character(grp), mean)))
      gx <- gx[order(-apply(zx, 1, sd)), , drop = FALSE]
      rownames(gx) <- sub(".*:", "", rownames(gx))

      ## ---------------- get GSVA values
      res <- getSingleSampleEnrichment()
      if (is.null(res)) {
        return(NULL)
      }

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
      cex2 <- ifelse(level == "gene", 1, 0.8)
      klrpal <- colorRampPalette(c("grey90", "grey60", "red3"))(16)

      nmax <- NULL
      if (input$markers_layout == "6x6") {
        nmax <- 35
        par(mfrow = c(6, 6), mar = c(0, 0.2, 0.5, 0.2), oma = c(2, 1, 2, 1) * 0.8)
      }
      if (input$markers_layout == "4x4") {
        nmax <- 15
        par(mfrow = c(4, 4), mar = c(0, 0.2, 0.5, 0.2), oma = c(2, 1, 2, 1) * 0.8)
      }

      top.gx <- head(gx, nmax)
      if (input$markers_sortby == "name") {
        top.gx <- top.gx[order(rownames(top.gx)), , drop = FALSE]
      }
      if (input$markers_sortby == "probability") {
        top.gx <- top.gx[order(-rowMeans(top.gx)), , drop = FALSE]
      }
      if (input$markers_sortby == "correlation") {
        rho <- cor(t(top.gx), fc1)[, 1]
        top.gx <- top.gx[order(-rho), , drop = FALSE]
      }

      for (i in 0:min(nmax, nrow(top.gx))) {
        jj <- 1:ncol(top.gx)
        if (i == 0) {
          klr1 <- BLUERED(16)[8 + round(7 * fc1)]
          tt <- "INPUT SIGNATURE"
          jj <- order(abs(fc1))
        } else {
          colvar <- pmax(top.gx[i, ], 0)
          colvar <- 1 + round(15 * (colvar / (0.7 * max(colvar) + 0.3 * max(top.gx))))
          klr1 <- klrpal[colvar]
          gene <- substring(sub(".*:", "", rownames(top.gx)[i]), 1, 80)
          tt <- breakstring(gene, n = 20, force = TRUE)
          jj <- order(abs(top.gx[i, ]))
        }
        klr1 <- paste0(gplots::col2hex(klr1), "99")

        base::plot(pos[jj, ],
          pch = 19, cex = cex1, col = klr1[jj],
          xlim = 1.2 * range(pos[, 1]), ylim = 1.2 * range(pos[, 2]),
          fg = gray(ifelse(i == 0, 0.1, 0.8)), bty = "o",
          xaxt = "n", yaxt = "n", xlab = "tSNE1", ylab = "tSNE2"
        )
        legend("topleft", tt,
          cex = cex2, col = "grey30", text.font = ifelse(i == 0, 2, 1),
          inset = c(-0.1, -0.05), bty = "n"
        )
      }
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = markers.RENDER,
      res = c(100, 95), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
