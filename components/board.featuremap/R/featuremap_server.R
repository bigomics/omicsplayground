##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

FeatureMapBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 800 ## full height of page
    rowH1 <- 220 ## row 1 height
    rowH2 <- 460 ## row 2 height

    infotext <- "Visually explore and compare expression signatures on UMAP plots. Feature-level clustering is based on pairwise co-expression between genes (or genesets). This is in contrast to sample-level clustering which clusters samples by similarity of their expression profile. Feature-level clustering allows one to detect gene modules, explore gene neighbourhoods, and identify potential drivers, to study the relationships between features.
<br><br>The tabs present Gene and Geneset UMAP dimensionality reduction plots and are computed for gene and geneset features, respectively. The clustering of features is computed using UMAP from either the normalized log-expression matrix (logCPM) or the log-foldchange matrix (logFC), with the covariance as distance metric. The UMAP from the logCPM is the default, but in cases of strong batch/tissue effects the UMAP from the logFC matrix is a better choice. We prefer the covariance distance metric instead of the correlation because it takes the size of the foldchange into account. Doing so, genes that are close together in corners in the outer rim are those with high pairwise covariance, i.e. have high correlation and high FC.
<br><br>The maps can be colored according to the foldchange signature of the group contrasts (i.e. comparisons), or colored by the average relative log-expression according to some phenotype condition. Multiple signatures can then be easily compared by visually inspection of the colors.
"


    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Feature Map Analysis</strong>"),
        shiny::HTML(infotext),
        easyClose = TRUE
      ))
    })

    shiny::observe({
      shiny::req(pgx)

      families <- names(FAMILIES)
      shiny::updateSelectInput(session, "filter_genes",
        choices = families,
        selected = "<all>"
      )

      gsetcats <- sort(unique(gsub(":.*", "", rownames(pgx$gsetX))))
      shiny::updateSelectInput(session, "filter_gsets",
        choices = gsetcats,
        selected = "H"
      )

      cvar <- pgx.getCategoricalPhenotypes(pgx$samples, max.ncat = 99)
      cvar0 <- grep("^[.]", cvar, invert = TRUE, value = TRUE)[1]
      shiny::updateSelectInput(session, "sigvar",
        choices = cvar,
        selected = cvar0
      )
    })

    ## ================================================================================
    ## ============================= FUNCTIONS ========================================
    ## ================================================================================

    ## hilight=hilight2=NULL;source="";plotlib='base';cex=0.9
    plotUMAP <- function(pos, var, hilight = NULL, nlabel = 20, title = "",
                         zlim = NULL, cex = 0.9, source = "", plotlib = "base") {
      if (!is.null(hilight)) {
        hilight <- intersect(hilight, rownames(pos))
        hilight <- intersect(hilight, names(var))
        hilight <- hilight[order(-var[hilight])]
        if (min(var, na.rm = TRUE) < 0) {
          hilight2 <- c(head(hilight, nlabel / 2), tail(hilight, nlabel / 2))
          hilight2 <- unique(hilight2)
        } else {
          hilight2 <- head(hilight, nlabel)
        }
      }
      if (length(hilight) > 0.33 * length(var)) hilight <- hilight2

      cexlab <- ifelse(length(hilight2) <= 20, 1, 0.85)
      cexlab <- ifelse(length(hilight2) <= 8, 1.15, cexlab)
      opacity <- ifelse(length(hilight2) > 0, 0.4, 0.90)
      ## cex = 0.9
      ## opacity = ifelse(length(hilight)>0, 0.15, 1)
      if (plotlib == "plotly") opacity <- sqrt(opacity) ## less opacity..

      p <- pgx.scatterPlotXY(
        pos,
        var = var,
        plotlib = plotlib,
        softmax = TRUE,
        cex.lab = 1.3 * cexlab,
        opacity = opacity,
        cex = cex,
        zsym = (min(var, na.rm = TRUE) < 0),
        zlim = zlim,
        hilight.cex = cex,
        ## hilight.col = 'red',
        hilight.col = NULL,
        hilight.lwd = 0.8,
        hilight = hilight,
        hilight2 = hilight2,
        title = title,
        ## legend.pos = 'bottomright',
        source = source,
        ## key = rownames(pos)
      )

      p
    }

    plotFeaturesPanel <- function(pos, F, ntop, nr, nc, sel, progress) {
      par(mar = c(1.6, 1.5, 0.5, 0), oma = c(1, 1, 0, 0) * 2)
      par(mar = c(1.1, 1.0, 0.5, 0), oma = c(1, 1, 0, 0) * 2)
      par(mgp = c(1.35, 0.5, 0), las = 0, cex.axis = 0.85, cex.lab = 0.9, xpd = TRUE)
      ## ntop <- 10
      cex <- ifelse(nc > 3, 0.5, 0.7)
      jj <- 1:nrow(F)
      if (ncol(F) > 4 && nrow(F) > 8000) jj <- sample(1:nrow(F), 8000) ## subsample for speed
      if (ncol(F) > 9 && nrow(F) > 4000) jj <- sample(1:nrow(F), 4000) ## subsample for speed
      if (ncol(F) > 16 && nrow(F) > 2000) jj <- sample(1:nrow(F), 2000) ## subsample for speed

      ## global zlim
      qq <- quantile(F, probs = c(0.05, 0.95), na.rm = TRUE)
      qq <- quantile(F, probs = c(0.01, 0.99), na.rm = TRUE)
      qq <- quantile(F, probs = c(0.002, 0.998), na.rm = TRUE)
      c(min(qq, na.rm = TRUE), max(qq, na.rm = TRUE))
      qq
      zlim <- qq
      ## zlim = NULL

      i <- 1
      for (i in 1:ncol(F)) {
        if (!interactive()) progress$inc(1 / ncol(F))
        var <- F[, i]
        var <- var[match(rownames(pos), names(var))]
        zsym <- ifelse(min(var, na.rm = TRUE) >= 0, FALSE, TRUE)
        hmarks <- NULL
        if (!is.null(sel)) {
          hmarks <- intersect(sel, names(var))
          hmarks <- head(hmarks[order(var[hmarks])], ntop)
        }
        opacity <- ifelse(is.null(hmarks), 0.9, 0.4)
        xlab <- "UMAP-x"
        ylab <- "UMAP-y"
        xaxs <- yaxs <- TRUE
        if (i %% nc != 1) {
          ylab <- ""
          yaxs <- FALSE
        }
        if ((i - 1) %/% nc < (nr - 1)) {
          xlab <- ""
          xaxs <- FALSE
        }

        pgx.scatterPlotXY.BASE(
          pos[jj, ],
          var = var[jj],
          zsym = zsym, zlim = zlim, set.par = FALSE, softmax = 1,
          cex = cex, cex.legend = 0.9, cex.lab = 1.2, bty = "n",
          col = "grey70", dlim = c(0.05, 0.05),
          hilight = hmarks, hilight2 = NULL,
          hilight.col = NULL, opacity = opacity,
          ## xlab = xlab, ylab = ylab,
          xlab = "", ylab = "",
          xaxs = xaxs, yaxs = yaxs,
          hilight.lwd = 0.5, hilight.cex = 1.3
        )

        cex1 <- ifelse(ncol(F) <= 16, 1.2, 1)
        title(colnames(F)[i], cex.main = cex1, line = -0.75)
        mtext(xlab, 1, line = 1.5, cex = 0.6)
        mtext(ylab, 2, line = 1.6, cex = 0.6)
      }
    }

    getGeneUMAP_FC <- shiny::reactive({
      ## buffered reactive
      shiny::withProgress(
        {
          F <- pgx.getMetaMatrix(pgx, level = "gene")$fc
          F <- scale(F, center = FALSE)
          pos <- pgx.clusterBigMatrix(t(F), methods = "umap", dims = 2)[[1]]
          pos <- pos.compact(pos)
        },
        message = "computing foldchange UMAP",
        value = 0.5
      )
      pos
    })

    getGeneUMAP <- shiny::reactive({
      if (input$umap_type == "logFC") {
        message("[getGeneUMAP] computing foldchange UMAP")
        pos <- getGeneUMAP_FC()
      } else {
        pos <- pgx$cluster.genes$pos[["umap2d"]]
      }
      pos
    })

    getGsetUMAP_FC <- shiny::reactive({
      ## buffered reactive
      shiny::withProgress(
        {
          F <- pgx.getMetaMatrix(pgx, level = "geneset")$fc
          F <- scale(F, center = FALSE)
          pos <- pgx.clusterBigMatrix(t(F), methods = "umap", dims = 2)[[1]]
          pos <- pos.compact(pos)
        },
        message = "computing foldchange UMAP (genesets)",
        value = 0.5
      )
      pos
    })

    getGsetUMAP <- shiny::reactive({
      if (input$umap_type == "logFC") {
        message("[getGsetUMAP] computing foldchange UMAP (genesets)")
        pos <- getGsetUMAP_FC()
      } else {
        pos <- pgx$cluster.gsets$pos[["umap2d"]]
      }
      pos
    })

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    WATERMARK <- FALSE

    # Gene Map

    featuremap_plot_gene_map_server(
      "gene_map",
      pgx    = pgx,
      getGeneUMAP  = getGeneUMAP,
      plotUMAP     = plotUMAP,
      sigvar       = shiny::reactive(input$sigvar),
      filter_genes = shiny::reactive(input$filter_genes),
      watermark    = WATERMARK
    )

    # Gene Signatures

    featuremap_plot_gene_sig_server(
      "gene_sig",
      pgx         = pgx,
      getGeneUMAP       = getGeneUMAP,
      sigvar            = shiny::reactive(input$sigvar),
      plotFeaturesPanel = plotFeaturesPanel,
      watermark         = WATERMARK
    )

    # Geneset map

    featuremap_plot_table_geneset_map_server(
      "gsetUMAP",
      pgx = pgx,
      getGsetUMAP = getGsetUMAP,
      plotUMAP = plotUMAP,
      filter_gsets = shiny::reactive(input$filter_gsets),
      sigvar = shiny::reactive(input$sigvar),
      watermark = WATERMARK
    )

    # Geneset signatures

    featuremap_plot_gset_sig_server(
      "gsetSigPlots",
      pgx         = pgx,
      getGsetUMAP       = getGsetUMAP,
      sigvar            = shiny::reactive(input$sigvar),
      plotFeaturesPanel = plotFeaturesPanel,
      watermark         = WATERMARK
    )
  }) ## end of serverModule
} ## end of Board

## ========================================================================
## ========================= END OF FILE ==================================
## ========================================================================
