##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

CorrelationBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 800 ## full height of page
    rowH <- 340 ## full height of page

    cor_infotext <- "The <strong>Correlation Analysis Board</strong> provides statistical correlation analysis on gene level with visualisations. During the visual analysis, users can filter out some samples or collapse the samples by predetermined groups. The dark shaded area in the barplot estimates the partial correlation."

    ## ------- observe functions -----------
    shiny::observeEvent(input$data_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Expression: Correlation analysis board</strong>"),
        shiny::HTML(cor_infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    COL <- RColorBrewer::brewer.pal(12, "Paired")[seq(1, 12, 2)]
    COL <- RColorBrewer::brewer.pal(9, "Set1")[c(2, 1, 3:9)]
    COL2 <- rev(grey.colors(2))
    COL2 <- RColorBrewer::brewer.pal(2, "Paired")[1:2]
    COL2 <- COL[1:2]

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$cor_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Correlation Analysis Board</strong>"),
        shiny::HTML(cor_infotext),
        easyClose = TRUE
      ))
    })

    ## update filter choices upon change of data set
    shiny::observe({
      req(pgx$X)

      genes <- sort(pgx$genes[rownames(pgx$X), ]$gene_name)
      sel <- genes[1] ## most var gene
      sel <- names(head(sort(-rowMeans(playbase::pgx.getMetaMatrix(pgx)$fc**2)), 1))
      shiny::updateSelectizeInput(session, "cor_gene", choices = genes, selected = sel, server = TRUE)

      fam <- playbase::pgx.getFamilies(pgx, lib.dir = FILES, nmin = 10, extended = FALSE)
      fam <- sort(c("<custom>", fam))
      shiny::updateSelectInput(session, "cor_features", choices = fam)

      px <- colnames(pgx$Y)
      s1 <- grep("^[.]", px, value = TRUE, invert = TRUE)[1]
      shiny::updateSelectInput(session, "cor_group", choices = px, selected = s1)
    })


    getFilteredExpression <- shiny::reactive({
      shiny::req(pgx, input$cor_gene)
      X <- pgx$X
      gene <- rownames(X)[1]
      gene <- input$cor_gene

      ## filter genes
      ft <- input$cor_features
      if (ft == "<custom>" && input$cor_customfeatures != "") {
        genes <- toupper(pgx$genes$gene_name)
        gg1 <- strsplit(input$cor_customfeatures, split = "[, ;\n\t]")[[1]]
        if (length(gg1) == 1) gg1 <- paste0(gg1, "*")
        gg1 <- gsub("[ \n\t]", "", gg1)
        starred <- grep("[*]", gg1)
        if (length(starred) > 0) {
          gg2 <- lapply(gg1[starred], function(a) {
            genes[grep(paste0("^", sub("[*]", "", a)), genes, ignore.case = TRUE)]
          })
          gg1 <- unique(c(gg1, unlist(gg2)))
        }
        gg1 <- gg1[which(toupper(gg1) %in% toupper(pgx$genes$gene_name))]
        psel <- playbase::filterProbes(pgx$genes, c(gg1, gene))
        psel <- intersect(psel, rownames(X))
        X <- X[psel, , drop = FALSE]
      } else if (ft != "<all>" && ft %in% names(iGSETS)) {
        ft <- input$cor_features
        psel <- playbase::filterProbes(pgx$genes, c(gene, unlist(getGSETS(ft))))
        ## psel = unique(c(gene, psel))
        psel <- intersect(psel, rownames(X))
        X <- X[psel, , drop = FALSE]
      }

      X <- X + 1e-3 * matrix(rnorm(length(X)), nrow(X), ncol(X))
      X
    })

    getPartialCorrelationMatrix <- shiny::reactive({
      shiny::req(pgx, input$cor_gene)

      gene <- rownames(pgx$X)[1]
      gene <- input$cor_gene

      ## filter gene expression matrix
      X <- getFilteredExpression()

      # shiny::showNotification(paste("Computing correlation...\n"))
      NTOP <- 50
      NTOP <- as.integer(input$pcor_ntop)
      ## res <- playbase::pgx.computePartialCorrelationAroundGene(
      ##    X, gene, method=methods, nmax=NTOP, fast=FALSE)
      res <- playbase::pgx.computeGlassoAroundGene(X, gene, nmax = NTOP)
      res$meta.pcor <- res$pcor

      j <- which(rownames(res$pcor) == gene)
      P <- res$pcor
      diag(P) <- 0
      rho1 <- min(head(sort(P, decreasing = TRUE), 200))
      max1 <- round(max(P), digits = 3)
      shiny::updateSliderInput(session, "cor_graph-cor_graph_threshold", value = rho1, max = max1)
      shiny::updateSliderInput(session, "dcga_graph_threshold", value = rho1, max = max1)

      res
    })

    getPartialCorrelation <- shiny::reactive({
      res <- getPartialCorrelationMatrix()
      gene <- rownames(res$cor)[1]
      gene <- input$cor_gene
      rho <- res$cor[gene, ]
      prho <- res$pcor[gene, ]
      df <- data.frame(cor = rho, pcor = prho)
      df
    })

    getGeneCorr <- shiny::reactive({
      shiny::req(pgx)
      gene <- input$cor_gene
      if (is.null(gene)) {
        return(NULL)
      }

      ## corr always in log.scale and restricted to selected samples subset
      zx <- pgx$X
      zx <- getFilteredExpression()
      dim(zx)
      zx.genes0 <- rownames(zx)
      ## rownames(zx) <- toupper(sub(".*:","",rownames(zx)))  ## NEED RETHINK!
      zx.genes <- as.character(pgx$genes[rownames(zx), ]$gene_name)
      rownames(zx) <- toupper(zx.genes)
      xref <- list(
        "cor" = 2**zx,
        "cor.HPA" = as.matrix(TISSUE),
        "cor.ImmProt" = as.matrix(IMMPROT)
      )
      gene0 <- toupper(gene) ## uppercase mouse

      R <- playbase::pgx.getGeneCorrelation(gene0, xref = xref)
      if (is.null(R)) {
        return(NULL)
      }
      R <- R[rownames(zx), , drop = FALSE]

      zx <- zx - rowMeans(zx, na.rm = TRUE)
      sdx <- sqrt(rowMeans(zx**2))
      R <- cbind(R, cov = R[, "cor"] * sdx * sdx[gene])

      rho.genes <- rownames(zx)
      if ("hgnc_symbol" %in% colnames(pgx$genes)) {
        rho.genes <- as.character(pgx$genes[zx.genes0, ]$hgnc_symbol)
      }
      R <- R[match(rho.genes, rownames(R)), , drop = FALSE]
      rownames(R) <- zx.genes0

      R <- R[order(R[, "cor"], decreasing = TRUE), , drop = FALSE]

      R
    })

    getFullGeneCorr <- shiny::reactive({
      shiny::req(pgx)

      gene <- input$cor_gene
      if (is.null(gene)) {
        return(NULL)
      }

      ## corr always in log.scale and restricted to selected samples subset
      zx <- pgx$X
      zx.genes0 <- rownames(zx)
      ## rownames(zx) <- toupper(sub(".*:","",rownames(zx)))  ## NEED RETHINK!
      zx.genes <- as.character(pgx$genes[rownames(zx), ]$gene_name)
      rownames(zx) <- toupper(zx.genes)
      xref <- list(
        "cor" = 2**zx,
        "cor.HPA" = as.matrix(TISSUE),
        "cor.ImmProt" = as.matrix(IMMPROT)
      )
      gene0 <- toupper(gene) ## uppercase mouse

      R <- playbase::pgx.getGeneCorrelation(gene0, xref = xref)
      if (is.null(R)) {
        return(NULL)
      }
      R <- R[rownames(zx), , drop = FALSE]

      zx <- zx - rowMeans(zx, na.rm = TRUE)
      sdx <- sqrt(rowMeans(zx**2))
      R <- cbind(R, cov = R[, "cor"] * sdx * sdx[gene])

      rho.genes <- rownames(zx)
      if ("hgnc_symbol" %in% colnames(pgx$genes)) {
        rho.genes <- as.character(pgx$genes[zx.genes0, ]$hgnc_symbol)
      }
      R <- R[match(rho.genes, rownames(R)), , drop = FALSE]
      rownames(R) <- zx.genes0
      R <- R[order(R[, "cor"], decreasing = TRUE), , drop = FALSE]

      R
    })


    ## ================================================================================
    ## ======================= PLOTTING FUNCTIONS =====================================
    ## ================================================================================

    ## ----------------------------------------------------------------------
    ##  Cumulative correlation plot (stacked)
    ## ----------------------------------------------------------------------

    cum_corplot_data <- shiny::reactive({
      shiny::req(pgx)
      R <- getGeneCorr()

      ## get top correlated genes
      ## jj = head(order(rowSums(R),decreasing=FALSE),35)
      rsum <- rowSums(R, na.rm = TRUE)
      jj <- head(order(abs(rsum), decreasing = TRUE), 35)
      jj <- head(order(-abs(rsum)), 30)
      jj <- c(head(order(rsum), 15), head(order(-rsum), 15))
      jj <- jj[order(-rsum[jj])]
      head(rsum[jj])
      Rtop <- R[jj, , drop = FALSE]
      rownames(Rtop) <- sub(".*:", "", rownames(Rtop))
      offset <- min(Rtop, na.rm = TRUE) * 0.95
      offset <- 0
      klr <- grey.colors(ncol(Rtop), start = 0.3, end = 0.7)

      ## --- color test -----##
      klr <- grey.colors(ncol(Rtop), start = 0.3, end = 0.7)
      klr <- colorRampPalette(c(rgb(0.2, 0.5, 0.8, 0.8), rgb(0.2, 0.5, 0.8, 0.2)), alpha = TRUE)(ncol(Rtop))

      res <- list(Rtop = Rtop, offset = offset, klr = klr)
      return(res)
    })

    cum_corplot.RENDER <- shiny::reactive({
      res <- cum_corplot_data()
      if (is.null(res)) {
        return(NULL)
      }

      par(mar = c(6, 4, 2, 1), mgp = c(2.2, 0.8, 0))
      mar <- MARGINS1
      par(mar = mar, mgp = c(1.5, 0.5, 0))

      barplot(t(res$Rtop) - res$offset,
        col = res$klr, border = NA, ## horiz=TRUE,
        las = 3, cex.names = 0.73, ## names.arg=rep(NA,nrow(R)),
        offset = res$offset, ylab = "cumulative correlation (r)"
        ## cex.main=1.2, main="cumulative correlation\nwith other data sets"
      )
      if (!is.null(colnames(res$Rtop))) {
        legend("topright",
          legend = rev(colnames(res$Rtop)), fill = rev(res$klr),
          cex = 0.8, y.intersp = 0.8
        )
      }
    })

    cum_corplot_text <- paste0("Top cumulative positively and negatively correlated genes with the selected gene in the current dataset as well as in public datasets such as ", a_ImmProt, " and ", a_HPA, ". The correlations of genes are colored by dataset.")

    shiny::callModule(
      plotModule, "cum_corplot",
      func = cum_corplot.RENDER,
      func2 = cum_corplot.RENDER,
      info.text = cum_corplot_text,
      height = 400,
      pdf.width = 8, pdf.height = 6,
      label = "f",
      title = "Cumulative correlation",
      add.watermark = WATERMARK
    )

    ## --------------------------------------------------------------------------------
    ## Correlation GSEA
    ## --------------------------------------------------------------------------------

    getCorrelationGSEA <- shiny::reactive({
      alertDataLoaded(session, pgx)
      shiny::req(pgx)

      pgx.showSmallModal("Calculating correlation GSEA...<br>please wait")

      gene <- "CD4"
      gene <- rownames(pgx$X)[1]
      gene <- input$cor_gene

      ## single gene correlation as rank metric
      gx <- pgx$X[gene, ]
      rho <- cor(t(pgx$X), gx, use = "pairwise")[, 1]
      names(rho) <- toupper(names(rho))

      ## gmt <- GSETS[colnames(pgx$GMT)]
      gmt <- getGSETS(colnames(pgx$GMT))
      ## gmt <- GSETS  ## all???
      gsea <- fgsea::fgsea(gmt, rho, minSize = 15, maxSize = 1000)
      gsea <- gsea[order(-gsea$NES), ]
      head(gsea)

      ## done!
      shiny::removeModal()
      beepr::beep(10) ## short beep

      res <- list(gsea = gsea, rho = rho)
      return(res)
    })


    ## -----------------------------------------------------------------------------
    ## DGCA Correlation scatter plots
    ## -----------------------------------------------------------------------------

    gene1 <- "MKI67"
    gene2 <- "NCAPD2"

    dgca.scatterplot <- function(X, gene1, gene2, grp, cex = 1, key = TRUE,
                                 rho = TRUE, col = c("grey40", "red2")) {
      grp.levels <- sort(unique(setdiff(as.character(grp), NA)))
      x1 <- X[gene1, which(grp == grp.levels[1])]
      y1 <- X[gene2, which(grp == grp.levels[1])]
      x2 <- X[gene1, which(grp == grp.levels[2])]
      y2 <- X[gene2, which(grp == grp.levels[2])]

      col1 <- paste0(gplots::col2hex(col), "AA") ## add opacity
      rgb2col <- function(cc) rgb(cc[1], cc[2], cc[3], maxColorValue = 255)
      col2 <- apply(0.66 * col2rgb(col), 2, rgb2col)

      ## par(mfrow=c(2,2), mar=c(4,4,2,2), oma=c(0,0,0,0))
      par(mgp = c(1.6, 0.6, 0))
      base::plot(x1, y1,
        col = col1[1], pch = 19, cex = cex,
        xlim = range(c(x1, x2)), ylim = range(c(y1, y2)),
        ## main = 'differential correlation',
        xlab = gene1, ylab = gene2
      )

      y1 <- y1 + 1e-3 * rnorm(length(y1))
      x1 <- x1 + 1e-3 * rnorm(length(x1))
      y2 <- y2 + 1e-3 * rnorm(length(y2))
      x2 <- x2 + 1e-3 * rnorm(length(x2))

      abline(lm(y1 ~ x1), col = col2[1], lty = 1, lwd = 1.5)
      points(x2, y2, col = col1[2], pch = 19, cex = cex)
      abline(lm(y2 ~ x2), col = col2[2], lty = 1, lwd = 1.5)

      if (key) {
        legend("topleft",
          legend = grp.levels,
          fill = col, bty = "n",
          cex = 0.9, y.intersp = 0.8
        )
      }
      if (rho) {
        r1 <- round(cor(x1, y1), 3)
        r2 <- round(cor(x2, y2), 3)
        rr <- paste0("R_", grp.levels, " = ", c(r1, r2))
        legend("bottomright",
          legend = rr,
          ## fill=c('red','blue'),
          bty = "n",
          cex = 0.9, y.intersp = 0.85
        )
      }
    }

    dgca_scatter.PLOTFUN <- shiny::reactive({
      shiny::req(input$cor_gene)
      res <- dgca.output()

      ii <- dgca_table$rows_all()
      if (is.null(ii) || length(ii) == 0) {
        return(NULL)
      }
      res1 <- res[ii, , drop = FALSE]
      ## res1 <- res1[rowSums(is.na(res1))==0,]

      ph <- "ER_STATUS"
      ph <- input$cor_group
      shiny::req(ph)

      grp <- factor(pgx$samples[, ph])
      ndim <- nrow(pgx$samples)
      sel1 <- which(as.integer(grp) == 1)
      sel2 <- which(as.integer(grp) == 2)
      NTOP <- 25

      nplots <- min(NTOP, nrow(res1))
      nr <- ceiling(sqrt(nplots))
      par(
        mfrow = c(nr, nr), mar = c(3, 3.5, 0.5, 0),
        mgp = c(1.9, 0.7, 0), oma = c(0, 0, 0.5, 0.5)
      )

      cex_levels <- c(1.2, 0.8, 0.5, 0.2)
      dim_cuts <- c(0, 40, 100, 200, Inf)
      cex <- cex_levels[findInterval(ndim, dim_cuts)]

      klrpal <- c("grey30", "red2")
      klrpal <- COL2
      klrpal <- COL ## more colors
      if (nr <= 2) cex <- cex * 2

      i <- 1
      for (i in 1:nplots) {
        gene1 <- res1$Gene2[i]
        gene2 <- res1$Gene1[i]
        X1 <- pgx$X[c(gene1, gene2), ]
        dgca.scatterplot(X1, gene1, gene2,
          grp = grp, cex = cex,
          key = 0, rho = 1, col = klrpal
        )
        if (i == 1) {
          tt <- c("   ", levels(grp))
          legend("topleft",
            legend = tt,
            fill = c(NA, klrpal), inset = c(0.01, 0.01),
            border = c(NA, "black", "black"),
            cex = 0.9, box.lwd = 0, pt.lwd = 0,
            x.intersp = 0.5, y.intersp = 0.8
          )
          legend("topleft", ph,
            x.intersp = -0.2, inset = c(0.01, 0.01),
            cex = 0.9, y.intersp = 0.45, bty = "n"
          )
        }
      }
    })

    dgca_scatter.opts <- shiny::tagList(
      ## checkboxInput(ns("dgcascatter.swapaxis"),"swap axes"),
      ## selectInput(ns("dgcascatter.colorby"),"color by:", choices=NULL),
    )

    dgca_scatter.info <- "<b>DGCA scatter plots.</b> Pairwise scatter plots for the co-expression of the gene pairs in two different conditions. Differentially correlated gene pairs will show different correlation values measured by the difference in their z-score ('zScoreDiff'). The straight lines correspond to the linear regression fits. "

    shiny::callModule(
      plotModule,
      id = "dgca_scatter", label = "c",
      func = dgca_scatter.PLOTFUN,
      func2 = dgca_scatter.PLOTFUN,
      info.text = dgca_scatter.info,
      options = dgca_scatter.opts,
      title = "DGCA correlation scatter plots",
      ## pdf.width = 14, pdf.height = 4,
      height = c(fullH - 80, 760),
      width = c("auto", 900),
      res = c(80, 95),
      add.watermark = WATERMARK
    )

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    WATERMARK <- FALSE

    correlation_plot_table_corr_server(
      "cor_barplot",
      getPartialCorrelation = getPartialCorrelation,
      getGeneCorr           = getGeneCorr,
      cor_table             = cor_table,
      pgx             = pgx,
      pcor_ntop             = input$pcor_ntop,
      scrollY               = "calc(35vh - 140px)",
      watermark             = WATERMARK
    )

    correlation_plot_scattercorr_server(
      "cor_scatter",
      getFilteredExpression = getFilteredExpression,
      pgx = pgx,
      getPartialCorrelationMatrix = getPartialCorrelationMatrix,
      getGeneCorr = getGeneCorr,
      cor_gene = input$cor_gene,
      COL = COL,
      watermark = WATERMARK
    )

    correlation_plot_correlation_UMAP_server(
      "cor_umap",
      pgx       = pgx,
      cor_gene        = input$cor_gene,
      getFullGeneCorr = getFullGeneCorr,
      getGeneCorr     = getGeneCorr,
      watermark       = WATERMARK
    )

    correlation_plot_cor_graph_server(
      "cor_graph",
      cor_gene = input$cor_gene,
      getPartialCorrelationMatrix = getPartialCorrelationMatrix,
      watermark = WATERMARK
    )

    # output$test <- visNetwork::renderVisNetwork(cor_graph.VISNETWORK())
  })
} ## end of Board
