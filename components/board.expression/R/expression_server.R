##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

ExpressionBoard <- function(id, inputData) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 780
    rowH <- 390 ## row height of panels
    imgH <- 340 ## height of images
    tabV <- "70vh" ## height of tables
    tabH <- 320 ## row height of panels

    gx_infotext <- "The <strong>Differential Expression Analysis</strong> module compares expression between two conditions (i.e. tumor versus control),
     which is one of the fundamental analysis in the transcriptomics data analytics workflow. For each comparison of two conditions (also called \'contrast\'),
     the analysis identifies which genes are significantly downregulated or overexpressed in one of the groups.<br><br>
     The <strong>Plots</strong> panel shows volcano and MA plots for the chosen contrast. It also shows the so-called \'signature\',
     i.e. the top downregulated and overexpressed genes, for that contrast. The <strong>Top genes</strong> panel shows the average expression
     plots across the samples for top differentially expressed genes within the selected comparison. A very useful feature of the platform is
     that it can display volcano plots for all comparisons simultaneously under the <strong>Volcano (all)</strong> panel.
     This provides users an overview of the statistics of all comparisons. The <strong>Table</strong> panel on the bottom shows the results
     of the statistical tests. The <strong>Foldchange (all)</strong> panel reports the gene fold changes for all contrasts.
     <br><br>EXPERT MODE ONLY: To compare the different statistical methods, the <strong>Volcano (methods)</strong> panel shows volcano plots of all methods.
     The <strong>FDR table</strong> panel reports the number of significant genes at different FDR thresholds for all contrasts.<br><br><br><br>
     <center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=3'
     frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>"

    GX.DEFAULTTEST <- "trend.limma"
    GX.DEFAULTTEST <- c("trend.limma", "edger.qlf", "deseq2.wald", "edger.lrt")

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$gx_info,
      {
        shiny::showModal(shiny::modalDialog(
          title = shiny::HTML("<strong>Differential Expression Analysis Board</strong>"),
          shiny::HTML(gx_infotext),
          easyClose = TRUE, size = "l"
        ))
      },
      ignoreInit = TRUE
    )

    ## update choices upon change of data set
    shiny::observe({
      ngs <- inputData()
      shiny::req(ngs)

      contr <- colnames(ngs$model.parameters$contr.matrix)
      shiny::updateSelectInput(session, "gx_contrast", choices = sort(contr))
      fam <- pgx.getFamilies(ngs, nmin = 10, extended = FALSE)
      shiny::updateSelectInput(session, "gx_features", choices = fam)

      ## available statistical methods
      gx.methods <- colnames(ngs$gx.meta$meta[[1]]$fc) ## available
      sel1 <- c(intersect(GX.DEFAULTTEST, gx.methods), gx.methods)
      sel1 <- head(unique(sel1), 3) ## maximum three!!

      shiny::updateCheckboxGroupInput(session, "gx_statmethod",
        choices = sort(gx.methods),
        selected = sel1
      )

      shiny::updateCheckboxInput(session, "gx_ungroup", value = (ncol(ngs$X) <= 8))
    })

    ## ================================================================================
    ## ========================= REACTIVE FUNCTIONS ===================================
    ## ================================================================================

    selected_gxmethods <- shiny::reactive({
      ngs <- inputData()
      req(ngs)
      gx.methods0 <- colnames(ngs$gx.meta$meta[[1]]$fc)
      test <- input$gx_statmethod
      test <- intersect(test, gx.methods0)
      test
    })

    ## ================================================================================
    ## ========================= FUNCTIONS ============================================
    ## ================================================================================

    comparison <- 1
    testmethods <- c("trend.limma")
    add.pq <- 0
    getDEGtable <- function(ngs, testmethods, comparison, add.pq,
                            lfc, fdr) {
      ## ngs = inputData()
      ## if(is.null(ngs)) return(NULL)
      shiny::req(ngs)

      if (is.null(testmethods)) {
        return(NULL)
      }
      if (is.null(comparison)) {
        return(NULL)
      }
      if (length(testmethods) == 0 || testmethods == "") {
        return(NULL)
      }
      if (length(comparison) == 0 || comparison == "") {
        return(NULL)
      }

      message("[getDEGtable] called")

      ## build meta table
      mx <- ngs$gx.meta$meta[[comparison]]
      if (is.null(mx)) {
        return(NULL)
      }
      mm <- colnames(unclass(mx$p))
      testmethods <- intersect(mm, testmethods)

      mx.p <- unclass(mx$p[, testmethods, drop = FALSE]) ## get rid of AsIs
      mx.q <- unclass(mx$q[, testmethods, drop = FALSE])
      mx.fc <- unclass(mx$fc[, testmethods, drop = FALSE])
      ## mx$score = mx$fc * (-log10(1e-100+mx$q) )
      rownames(mx.p) <- rownames(mx)
      rownames(mx.q) <- rownames(mx)
      rownames(mx.fc) <- rownames(mx)

      mx.fc[is.infinite(mx.fc) | is.nan(mx.fc)] <- NA
      mx.p[is.infinite(mx.p) | is.nan(mx.p)] <- NA
      mx.q[is.infinite(mx.q) | is.nan(mx.q)] <- NA

      ## !!!!!!!!!!!!!!!!!!!!!!!! NEED RETHINK !!!!!!!!!!!!!!!!!!!!!!!!
      ## must recompute meta parameters (maxQ method)
      mx$meta.p <- apply(mx.p, 1, max, na.rm = TRUE)
      mx$meta.q <- apply(mx.q, 1, max, na.rm = TRUE)
      mx$meta.fx <- rowMeans(mx.fc, na.rm = TRUE)
      ## mx$meta.score = rowMeans(mx$score,na.rm=TRUE)

      stars.fdr <- fdr
      ## stars.fdr = 0.05  ## fixed otherwisetable will always have three stars..
      ## star.symbols = sapply(1:20,function(i) paste(rep("\u2605",i),collapse=""))
      ## is.sig <- (mx.q <= stars.fdr & abs(mx.fc) >= lfc)
      is.sig <- 1 * (mx.q <= stars.fdr) * (abs(mx$meta.fx) >= lfc)
      ## stars = c("",star.symbols)[ 1 + rowSums(is.sig, na.rm=TRUE)]
      stars <- sapply(rowSums(is.sig, na.rm = TRUE), star.symbols)

      ## recalculate group averages???
      y0 <- ngs$model.parameters$exp.matrix[, comparison]
      names(y0) <- rownames(ngs$model.parameters$exp.matrix)
      AveExpr1 <- rowMeans(ngs$X[rownames(mx), names(which(y0 > 0)), drop = FALSE])
      AveExpr0 <- rowMeans(ngs$X[rownames(mx), names(which(y0 < 0)), drop = FALSE])

      ## logFC <- unclass(ngs$gx.meta$meta[[comparison]][,"fc"])[,"trend.limma"]
      ## logFC <- ngs$gx.meta$meta[[comparison]][,"meta.fx"]
      logFC <- mx$meta.fx
      ## logFC <- (AveExpr1 - AveExpr0)  ## override ??? yes: see "contrast in R" Rose Maier 2015...
      ## [hack] adjust averages to match logFC...
      mean0 <- (AveExpr0 + AveExpr1) / 2
      AveExpr1 <- mean0 + logFC / 2
      AveExpr0 <- mean0 - logFC / 2

      message("[getDEGtable] creating results table")

      ## gene.annot = mx[,grep("^gene|^chr",colnames(mx)),drop=FALSE]
      aa <- intersect(c("gene_name", "gene_title", "chr"), colnames(ngs$genes))
      gene.annot <- ngs$genes[rownames(mx), aa]
      gene.annot$chr <- sub("_.*", "", gene.annot$chr) ## strip any alt postfix
      res <- data.frame(gene.annot,
        logFC = logFC,
        stars = stars, meta.q = mx$meta.q,
        AveExpr0, AveExpr1, check.names = FALSE
      )
      rownames(res) <- rownames(mx)

      if (add.pq) {
        message("[getDEGtable] adding PQ table")
        ## add extra columns
        ## res <- cbind( res, q=mx$q, p=mx$p)
        colnames(mx.q) <- paste0("q.", colnames(mx.q))
        res <- cbind(res, mx.q[rownames(mx), , drop = FALSE])
      }
      return(res)
    }

    fullDiffExprTable <- shiny::reactive({
      ## return the full DE table
      ngs <- inputData()
      if (is.null(ngs)) {
        return(NULL)
      }
      comp <- 1
      tests <- "trend.limma"
      fdr <- 1
      lfc <- 0
      comp <- input$gx_contrast
      tests <- input$gx_statmethod
      fdr <- as.numeric(input$gx_fdr)
      lfc <- as.numeric(input$gx_lfc)

      if (is.null(comp)) {
        return(NULL)
      }
      if (is.null(tests)) {
        return(NULL)
      }
      res <- getDEGtable(ngs,
        testmethods = tests, comparison = comp,
        add.pq = TRUE, lfc = lfc, fdr = fdr
      )

      ## Filter on features/genesets
      psel <- rownames(res)
      gx_features <- 1
      gx_features <- input$gx_features
      if (gx_features != "<all>") {
        ## gset <- GSETS[[gx_features]]
        gset <- unlist(getGSETS(gx_features))
        psel <- filterProbes(ngs$genes, gset)
      }
      res <- res[which(rownames(res) %in% psel), , drop = FALSE]
      dim(res)

      res <- res[order(-abs(res$logFC)), , drop = FALSE]
      return(res)
    })

    filteredDiffExprTable <- shiny::reactive({
      ##
      ## DE table filtered by FDR and gene family
      ##
      ##
      ngs <- inputData()
      ## if(is.null(ngs)) return(NULL)
      shiny::req(ngs, input$gx_features, input$gx_fdr, input$gx_lfc)

      comp <- 1
      test <- "trend.limma"
      comp <- input$gx_contrast
      tests <- input$gx_statmethod
      fdr <- as.numeric(input$gx_fdr)
      lfc <- as.numeric(input$gx_lfc)

      ## res = getDEGtable(ngs, testmethods="trend.limma", comparison=1,add.pq=FALSE)
      ## res = getDEGtable(ngs, testmethods=tests, comparison=comp,
      ## add.pq=TRUE, lfc=lfc, fdr=fdr, filter.sig=FALSE)
      res <- fullDiffExprTable()
      if (is.null(res) || nrow(res) == 0) {
        return(NULL)
      }

      fx.col <- grep("mean.diff|logfc|foldchange|meta.fx", colnames(res), ignore.case = TRUE)[1]
      res <- res[order(-abs(res[, fx.col])), ]

      ## just show significant genes
      if (!is.null(input$gx_showall) && !input$gx_showall) {
        n <- length(input$gx_statmethod)
        sel <- which(res$stars == star.symbols(n))
        res <- res[sel, , drop = FALSE]
      }

      ## just show top 10
      if (length(input$gx_top10) && input$gx_top10) {
        fx <- as.numeric(res[, fx.col])
        names(fx) <- rownames(res)
        pp <- unique(c(
          head(names(sort(-fx[which(fx > 0)])), 10),
          head(names(sort(fx[which(fx < 0)])), 10)
        ))
        res <- res[pp, , drop = FALSE]
        res <- res[order(-res[, fx.col]), , drop = FALSE]
      }

      if (nrow(res) == 0) {
        shiny::validate(shiny::need(nrow(res) > 0, "warning. no genes passed current filters."))
        return(NULL)
      }

      ## limit number of rows???
      ## res <- head(res, 1000)
      return(res)
    })
    ## %>%
    ## bindCache(input$gx_features, input$gx_fdr, input$gx_lfc)

    # Plotting ###

    # tab differential expression > Plots ####

    expression_plot_volcano_server(
      id = "plots_volcano",
      pgx_fdr = reactive(input$gx_fdr),
      pgx_contrast = reactive(input$gx_contrast),
      pgx_lfc = reactive(input$gx_lfc),
      pgx_features = reactive(input$gx_features),
      res = fullDiffExprTable,
      sel1 = genetable$rows_selected,
      df1 = filteredDiffExprTable,
      sel2 = gsettable$rows_selected,
      df2 = gx_related_genesets,
      fam.genes = res$gene_name
    )

    expression_plot_maplot_server(
      id = "plots_maplot",
      inputData = inputData,
      gx_fdr = reactive(input$gx_fdr),
      gx_contrast = reactive(input$gx_contrast),
      gx_lfc = reactive(input$gx_lfc),
      gx_features = reactive(input$gx_features),
      res = fullDiffExprTable,
      sel1 = genetable$rows_selected,
      df1 = filteredDiffExprTable,
      sel2 = gsettable$rows_selected,
      df2 = gx_related_genesets,
      fam.genes = res$gene_name,
      watermark = FALSE
    )

    expression_plot_boxplot_server(
      id = "plots_barplot",
      inputData = inputData,
      sel = genetable$rows_selected,
      res = filteredDiffExprTable,
      watermark = FALSE
    )

    expression_plot_topfoldchange_server(
      id = "plots_topfoldchange",
      inputData = inputData,
      sel = genetable$rows_selected,
      res = filteredDiffExprTable,
      watermark = FALSE
    )

    # tab differential expression > Top genes ####

    getAllContrasts <- shiny::reactive({
      ngs <- inputData()
      if (is.null(ngs)) {
        return(NULL)
      }
      comp <- names(ngs$gx.meta$meta)
      if (length(comp) == 0) {
        return(NULL)
      }
      ## if(is.null(input$gx_features)) return(NULL)

      ## fdr=1;lfc=0
      ## fdr = as.numeric(input$gx_fdr)
      ## lfc = as.numeric(input$gx_lfc)
      tests <- colnames(ngs$gx.meta$meta[[1]]$p)
      tests <- input$gx_statmethod
      if (is.null(tests)) {
        return(NULL)
      }

      ## comp <- head(comp,75)  ## maximum 75!!!
      i <- 1
      F <- list()
      Q <- list()
      shiny::withProgress(message = "computing contrasts ...", value = 0, {
        for (i in 1:length(comp)) {
          res <- getDEGtable(ngs,
            testmethods = tests, comparison = comp[i],
            add.pq = FALSE, lfc = 0, fdr = 1
          )
          fc.gene <- res[, grep("^gene$|^gene_name$", colnames(res))]
          ## pv.col = grep("p.val|pval|meta.p",colnames(res),ignore.case=TRUE)[1]
          qv.col <- grep("qval|adj.p|padj|fdr|meta.q", colnames(res), ignore.case = TRUE)[1]
          fx.col <- grep("mean.diff|logfc|foldchange|meta.fx", colnames(res), ignore.case = TRUE)[1]
          qval <- res[, qv.col]
          fx <- res[, fx.col]
          names(qval) <- names(fx) <- fc.gene
          F[[i]] <- fx
          Q[[i]] <- qval

          if (!interactive()) shiny::incProgress(1 / length(comp))
        }
      })
      names(Q) <- names(F) <- comp

      ct <- list(Q = Q, F = F)
      ct
    })

    expression_plot_topgenes_server(
      id = "topgenes",
      inputData = inputData,
      res = filteredDiffExprTable,
      ii = genetable$rows_current
    )

    ## ================================================================================
    ## Volcano (all contrasts)
    ## ================================================================================


    volcanoAll.RENDER <- shiny::reactive({
      ## volcanoAll.RENDER <- shiny::reactive({

      ngs <- inputData()
      if (is.null(ngs)) {
        return(NULL)
      }
      ct <- getAllContrasts()
      F <- ct$F
      Q <- ct$Q

      ## comp = names(ngs$gx.meta$meta)
      comp <- names(F)
      if (length(comp) == 0) {
        return(NULL)
      }
      if (is.null(input$gx_features)) {
        return(NULL)
      }

      fdr <- 1
      lfc <- 0
      fdr <- as.numeric(input$gx_fdr)
      lfc <- as.numeric(input$gx_lfc)

      sel.genes <- rownames(ngs$X)
      if (input$gx_features != "<all>") {
        gset <- getGSETS(input$gx_features)
        sel.genes <- unique(unlist(gset))
      }

      ## -------------------------------------------------
      ## plot layout
      ## -------------------------------------------------
      ng <- length(comp)
      nn <- c(2, max(ceiling(ng / 2), 5))
      ## if(ng>12) nn = c(3,8)
      par(mfrow = nn, mar = c(1, 1, 1, 1) * 0.2, mgp = c(2.6, 1, 0), oma = c(1, 1, 0, 0) * 2)
      nr <- 2
      nc <- ceiling(sqrt(ng))
      if (ng > 24) {
        nc <- max(ceiling(ng / 3), 6)
        nr <- 3
      } else if (TRUE && ng <= 4) {
        nc <- 4
        nr <- 1
      } else {
        nc <- max(ceiling(ng / 2), 6)
        nr <- 2
      }
      nr
      nc
      par(mfrow = c(nr, nc))

      ymax <- 15
      nlq <- -log10(1e-99 + unlist(Q))
      ymax <- max(1.3, 1.2 * quantile(nlq, probs = 0.999, na.rm = TRUE)[1]) ## y-axis
      xmax <- max(1, 1.2 * quantile(abs(unlist(F)), probs = 0.999, na.rm = TRUE)[1]) ## x-axis

      shiny::withProgress(message = "rendering volcano plots ...", value = 0, {
        plt <- list()
        i <- 1
        for (i in 1:length(comp)) {
          qval <- Q[[i]]
          fx <- F[[i]]
          fc.gene <- names(qval)
          is.sig <- (qval <= fdr & abs(fx) >= lfc)
          sig.genes <- fc.gene[which(is.sig)]
          genes1 <- sig.genes[which(toupper(sig.genes) %in% toupper(sel.genes))]
          genes2 <- head(genes1[order(-abs(fx[genes1]) * (-log10(qval[genes1])))], 10)
          xy <- data.frame(x = fx, y = -log10(qval))
          is.sig2 <- factor(is.sig, levels = c(FALSE, TRUE))

          plt[[i]] <- pgx.scatterPlotXY.GGPLOT(
            xy,
            title = comp[i], cex.title = 0.85,
            var = is.sig2, type = "factor",
            col = c("#bbbbbb", "#1e60bb"),
            legend.pos = "none", ## plotlib="ggplot",
            hilight = NULL, hilight2 = genes2,
            xlim = xmax * c(-1, 1), ylim = c(0, ymax),
            xlab = "difference  (log2FC)",
            ylab = "significance  (-log10q)",
            hilight.lwd = 0, hilight.col = "#1e60bb", hilight.cex = 1.5,
            cex = 0.45, cex.lab = 0.62
          )
          ## ggplot2::theme(legend.position='none')
          ## ggplot2::theme_bw(base_size=11)

          if (!interactive()) shiny::incProgress(1 / length(comp))
        }
      }) ## progress


      ## patchwork::wrap_plots(plt, nrow=nr, ncol=nc) &
      ##    ggplot2::theme_bw(base_size=11) &
      ##    ggplot2::theme(legend.position='none')


      gridExtra::grid.arrange(grobs = plt, nrow = nr, ncol = nc)
    })

    volcanoAll_text <- "Under the <strong>Volcano (all)</strong> tab, the platform simultaneously displays multiple volcano plots for genes across all contrasts. This provides users an overview of the statistics for all comparisons. By comparing multiple volcano plots, the user can immediately see which comparison is statistically weak or strong."

    shiny::callModule(plotModule,
      id = "volcanoAll",
      func = volcanoAll.RENDER,
      func2 = volcanoAll.RENDER,
      info.text = volcanoAll_text,
      pdf.width = 16, pdf.height = 5,
      height = c(imgH, 500), width = c("auto", 1600),
      res = c(70, 90),
      title = "Volcano plots for all contrasts",
      add.watermark = WATERMARK
    )


    ## ================================================================================
    ## Volcano (all2 contrasts)
    ## ================================================================================

    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ## PLOTS SEEMS NOT TO REFRESH/DRAW CORRECTLY. Maybe viz.Contrast is isolated????
    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ## volcanoAll2.RENDER <- shiny::reactive({
    volcanoAll2.RENDER <- shiny::reactive({
      ngs <- inputData()
      if (is.null(ngs)) {
        return(NULL)
      }


      fdr <- 1
      lfc <- 0
      fdr <- as.numeric(input$gx_fdr)
      lfc <- as.numeric(input$gx_lfc)

      sel.genes <- rownames(ngs$X)
      if (input$gx_features != "<all>") {
        ## gset <- GSETS[input$gx_features]
        gset <- getGSETS(input$gx_features)
        sel.genes <- unique(unlist(gset))
      }

      ## -------------------------------------------------
      ## plot layout
      ## -------------------------------------------------
      comp <- names(ngs$gx.meta$meta)
      ng <- length(comp)
      nn <- c(2, max(ceiling(ng / 2), 5))
      ## if(ng>12) nn = c(3,8)
      par(mfrow = nn, mar = c(1, 1, 1, 1) * 0.2, mgp = c(2.6, 1, 0), oma = c(1, 1, 0, 0) * 2)
      nr <- 2
      nc <- ceiling(sqrt(ng))
      if (ng > 24) {
        nc <- max(ceiling(ng / 3), 6)
        nr <- 3
      } else if (TRUE && ng <= 4) {
        nc <- 4
        nr <- 1
      } else {
        nc <- max(ceiling(ng / 2), 6)
        nr <- 2
      }
      nr
      nc
      ## par(mfrow=c(nr,nc))

      tests <- "meta"
      methods <- NULL
      methods <- selected_gxmethods()
      plist <- viz.Contrasts(
        pgx = ngs, ## pgxRT=inputData,
        methods = methods, type = "volcano", fixed.axis = TRUE,
        psig = fdr, fc = lfc, ntop = 10, cex = 0.5, cex.lab = 0.7,
        plots.only = TRUE, title = NULL, subtitle = NULL, caption = NULL
      )

      fig <- viz.showFigure(plist) + plot_layout(nrow = nr, ncol = nc) &
        ggplot2::theme_bw(base_size = 11) &
        ## ggplot2::theme_bw(base_size=16) &
        ggplot2::theme(legend.position = "none")

      fig
    })

    volcanoAll2_text <- "Under the <strong>Volcano (all)</strong> tab, the platform simultaneously displays multiple volcano plots for genes across all contrasts. This provides users an overview of the statistics for all comparisons. By comparing multiple volcano plots, the user can immediately see which comparison is statistically weak or strong."

    volcanoAll2_caption <- "<b>Volcano plot for all contrasts.</b> Simultaneous visualisation of volcano plots of genes for all contrasts. Experimental contrasts with better statistical significance will show volcano plots with 'higher' wings."

    shiny::callModule(
      plotModule,
      id = "volcanoAll2",
      func = volcanoAll2.RENDER,
      func2 = volcanoAll2.RENDER,
      ## plotlib = 'ggplot',
      info.text = volcanoAll2_text,
      ## caption = volcanoAll_caption,
      pdf.width = 16, pdf.height = 5,
      ## height = imgH, res=75,
      height = c(imgH, 500), width = c("auto", 1600),
      res = c(75, 95),
      title = "Volcano plots for all contrasts",
      add.watermark = WATERMARK
    )


    output$volcanoAll2_UI <- shiny::renderUI({
      shiny::fillCol(
        ## id = ns("topgenes"),
        height = rowH,
        flex = c(1, NA, NA), ## height = 370,
        plotWidget(ns("volcanoAll2")),
        shiny::br(),
        shiny::div(shiny::HTML(volcanoAll_caption), class = "caption")
      )
    })

    ## ================================================================================
    ## Volcano (all methods)
    ## ================================================================================

    volcanoMethods.RENDER <- shiny::reactive({
      comp <- input$gx_contrast
      if (is.null(comp)) {
        return(NULL)
      }
      ngs <- inputData()
      shiny::req(ngs)
      if (is.null(input$gx_features)) {
        return(NULL)
      }

      fdr <- 1
      lfc <- 1
      comp <- names(ngs$gx.meta$meta)[1]
      fdr <- as.numeric(input$gx_fdr)
      lfc <- as.numeric(input$gx_lfc)
      genes <- NULL

      gset <- getGSETS(input$gx_features)
      sel.genes <- unique(unlist(gset))

      ## meta tables
      mx <- ngs$gx.meta$meta[[comp]]
      fc <- unclass(mx$fc)
      ## pv = unclass(mx$p)
      qv <- unclass(mx$q)
      nlq <- -log10(1e-99 + qv)
      ymax <- max(3, 1.2 * quantile(nlq, probs = 0.999, na.rm = TRUE)[1]) ## y-axis
      xlim <- c(-1.1, 1.1) * max(abs(fc))
      xlim <- 1.3 * c(-1, 1) * quantile(abs(fc), probs = 0.999)
      fc.genes <- ngs$genes[rownames(mx), "gene_name"]
      nplots <- min(24, ncol(qv))

      ## methods = names(ngs$gx.meta$output)
      methods <- colnames(ngs$gx.meta$meta[[1]]$fc)
      nc <- 6
      par(mfrow = c(2, 6), mar = c(4, 4, 2, 2) * 0, oma = c(1, 1, 0, 0) * 2)
      if (nplots > 12) {
        nplots <- min(nplots, 24)
        par(mfrow = c(3, 8), mar = c(4, 4, 2, 2) * 0)
        nc <- 8
      }

      shiny::withProgress(message = "computing volcano plots ...", value = 0, {
        i <- 1
        for (i in 1:nplots) {
          fx <- fc[, i]
          ## pval = pv[,i]
          qval <- qv[, i]
          sig.genes <- fc.genes[which(qval <= fdr & abs(fx) >= lfc)]
          ## genes1 = intersect(sig.genes, sel.genes)
          genes1 <- sig.genes[which(toupper(sig.genes) %in% toupper(sel.genes))]
          gx.volcanoPlot.XY(
            x = fx, pv = qval, gene = fc.genes,
            render = "canvas", n = 5000, nlab = 5,
            xlim = xlim, ylim = c(0, ymax), axes = FALSE,
            use.fdr = TRUE, p.sig = fdr, lfc = lfc,
            ## main=comp[i],
            ## ma.plot=TRUE, use.rpkm=TRUE,
            cex = 0.6, lab.cex = 1.5, highlight = genes1
          )

          is.first <- (i %% nc == 1)
          last.row <- ((i - 1) %/% nc == (nplots - 1) %/% nc)
          is.first
          last.row
          if (is.first) axis(2, mgp = c(2, 0.7, 0), cex.axis = 0.8)
          if (last.row) axis(1, mgp = c(2, 0.7, 0), cex.axis = 0.8)
          graphics::box(lwd = 1, col = "black", lty = "solid")
          legend("top",
            legend = colnames(fc)[i], cex = 1.2,
            bg = "white", box.lty = 0, inset = c(0, 0.01),
            x.intersp = 0.1, y.intersp = 0.1
          )
          shiny::incProgress(1 / length(nplots))
        }
      })
    })

    volcanoMethods_text <- "Under the <strong>Volcano (methods)</strong> tab, the platform displays the volcano plots provided by multiple differential expression calculation methods for the selected contrast. This provides users an overview of the statistics of all methods at the same time."

    shiny::callModule(plotModule,
      id = "volcanoMethods",
      func = volcanoMethods.RENDER,
      func2 = volcanoMethods.RENDER,
      title = "Volcano plots for all methods",
      info.text = volcanoMethods_text,
      ## caption = volcanoMethods_caption,
      height = c(imgH, 450), width = c("auto", 1600),
      res = c(75, 95),
      pdf.width = 18, pdf.height = 6,
      add.watermark = WATERMARK
    )

    ## ================================================================================
    ## Statistics Table
    ## ================================================================================

    gene_selected <- shiny::reactive({
      i <- as.integer(genetable$rows_selected())
      if (is.null(i) || length(i) == 0) {
        return(NULL)
      }
      res <- filteredDiffExprTable()
      gene <- rownames(res)[i]
      return(gene)
    })

    genetable.RENDER <- shiny::reactive({

      res <- filteredDiffExprTable()
      ## res <- fullDiffExprTable()

      if (is.null(res) || nrow(res) == 0) {
        return(NULL)
      }

      fx.col <- grep("fc|fx|mean.diff|logfc|foldchange", tolower(colnames(res)))[1]
      fx.col
      fx <- res[, fx.col]

      if ("gene_title" %in% colnames(res)) res$gene_title <- shortstring(res$gene_title, 50)
      rownames(res) <- sub(".*:", "", rownames(res))

      if (!DEV) {
        kk <- grep("meta.fx|meta.fc|meta.p", colnames(res), invert = TRUE)
        res <- res[, kk, drop = FALSE]
      }
      if (!input$gx_showqvalues) {
        kk <- grep("^q[.]", colnames(res), invert = TRUE)
        res <- res[, kk, drop = FALSE]
      }

      numeric.cols <- which(sapply(res, is.numeric))
      numeric.cols <- colnames(res)[numeric.cols]

      DT::datatable(res,
        rownames = FALSE,
        ## class = 'compact cell-border stripe hover',
        class = "compact hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          paging = TRUE,
          pageLength = 16, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE,
          scrollY = FALSE,
          scroller = FALSE,
          deferRender = TRUE,
          search = list(
            regex = TRUE,
            caseInsensitive = TRUE
            ## , search = 'M[ae]'
          )
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(colnames(res)[fx.col],
          ## background = DT::styleColorBar(c(0,3), 'lightblue'),
          background = color_from_middle(fx, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }) %>%
      bindCache(filteredDiffExprTable(), input$gx_showqvalues)

    genetable_text <- "Table <strong>I</strong> shows the results of the statistical tests. To increase the statistical reliability of the Omics Playground, we perform the DE analysis using four commonly accepted methods in the literature, namely, T-test (standard, Welch), <a href='https://www.ncbi.nlm.nih.gov/pubmed/25605792'> limma</a> (no trend, trend, voom), <a href='https://www.ncbi.nlm.nih.gov/pubmed/19910308'> edgeR</a> (QLF, LRT), and <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4302049'> DESeq2</a> (Wald, LRT), and merge the results.
<br><br>For a selected comparison under the <code>Contrast</code> setting, the results of the selected methods are combined and reported under the table, where <code>meta.q</code> for a gene represents the highest <code>q</code> value among the methods and the number of stars for a gene indicate how many methods identified significant <code>q</code> values (<code>q < 0.05</code>). The table is interactive (scrollable, clickable); users can sort genes by <code>logFC</code>, <code>meta.q</code>, or average expression in either conditions. Users can filter top N = {10} differently expressed genes in the table by clicking the <code>top 10 genes</code> from the table <i>Settings</i>."

    genetable_opts <- shiny::tagList(
      withTooltip(shiny::checkboxInput(ns("gx_top10"), "top 10 up/down genes", FALSE),
        "Display only top 10 differentially (positively and negatively) expressed genes in the table.",
        placement = "top", options = list(container = "body")
      ),
      withTooltip(shiny::checkboxInput(ns("gx_showqvalues"), "show indivivual q-values", FALSE),
        "Show q-values of each indivivual statistical method in the table.",
        placement = "top", options = list(container = "body")
      )
    )

    genetable <- shiny::callModule(
      tableModule,
      id = "genetable",
      func = genetable.RENDER,
      info.text = genetable_text,
      label = "I", info.width = "500px",
      options = genetable_opts,
      server = TRUE,
      title = "Differential expression analysis",
      height = c(tabH - 10, 700)
    )
    ## output$genetable <- genetable_module$render

    ## NEED RETHINK: reacts too often
    gx_related_genesets <- shiny::reactive({

      ngs <- inputData()
      res <- filteredDiffExprTable()
      if (is.null(res) || nrow(res) == 0) {
        return(NULL)
      }
      contr <- input$gx_contrast
      if (is.null(contr)) {
        return(NULL)
      }
      ## get table
      sel.row <- 1
      ## sel.row = input$genetable_rows_selected
      sel.row <- genetable$rows_selected()
      if (is.null(sel.row)) {
        return(NULL)
      }
      gene0 <- rownames(res)[sel.row]
      gene1 <- toupper(sub(".*:", "", gene0)) ## always uppercase...

      j <- which(toupper(rownames(ngs$GMT)) == gene1)
      gset <- names(which(ngs$GMT[j, ] != 0))
      gset <- intersect(gset, rownames(ngs$gsetX))
      if (length(gset) == 0) {
        return(NULL)
      }

      fx <- ngs$gset.meta$meta[[contr]]$meta.fx
      names(fx) <- rownames(ngs$gset.meta$meta[[contr]])
      fx <- round(fx[gset], digits = 4)

      rho <- cor(t(ngs$gsetX[gset, ]), ngs$X[gene0, ])[, 1]
      rho <- round(rho, digits = 3)
      gset1 <- substring(gset, 1, 60)

      df <- data.frame(geneset = gset1, rho = rho, fx = fx, check.names = FALSE)
      rownames(df) <- gset
      df <- df[order(-abs(df$fx)), ]

      return(df)
    })

    gsettable.RENDER <- shiny::reactive({
      df <- gx_related_genesets()
      if (is.null(df)) {
        return(NULL)
      }

      df$geneset <- wrapHyperLink(df$geneset, rownames(df))

      DT::datatable(df,
        ## class = 'compact cell-border stripe',
        class = "compact",
        rownames = FALSE, escape = c(-1, -2),
        extensions = c("Scroller"),
        fillContainer = TRUE,
        options = list(
          ## dom = 'lfrtip',
          dom = "frtip",
          paging = TRUE,
          pageLength = 16, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE,
          ## scrollY = tabV,
          scrollY = FALSE,
          scroller = FALSE,
          deferRender = TRUE,
          search = list(
            regex = TRUE,
            caseInsensitive = TRUE
            ## search = 'GOBP:'
          )
        ), ## end of options.list
        selection = list(mode = "single", target = "row", selected = NULL)
      ) %>%
        ## formatSignif(1:ncol(df),4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("fx", background = color_from_middle(df$fx, "lightblue", "#f5aeae"))
      # }, server=FALSE)
    })

    gsettable_text <- "By clicking on a gene in the Table <code>I</code>, it is possible to see which genesets contain that gene in this table, and check the differential expression status in other comparisons from the <code>Gene in contrasts</code> plot under the <code>Plots</code> tab."

    gsettable <- shiny::callModule(
      tableModule,
      id = "gsettable",
      func = gsettable.RENDER,
      info.text = gsettable_text, label = "II",
      title = "Gene sets with gene",
      height = c(tabH - 10, 700), width = c("100%", 800)
    )

    ## ================================================================================
    ## Foldchange (all)
    ## ================================================================================

    fctable.RENDER <- shiny::reactive({
      ngs <- inputData()
      res <- filteredDiffExprTable()
      if (is.null(res) || nrow(res) == 0) {
        return(NULL)
      }

      ## F <- sapply(ngs$gx.meta$meta, function(x) unclass(x$fc)[,"trend.limma"])
      ## Q <- sapply(ngs$gx.meta$meta, function(x) x$meta.q)
      ## F <- sapply(ngs$gx.meta$meta, function(x) x$meta.fx)
      ## rownames(F)=rownames(Q)=rownames(ngs$gx.meta$meta[[1]])
      F <- metaFC()
      Q <- metaQ()

      fc.rms <- sqrt(F[, 1]**2)
      if (NCOL(F) > 1) {
        fc.rms <- round(sqrt(rowMeans(F**2)), digits = 4)
      }

      show.q <- TRUE
      show.q <- input$fctable_showq
      df <- NULL
      if (show.q) {
        F1 <- do.call(cbind, lapply(1:ncol(F), function(i) cbind(F[, i], Q[, i])))
        colnames(F1) <- as.vector(rbind(paste0("FC.", colnames(F)), paste0("q.", colnames(Q))))
        ## colnames(F1) <- sub("q.*","q",colnames(F1))
        df <- data.frame(gene = rownames(F), rms.FC = fc.rms, F1, check.names = FALSE)
      } else {
        F1 <- F
        colnames(F1) <- paste0("FC.", colnames(F))
        df <- data.frame(gene = rownames(F), rms.FC = fc.rms, F1, check.names = FALSE)
      }

      df <- df[intersect(rownames(df), rownames(res)), ] ## take intersection of current comparison
      df <- df[order(-df$rms.FC), ]
      colnames(df) <- gsub("_", " ", colnames(df)) ## so it allows wrap line
      colnames(F1) <- gsub("_", " ", colnames(F1)) ## so it allows wrap line
      qv.cols <- grep("^q", colnames(F1))
      fc.cols <- setdiff(which(colnames(df) %in% colnames(F1)), qv.cols)
      ## if(length(qv.cols)==0) qv = 0

      dt <- DT::datatable(df,
        rownames = FALSE,
        # class = 'compact cell-border stripe hover',
        class = "compact hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = c(1)),
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE,
          scrollY = tabV,
          scroller = TRUE, deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatSignif(columns = fc.cols, digits = 3) %>%
        DT::formatStyle("rms.FC",
          ## background = DT::styleColorBar(c(0,3), 'lightblue'),
          background = color_from_middle(fc.rms, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) %>%
        DT::formatStyle(fc.cols,
          ## background = DT::styleColorBar(c(0,3), 'lightblue'),
          background = color_from_middle(F, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )

      if (length(qv.cols) > 0) {
        dt <- dt %>%
          DT::formatSignif(columns = qv.cols, digits = 3)
      }

      dt
    })

    fctable_text <- "The <strong>Foldchange (all)</strong> tab reports the gene fold changes for all contrasts in the selected dataset."

    fctable_caption <- "<b>Differential expression (fold-change) across all contrasts.</b> The column `rms.FC` corresponds to the root-mean-square fold-change across all contrasts."

    fctable_opts <- shiny::tagList(
      withTooltip(shiny::checkboxInput(ns("fctable_showq"), "show q-values", TRUE),
        "Show q-values next to FC values.",
        placement = "right", options = list(container = "body")
      )
    )

    shiny::callModule(
      tableModule,
      id = "fctable",
      func = fctable.RENDER,
      title = "Gene fold changes for all contrasts",
      info.text = fctable_text,
      options = fctable_opts,
      caption = fctable_caption,
      height = c(tabH, 700)
    )

    ## ================================================================================
    ## FDR table
    ## ================================================================================

    FDRtable.RENDER <- shiny::reactive({

      methods <- GX.DEFAULTTEST
      methods <- input$gx_statmethod
      ## methods = input$gx_statmethod
      if (is.null(methods)) {
        return(NULL)
      }

      ## comp <- input$gx_contrast
      ngs <- inputData()

      kk <- rownames(ngs$gx.meta$sig.counts[[1]][[1]])
      kk <- intersect(methods, rownames(ngs$gx.meta$sig.counts[[1]][[1]]))
      counts.up <- ngs$gx.meta$sig.counts$up
      counts.down <- ngs$gx.meta$sig.counts$down
      counts.up <- lapply(counts.up, function(x) x[kk, , drop = FALSE])
      counts.down <- lapply(counts.down, function(x) x[kk, , drop = FALSE])
      for (i in 1:length(counts.up)) {
        rownames(counts.up[[i]]) <- paste0(names(counts.up)[i], "::", rownames(counts.up[[i]]))
        rownames(counts.down[[i]]) <- paste0(names(counts.down)[i], "::", rownames(counts.down[[i]]))
      }
      sig.up <- do.call(rbind, counts.up)
      sig.down <- do.call(rbind, counts.down)

      sig.up <- sig.up[order(rownames(sig.up)), , drop = FALSE]
      sig.down <- sig.down[order(rownames(sig.down)), , drop = FALSE]
      colnames(sig.up)[1] <- paste("UP   FDR = ", colnames(sig.up)[1])
      colnames(sig.down)[1] <- paste("DOWN   FDR = ", colnames(sig.down)[1])
      colnames(sig.down) <- paste0("  ", colnames(sig.down))
      sigcount <- cbind(sig.down, sig.up[rownames(sig.down), , drop = FALSE])
      dim(sigcount)
      maxsig <- 0.99 * max(sigcount, na.rm = TRUE)

      contr <- sub("::.*", "", rownames(sigcount))
      ## contr = rownames(sigcount)
      metd <- sub(".*::", "", rownames(sigcount))
      D <- data.frame(method = metd, contrast = contr, sigcount, check.names = FALSE)

      DT::datatable(D,
        rownames = FALSE,
        #                      class = 'compact cell-border stripe hover',
        class = "compact hover",
        fillContainer = TRUE,
        extensions = c("Scroller"),
        options = list(
          dom = "lfrtip",
          pageLength = 999, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE,
          scrollY = tabV,
          scroller = TRUE, deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(colnames(sig.up),
          background = DT::styleColorBar(c(0, maxsig), "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) %>%
        DT::formatStyle(colnames(sig.down),
          background = DT::styleColorBar(c(0, maxsig), "lightblue"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    FDRtable_text <- "The <strong>FDR table</strong> tab reports the number of significant genes at different FDR thresholds for all contrasts within the dataset."

    FDRtable_caption <- "<b>Number of significant genes versus FDR.</b> This table reports the number of significant genes at different FDR thresholds for all contrasts and methods. This enables to quickly see which methods are more sensitive. The left part of the table (in blue) correspond to the number of significant down-regulated genes, the right part (in red) correspond to the number of significant overexpressed genes."

    shiny::callModule(
      tableModule,
      id = "FDRtable",
      func = FDRtable.RENDER,
      info.text = FDRtable_text,
      title = "Number of significant genes",
      caption = FDRtable_caption,
      height = c(tabH, 700)
    )

    ## ----------------------------------------------------------------------
    ## reactive values to return to parent environment
    ## ----------------------------------------------------------------------

    metaQ <- shiny::reactive({
      ngs <- inputData()
      req(ngs)
      methods <- selected_gxmethods()
      metaQ <- sapply(ngs$gx.meta$meta, function(m) apply(m$q[, methods, drop = FALSE], 1, max, na.rm = TRUE))
      rownames(metaQ) <- rownames(ngs$gx.meta$meta[[1]])
      metaQ
    })

    metaFC <- shiny::reactive({
      ngs <- inputData()
      req(ngs)
      methods <- selected_gxmethods()
      ## metaFC <- sapply(ngs$gx.meta$meta, function(m) rowMeans(m$fc[,methods,drop=FALSE]))
      metaFC <- sapply(ngs$gx.meta$meta, function(m) m$meta.fx)
      rownames(metaFC) <- rownames(ngs$gx.meta$meta[[1]])
      metaFC
    })


    outx <- list(
      selected_gxmethods = selected_gxmethods
    )
    return(outx)
  }) ## end of moduleServer
} ## end-of-ExpressionBoard
