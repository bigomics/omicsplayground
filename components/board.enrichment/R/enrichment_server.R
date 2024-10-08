##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

EnrichmentBoard <- function(id, pgx, selected_gxmethods = reactive(colnames(pgx$gx.meta$meta[[1]]$fc))) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 800
    rowH <- 420 ## row height of panels
    imgH <- 340 ## height of images
    tabV <- "70vh" ## height of tables
    tabH <- 340 ## row height of panels
    tabH <- "80vh" ## height of tables

    gs_infotext <- tspan(paste("Similar to the differential gene expression analysis, users can perform differential
        expression analysis on a geneset level in this page, which is also referred as gene set enrichment (GSE) analysis.
        The platform has more than 50.000 genesets (or pathways) in total that are divided into 30 geneset collections
        such as ", a_Hallmark, ", ", a_MSigDB, " and ", a_GO, ". Users have to specify which comparison they want to
        visually analyze employing a certain geneset collection.<br><br>
        To ensure the statistical reliability, the platform performs Enrichment Analyses using multiple methods.
        The result from the statistical methods is displayed in <strong>Enrichment table</strong> panel.
        In the <strong>Top enriched</strong> panel, the top 10 differentially enriched geneses (pathways) are displayed.
        In the <strong>Plots</strong> panel, a volcano plot of genes contained in the selected geneset and a barplot
        of expressions per sample group are displayed. In the <strong>Compare</strong> panel, users can compare the
        differential expression status of that geneset for all other comparisons. Finally, volcano plots of genesets
        for all comparisons are displayed under the <strong>Volcano (all) </strong> tab. This allows users to have
        an overall picture across comparisons at the same time.<br><br>
        EXPERT MODE ONLY: To compare the different statistical methods, the <strong>Volcano (methods)</strong>
        panel shows volcano plots of all methods. The <strong>FDR table</strong> panel reports the number of
        significant gene sets at different FDR thresholds for all contrasts.<br><br><br><br>
        <center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=4'
        frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture'
        allowfullscreen></iframe></center>"), js = FALSE)

    GSET.DEFAULTMETHODS <- c("gsva", "camera", "fgsea", "fisher")

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$gs_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Enrichment Analysis Board</strong>"),
        shiny::HTML(gs_infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    shiny::observe({
      shiny::req(pgx$X)
      meta <- pgx$gset.meta$meta
      comparisons <- colnames(pgx$model.parameters$contr.matrix)
      comparisons <- sort(intersect(comparisons, names(meta)))
      shiny::updateSelectInput(session, "gs_contrast", choices = comparisons)

      ## get the computed geneset methods
      gset.methods <- sort(colnames(meta[[1]]$fc))
      sel2 <- c(intersect(GSET.DEFAULTMETHODS, gset.methods), gset.methods)
      sel2 <- head(unique(sel2), 3)

      shiny::updateCheckboxGroupInput(session, "gs_statmethod",
        choices = sort(gset.methods),
        selected = sel2
      )
    })

    shiny::observe({
      shiny::req(pgx$X)
      gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
      nn <- sapply(gset_collections, function(k) sum(k %in% rownames(pgx$gsetX)))
      gsets.groups <- names(gset_collections)[which(nn >= 5)]
      gsets.groups <- c("<all>", sort(gsets.groups))
      sel <- "<all>"
      hmark <- grep("^H$|hallmark|", gsets.groups, ignore.case = TRUE, value = TRUE)
      if (length(hmark) > 0) sel <- hmark[1]
      shiny::updateSelectInput(session, "gs_features", choices = gsets.groups, selected = sel)
    })

    ## ================================================================================
    ## ========================= REACTIVE FUNCTIONS ===================================
    ## ================================================================================

    selected_gsetmethods <- shiny::reactive({
      shiny::req(pgx$X)
      gset.methods0 <- colnames(pgx$gset.meta$meta[[1]]$fc)
      test <- input$gs_statmethod
      test <- intersect(test, gset.methods0)
      test
    })

    calcGsetMeta <- function(comparison, methods, pgx) {
      mx <- pgx$gset.meta$meta[[comparison]]
      shiny::req(mx)
      mx.methods <- colnames(unclass(mx$fc))
      methods <- intersect(methods, mx.methods)
      if (is.null(methods) || length(methods) == 0) {
        cat("ERROR: calcGsetMeta:: no valid methods\n")
        return(NULL)
      }

      ## recalculate meta values
      pv <- unclass(mx$p)[, methods, drop = FALSE]
      qv <- unclass(mx$q)[, methods, drop = FALSE]
      fc <- unclass(mx$fc)[, methods, drop = FALSE]

      ## !!!! Because the methods have all very difference "fold-change" !!!!
      ## estimators, we use the meta.fx (average of all genes in gset)
      fc <- do.call(cbind, rep(list(mx$meta.fx), length(methods)))
      colnames(fc) <- methods

      pv[is.na(pv)] <- 1
      qv[is.na(qv)] <- 1
      fc[is.na(fc)] <- 0
      score <- fc * (-log10(qv))
      if (NCOL(pv) > 1) {
        ss.rank <- function(x) scale(sign(x) * rank(abs(x)), center = FALSE)
        fc <- rowMeans(fc, na.rm = TRUE) ## NEED RETHINK!!!
        pv <- apply(pv, 1, max, na.rm = TRUE)
        qv <- apply(qv, 1, max, na.rm = TRUE)
        score <- rowMeans(apply(score, 2, ss.rank), na.rm = TRUE)
      }

      meta <- cbind(score = score, fc = fc, pv = pv, qv = qv)
      rownames(meta) <- rownames(mx)
      colnames(meta) <- c("score", "fc", "pv", "qv") ## need
      return(meta)
    }

    getFullGeneSetTable <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(input$gs_contrast)
      shiny::req(input$gs_features)
      comp <- input$gs_contrast
      gene_symbols <- pgx$genes[rownames(pgx$gx.meta$meta[[comp]]), "symbol"]
      names(gene_symbols) <- rownames(pgx$gx.meta$meta[[comp]])

      if (!(comp %in% names(pgx$gset.meta$meta))) {
        return(NULL)
      }
      mx <- pgx$gset.meta$meta[[comp]]

      outputs <- NULL
      gsmethod <- colnames(unclass(mx$fc))
      gsmethod <- input$gs_statmethod
      if (is.null(gsmethod) || length(gsmethod) == 0) {
        return(NULL)
      }

      lfc <- as.numeric(input$gs_lfc)
      fdr <- as.numeric(input$gs_fdr)

      ## filter gene sets for table
      gsfeatures <- "<all>"
      gsfeatures <- input$gs_features
      gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
      if (1 && !(gsfeatures %in% c(NA, "", "*", "<all>")) &&
        gsfeatures %in% names(gset_collections)) {
        sel <- intersect(rownames(mx), gset_collections[[gsfeatures]])
        mx <- mx[sel, , drop = FALSE]
      }

      rpt <- NULL

      if (is.null(outputs) || length(gsmethod) > 1) {
        ## show meta-statistics table (multiple methods)
        pv <- unclass(mx$p)[, gsmethod, drop = FALSE]
        qv <- unclass(mx$q)[, gsmethod, drop = FALSE]
        fx <- unclass(mx$fc)[, gsmethod, drop = FALSE]

        ## !!!! Because the methods have all very difference "fold-change" !!!!
        ## estimators, we use the meta.fx (average of all genes in gset)
        fx <- do.call(cbind, rep(list(mx$meta.fx), length(gsmethod)))
        colnames(fx) <- gsmethod

        pv[is.na(pv)] <- 1
        qv[is.na(qv)] <- 1
        fx[is.na(fx)] <- 0

        is.sig <- (qv <= fdr & abs(fx) >= lfc)
        stars <- sapply(rowSums(is.sig, na.rm = TRUE), playbase::star.symbols, pch = "\u2605")
        names(stars) <- rownames(mx)

        ## ------------ calculate META parameters ----------------
        meta <- calcGsetMeta(comp, gsmethod, pgx = pgx)
        meta <- meta[rownames(mx), , drop = FALSE]

        if (!is.data.frame(pgx$gset.meta$info)) {
          # below is legacy code, before branch dev-gmt 2023-05
          # in new pgx code, this will be automatically calculated, with more details
          gset.size <- Matrix::colSums(pgx$GMT[, rownames(mx), drop = FALSE] != 0)
          names(gset.size) <- rownames(mx)
        }

        ## ---------- report *average* group expression FOLD CHANGE
        ## THIS SHOULD BETTER GO DIRECTLY WHEN CALCULATING GSET TESTS
        ##
        s1 <- names(which(pgx$model.parameters$exp.matrix[, comp] > 0))
        s0 <- names(which(pgx$model.parameters$exp.matrix[, comp] < 0))
        jj <- rownames(mx)
        jj <- intersect(jj, colnames(pgx$GMT))
        rnaX <- pgx$X
        rnaX <- playbase::rename_by(rnaX, pgx$genes, "symbol")
        gsdiff.method <- "fc" ## OLD default

        if (gsdiff.method == "gs") {
          AveExpr1 <- rowMeans(pgx$gsetX[jj, s1], na.rm = TRUE)
          AveExpr0 <- rowMeans(pgx$gsetX[jj, s0], na.rm = TRUE)
          meta.fc <- AveExpr1 - AveExpr0
        } else {
          ## WARNING!!! THIS STILL ASSUMES GENES AS rownames(pgx$X)
          ## and rownames(GMT)
          fc <- pgx$gx.meta$meta[[comp]]$meta.fx
          names(fc) <- rownames(pgx$gx.meta$meta[[comp]])
          pp <- intersect(rownames(pgx$GMT), names(fc))

          # if pp is null, use human ortholog.
          # pp is null when collapse by gene is false, but we dont have a parameter for that.
          if (length(pp) == 0 || is.null(pp)) {
            names(fc) <- pgx$genes[names(fc), "symbol"]
            pp <- intersect(pgx$genes$symbol, names(fc))
            pp <- intersect(rownames(pgx$GMT), names(fc))

            pp <- pp[!is.na(pp)]
          }

          ## check if multi-omics (TEMPORARILY FALSE)
          is.multiomics <- FALSE # any(grepl("\\[gx\\]|\\[mrna\\]", names(fc)))
          if (is.multiomics) {
            ii <- grep("\\[gx\\]|\\[mrna\\]", names(fc))
            fc <- fc[ii]
            rnaX <- rnaX[names(fc), ]
            names(fc) <- sub(".*:|.*\\]", "", names(fc))
            rownames(rnaX) <- sub(".*:|.*\\]", "", rownames(rnaX))
            pp <- intersect(rownames(pgx$GMT), names(fc))
          }

          G <- Matrix::t(pgx$GMT[pp, jj] != 0)
          ngenes <- Matrix::rowSums(G, na.rm = TRUE)
          meta.fc <- pgx$gset.meta$meta[[comp]]$meta.fx
          names(meta.fc) <- rownames(pgx$gset.meta$meta[[comp]])

          # subset rnaX by pp
          AveExpr1 <- Matrix::rowMeans(G %*% rnaX[pp, s1], na.rm = TRUE) / ngenes
          AveExpr0 <- Matrix::rowMeans(G %*% rnaX[pp, s0], na.rm = TRUE) / ngenes
          remove(rnaX)
        }

        ## TWIDDLE means to reflect foldchange...
        mean0 <- (AveExpr0 + AveExpr1) / 2
        AveExpr1 <- mean0 + meta.fc / 2
        AveExpr0 <- mean0 - meta.fc / 2

        # Subset meta.fc with non-na values
        meta.fc <- meta.fc[which(!is.na(meta.fc) & !is.na(mean0))]
        gs <- intersect(names(meta.fc), rownames(meta))

        if (!is.data.frame(pgx$gset.meta$info)) {
          rpt <- data.frame(
            logFC = meta.fc[gs],
            meta.q = meta[gs, "qv"],
            matched.genes = gset.size[gs],
            stars = stars[gs],
            AveExpr0 = AveExpr0[gs],
            AveExpr1 = AveExpr1[gs]
          )
        }

        if (is.data.frame(pgx$gset.meta$info)) {
          rpt <- data.frame(
            logFC = meta.fc[gs],
            meta.q = meta[gs, "qv"],
            matched.genes = pgx$gset.meta$info[gs, "gset.size"],
            total.genes = pgx$gset.meta$info[gs, "gset.size.raw"],
            fraction.genes.covered = pgx$gset.meta$info[gs, "gset.fraction"],
            stars = stars[gs],
            AveExpr0 = AveExpr0[gs],
            AveExpr1 = AveExpr1[gs]
          )
        }

        ## add extra p/q value columns
        jj <- match(gs, rownames(mx))
        rpt <- cbind(rpt, q = qv[jj, ])

        #
      } else {
        ## show original table (single method)
        rpt <- outputs[[gsmethod]]
      }

      rpt <- rpt[order(rpt$meta.q, -rpt$logFC), ] ## positive
      rpt <- data.frame(rpt)

      return(rpt)
    })

    getFilteredGeneSetTable <- shiny::reactive({
      req(getFullGeneSetTable())
      if (is.null(input$gs_showall) || length(input$gs_showall) == 0) {
        return(NULL)
      }
      if (is.null(input$gs_top10) || length(input$gs_top10) == 0) {
        return(NULL)
      }

      res <- getFullGeneSetTable()

      ## just show significant genes
      if (!input$gs_showall && nrow(res) > 0) {
        lfc <- as.numeric(input$gs_lfc)
        fdr <- as.numeric(input$gs_fdr)
        is.sig <- (abs(res$logFC) >= lfc & res$meta.q <= fdr)
        res <- res[is.sig, , drop = FALSE]
      }

      ## just show top 10
      if (input$gs_top10 && nrow(res) > 10 && length(input$gs_top10)) {
        ## get the meta-score column
        fx.col <- grep("score|fx|fc|sign|NES|logFC", colnames(res), value = TRUE)[1]
        fx <- as.numeric(res[, fx.col])
        names(fx) <- rownames(res)
        pp <- unique(c(
          head(names(sort(-fx[which(fx > 0)])), 10),
          head(names(sort(fx[which(fx < 0)])), 10)
        ))
        res <- res[pp, , drop = FALSE]
        fx <- as.numeric(res[, fx.col])
        res <- res[order(-fx), , drop = FALSE]
      }

      res <- data.frame(res)

      if (nrow(res) == 0) {
        shiny::validate(shiny::need(nrow(res) > 0, tspan("No genesets passed the statistical thresholds. Please update the thresholds on the settings sidebar.", js = FALSE)))
        return(NULL)
      }
      return(res)
    })


    metaQ <- shiny::reactive({
      req(pgx$X)
      methods <- selected_gsetmethods()
      metaQ <- sapply(pgx$gset.meta$meta, function(m) {
        apply(m$q[, methods, drop = FALSE], 1, max, na.rm = TRUE)
      })
      rownames(metaQ) <- rownames(pgx$gset.meta$meta[[1]])
      metaQ
    })

    metaFC <- shiny::reactive({
      req(pgx$X)
      #
      metaFC <- sapply(pgx$gset.meta$meta, function(m) m$meta.fx)
      rownames(metaFC) <- rownames(pgx$gset.meta$meta[[1]])
      metaFC
    })

    ## ================================================================================
    ## Enrichment table
    ## ================================================================================

    gset_selected <- shiny::reactive({
      i <- as.integer(gseatable$rows_selected())
      if (is.null(i) || length(i) == 0) {
        return(NULL)
      }
      rpt <- getFilteredGeneSetTable()
      gs <- rownames(rpt)[i]

      return(gs)
    })

    geneDetails <- shiny::reactive({
      ## return details of the genes in the selected gene set
      shiny::req(pgx$X, input$gs_contrast)
      gs <- 1
      comp <- 1
      comp <- input$gs_contrast
      gs <- gset_selected()
      if (is.null(gs) || length(gs) == 0) {
        return(NULL)
      }

      mx <- pgx$gx.meta$meta[[comp]]
      is.multiomics <- any(grepl("\\[gx\\]|\\[mrna\\]", rownames(mx)))
      if (is.multiomics) {
        ii <- grep("\\[gx\\]|\\[mrna\\]", rownames(mx))
        mx <- mx[ii, ]
      }

      gxmethods <- selected_gxmethods() ## from module-expression
      shiny::req(gxmethods)
      limma1.fc <- mx$meta.fx
      limma1.pq <- sapply(mx[, c("p", "q")], function(x) {
        apply(x[, gxmethods, drop = FALSE], 1, max, na.rm = TRUE)
      })
      limma1 <- cbind(fc = limma1.fc, limma1.pq)
      rownames(limma1) <- rownames(mx)

      ## in multi-mode we select *common* genes
      ns <- length(gs)
      gmt1 <- pgx$GMT[, gs, drop = FALSE]
      genes <- rownames(gmt1)[which(Matrix::rowSums(gmt1 != 0) == ns)]
      # check which columns are in pgx$genes
      cols_in_pgx <- c("feature", "symbol", "human_ortholog", "gene_title")
      cols_in_pgx <- cols_in_pgx[which(cols_in_pgx %in% colnames(pgx$genes))]

      genes_user <- pgx$genes[rownames(limma1), cols_in_pgx]
      empty_cols <- apply(genes_user, 2, function(x) all(is.na(x)))
      genes_user <- genes_user[, !empty_cols, drop = FALSE]
      genes <- genes_user[genes_user$symbol %in% genes, ]
      limma1 <- limma1[rownames(genes), , drop = FALSE] ## align limma1

      genes <- cbind(genes, limma1)
      genes <- genes[which(!is.na(genes$fc) & !is.na(rownames(genes))), , drop = FALSE]

      if (nrow(genes) > 0) {
        genes <- genes[order(-abs(genes$fc)), , drop = FALSE]
      }
      return(genes)
    })

    gene_selected <- shiny::reactive({
      shiny::req(pgx$X)
      i <- as.integer(genetable$rows_selected())
      rpt <- geneDetails()
      not.i <- (is.null(i) || is.na(i) || length(i) == 0)
      if (not.i || is.null(rpt) || nrow(rpt) == 0) {
        ## return(list(gene = NA, probe = NA))
        return(NULL)
      }
      sel.row <- rownames(rpt)[i]
      ##      res <- unlist(pgx$genes[sel.row, c("symbol", "feature")])
      ##      return(list(rn = sel.row, gene = res["symbol"], probe = res["feature"]))
      symbol <- pgx$genes[sel.row, "symbol"]
      return(list(rn = sel.row, gene = symbol))
    })

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    subplot.MAR <- c(2.8, 4, 4, 0.8)

    # Top enriched gene sets

    enrichment_plot_top_enrich_gsets_server(
      "topEnriched",
      pgx = pgx,
      getFilteredGeneSetTable = getFilteredGeneSetTable,
      gs_contrast = shiny::reactive(input$gs_contrast),
      gseatable = gseatable,
      watermark = WATERMARK
    )

    # Frequency in top gene sets

    enrichment_plot_freq_top_gsets_server(
      "topEnrichedFreq",
      pgx = pgx,
      getFilteredGeneSetTable = getFilteredGeneSetTable,
      gs_contrast = shiny::reactive(input$gs_contrast),
      gseatable = gseatable,
      watermark = WATERMARK
    )

    # Volcano plot

    enrichment_plot_volcano_server(
      "subplot_volcano",
      pgx = pgx,
      gs_contrast = shiny::reactive(input$gs_contrast),
      selected_gxmethods = selected_gxmethods,
      gset_selected = gset_selected,
      gs_fdr = shiny::reactive(input$gs_fdr),
      gs_lfc = shiny::reactive(input$gs_lfc),
      subplot.MAR = subplot.MAR,
      geneDetails = geneDetails,
      watermark = WATERMARK
    )

    # Enrichment barplot

    enrichment_plot_barplot_server(
      "subplot_barplot",
      pgx = pgx,
      gset_selected = gset_selected,
      gs_contrast = shiny::reactive(input$gs_contrast),
      subplot.MAR = subplot.MAR,
      watermark = WATERMARK
    )

    # Expression geneplot

    enrichment_plot_geneplot_server(
      "subplot_geneplot",
      pgx = pgx,
      gs_contrast = shiny::reactive(input$gs_contrast),
      gene_selected = gene_selected,
      subplot.MAR = subplot.MAR,
      watermark = WATERMARK
    )

    # Enrichment vs. expression

    enrichment_plot_scatter_server(
      "subplot_scatter",
      pgx = pgx,
      gene_selected = gene_selected,
      gs_contrast = shiny::reactive(input$gs_contrast),
      subplot.MAR = subplot.MAR,
      gset_selected = gset_selected,
      watermark = WATERMARK
    )

    # Enrichment of geneset across multiple contrasts

    enrichment_plot_compare_server(
      "compare",
      pgx = pgx,
      gs_contrast = shiny::reactive(input$gs_contrast),
      gset_selected = gset_selected,
      selected_gsetmethods = selected_gsetmethods,
      watermark = WATERMARK
    )

    # Volcano plots for all contrasts

    enrichment_plot_volcanoall_server(
      "volcanoAll",
      pgx = pgx,
      gs_features = shiny::reactive(input$gs_features),
      gs_statmethod = shiny::reactive(input$gs_statmethod),
      gs_fdr = shiny::reactive(input$gs_fdr),
      gs_lfc = shiny::reactive(input$gs_lfc),
      calcGsetMeta = calcGsetMeta,
      gset_selected = gset_selected,
      watermark = WATERMARK
    )

    # Volcano plots for all methods

    enrichment_plot_volcanomethods_server(
      "volcanoMethods",
      pgx = pgx,
      gs_features = shiny::reactive(input$gs_features),
      gs_contrast = shiny::reactive(input$gs_contrast),
      gs_fdr = shiny::reactive(input$gs_fdr),
      gs_lfc = shiny::reactive(input$gs_lfc),
      gset_selected = gset_selected,
      watermark = WATERMARK
    )

    # Enrichment analysis

    gseatable <- enrichment_table_enrichment_analysis_server(
      "gseatable",
      getFilteredGeneSetTable = getFilteredGeneSetTable
    )

    # Genes in gene set

    genetable <- enrichment_table_genes_in_geneset_server(
      "genetable",
      organism = pgx$organism,
      geneDetails = geneDetails
    )

    # Gene set enrichment for all contrasts

    enrichment_table_gset_enrich_all_contrasts_server(
      "fctable",
      pgx = pgx,
      getFilteredGeneSetTable = getFilteredGeneSetTable,
      metaFC = metaFC,
      metaQ = metaQ
    )

    # Number of significant gene sets

    enrichment_table_n_sig_gsets_server(
      "FDRtable",
      pgx = pgx,
      gs_statmethod = shiny::reactive(input$gs_statmethod)
    )

    ## reactive values to return to parent environment
    outx <- list(selected_gsetmethods = selected_gsetmethods)
    return(outx)
  }) ## end of moduleServer
} ## end-of-Board
