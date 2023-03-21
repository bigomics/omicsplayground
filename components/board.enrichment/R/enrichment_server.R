##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

EnrichmentBoard <- function(id, inputData, selected_gxmethods) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 800
    rowH <- 420 ## row height of panels
    imgH <- 340 ## height of images
    tabV <- "70vh" ## height of tables
    tabH <- 340 ## row height of panels
    tabH <- "80vh" ## height of tables

    gs_infotext <- paste("Similar to the differential gene expression analysis, users can perform differential
        expression analysis on a geneset level in this page, which is also referred as gene set enrichment (GSE) analysis.
        The platform has more than 50.000 genesets (or pathways) in total that are divided into 30 geneset collections
        such as ", a_Hallmark, ", ", a_MSigDB, ", ", a_KEGG, " and ", a_GO, ". Users have to specify which comparison they want to
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
        allowfullscreen></iframe></center>")

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
      ngs <- inputData()
      shiny::req(ngs)
      meta <- ngs$gset.meta$meta
      comparisons <- colnames(ngs$model.parameters$contr.matrix)
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
      ngs <- inputData()
      shiny::req(ngs)
      nn <- sapply(COLLECTIONS, function(k) sum(k %in% rownames(ngs$gsetX)))
      gsets.groups <- names(COLLECTIONS)[which(nn >= 5)]
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
      ngs <- inputData()
      shiny::req(ngs)
      gset.methods0 <- colnames(ngs$gset.meta$meta[[1]]$fc)
      ## test = head(intersect(GSET.DEFAULTMETHODS,gset.methods0),3) ## maximum three
      test <- input$gs_statmethod
      test <- intersect(test, gset.methods0) ## maximum three
      test
    })

    calcGsetMeta <- function(comparison, methods, ngs) {
      ## ngs <- inputData()
      mx <- ngs$gset.meta$meta[[comparison]]
      if (is.null(mx)) {
        return(NULL)
      }
      mx.methods <- colnames(unclass(mx$fc))
      mx.methods
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
      dim(pv)
      if (NCOL(pv) > 1) {
        ss.rank <- function(x) scale(sign(x) * rank(abs(x)), center = FALSE)
        ## fc = rowMeans(scale(fc,center=FALSE),na.rm=TRUE)  ## REALLY???
        fc <- rowMeans(fc, na.rm = TRUE) ## NEED RETHINK!!!
        ## pv = apply(pv,1,function(x) metap::allmetap(x,method="sumz")$p[[1]])
        ## pv = apply(pv,1,vec.combinePvalues,method="stouffer")
        ## qv = p.adjust(pv, method="fdr")
        pv <- apply(pv, 1, max, na.rm = TRUE)
        qv <- apply(qv, 1, max, na.rm = TRUE)
        ## score = rowMeans(scale(score,center=FALSE),na.rm=TRUE)
        score <- rowMeans(apply(score, 2, ss.rank), na.rm = TRUE)
      }

      meta <- cbind(score = score, fc = fc, pv = pv, qv = qv)
      rownames(meta) <- rownames(mx)
      colnames(meta) <- c("score", "fc", "pv", "qv") ## need
      return(meta)
    }

    getFullGeneSetTable <- shiny::reactive({
      ngs <- inputData()
      shiny::req(ngs)
      comp <- 1
      comp <- input$gs_contrast
      if (is.null(comp)) {
        return(NULL)
      }
      if (!(comp %in% names(ngs$gset.meta$meta))) {
        return(NULL)
      }
      mx <- ngs$gset.meta$meta[[comp]]
      dim(mx)

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
      if (is.null(input$gs_features)) {
        return(NULL)
      }
      if (1 && !(gsfeatures %in% c(NA, "", "*", "<all>")) &&
        gsfeatures %in% names(COLLECTIONS)) {
        sel <- intersect(rownames(mx), COLLECTIONS[[gsfeatures]])
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
        stars <- sapply(rowSums(is.sig, na.rm = TRUE), star.symbols, pch = "\u2605")
        names(stars) <- rownames(mx)

        ## ------------ calculate META parameters ----------------
        meta <- calcGsetMeta(comp, gsmethod, ngs = ngs)
        meta <- meta[rownames(mx), , drop = FALSE]
        dim(meta)
        gset.size <- Matrix::colSums(ngs$GMT[, rownames(mx), drop = FALSE] != 0)
        names(gset.size) <- rownames(mx)

        ## ---------- report *average* group expression FOLD CHANGE
        ## THIS SHOULD BETTER GO DIRECTLY WHEN CALCULATING GSET TESTS
        ##
        s1 <- names(which(ngs$model.parameters$exp.matrix[, comp] > 0))
        s0 <- names(which(ngs$model.parameters$exp.matrix[, comp] < 0))
        jj <- colnames(ngs$GMT)
        jj <- rownames(mx)

        gsdiff.method <- "fc" ## OLD default
        if (gsdiff.method == "gs") {
          AveExpr1 <- rowMeans(ngs$gsetX[jj, s1])
          AveExpr0 <- rowMeans(ngs$gsetX[jj, s0])
          meta.fc <- AveExpr1 - AveExpr0
        } else {
          ## WARNING!!! THIS STILL ASSUMES GENES AS rownames(ngs$X)
          ## and rownames(GMT)
          fc <- ngs$gx.meta$meta[[comp]]$meta.fx ## stable
          names(fc) <- rownames(ngs$gx.meta$meta[[comp]])
          pp <- intersect(rownames(ngs$GMT), names(fc))
          rnaX <- ngs$X

          ## check if multi-omics
          is.multiomics <- any(grepl("\\[gx\\]|\\[mrna\\]", names(fc)))
          is.multiomics
          if (is.multiomics) {
            ii <- grep("\\[gx\\]|\\[mrna\\]", names(fc))
            fc <- fc[ii]
            rnaX <- ngs$X[names(fc), ]
            names(fc) <- sub(".*:|.*\\]", "", names(fc))
            rownames(rnaX) <- sub(".*:|.*\\]", "", rownames(rnaX))
            pp <- intersect(rownames(ngs$GMT), names(fc))
            length(pp)
          }

          G <- Matrix::t(ngs$GMT[pp, jj] != 0)
          ngenes <- Matrix::rowSums(G)
          meta.fc <- ngs$gset.meta$meta[[comp]]$meta.fx
          names(meta.fc) <- rownames(ngs$gset.meta$meta[[comp]])

          AveExpr1 <- Matrix::rowMeans(G %*% rnaX[pp, s1]) / ngenes
          AveExpr0 <- Matrix::rowMeans(G %*% rnaX[pp, s0]) / ngenes
          remove(rnaX)
        }

        ## TWIDDLE means to reflect foldchange...
        mean0 <- (AveExpr0 + AveExpr1) / 2
        AveExpr1 <- mean0 + meta.fc / 2
        AveExpr0 <- mean0 - meta.fc / 2

        ##
        gs <- intersect(names(meta.fc), rownames(meta))
        length(gs)

        rpt <- data.frame(
          size = gset.size[gs],
          logFC = meta.fc[gs],
          meta.q = meta[gs, "qv"],
          stars = stars[gs],
          AveExpr0 = AveExpr0[gs],
          AveExpr1 = AveExpr1[gs]
        )

        ## add extra p/q value columns
        jj <- match(gs, rownames(mx))
        rpt <- cbind(rpt, q = qv[jj, ])

        ## rownames(rpt) = gs
      } else {
        ## show original table (single method)
        rpt <- outputs[[gsmethod]]
      }

      rpt <- rpt[order(-rpt$logFC), ] ## positive
      rpt <- data.frame(rpt)

      return(rpt)
    })


    getFilteredGeneSetTable <- shiny::reactive({
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
        shiny::validate(shiny::need(nrow(res) > 0, "warning. no genesets passed current filters."))
        return(NULL)
      }
      return(res)
    })

    ## ----------------------------------------------------------------------
    ## 0: Volcano plot in gene space
    ## ----------------------------------------------------------------------
    subplot.MAR <- c(3, 3.5, 1.5, 0.5)
    subplot.MAR <- c(2.8, 4, 4, 0.8)

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
      ##

      ngs <- inputData()
      shiny::req(ngs, input$gs_contrast)
      gs <- 1
      comp <- 1

      comp <- input$gs_contrast
      gs <- gset_selected()
      if (is.null(gs) || length(gs) == 0) {
        return(NULL)
      }

      mx <- ngs$gx.meta$meta[[comp]]
      is.multiomics <- any(grepl("\\[gx\\]|\\[mrna\\]", rownames(mx)))
      is.multiomics
      if (is.multiomics) {
        ii <- grep("\\[gx\\]|\\[mrna\\]", rownames(mx))
        mx <- mx[ii, ]
        ## rownames(mx) <- sub(".*:|.*\\]","",rownames(mx))
      }

      ## gxmethods <- c("trend.limma","ttest.welch")
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
      gmt1 <- ngs$GMT[, gs, drop = FALSE]
      genes <- rownames(gmt1)[which(Matrix::rowSums(gmt1 != 0) == ns)]
      genes <- intersect(genes, ngs$genes[rownames(limma1), "gene_name"])
      genes <- setdiff(genes, c("", NA, "NA", " "))

      title <- rep(NA, length(genes))
      title <- as.character(GENE.TITLE[genes])
      title[is.na(title)] <- " "

      rpt <- data.frame("gene_name" = genes, "gene_title" = as.character(title))
      genes <- rpt[, "gene_name"]
      genes1 <- ngs$genes[rownames(limma1), "gene_name"]
      limma1 <- limma1[match(genes, genes1), , drop = FALSE] ## align limma1
      rpt <- cbind(rpt, limma1)
      rpt <- rpt[which(!is.na(rpt$fc) & !is.na(rownames(rpt))), , drop = FALSE]

      if (nrow(rpt) > 0) {
        rpt <- rpt[order(-abs(rpt$fc)), , drop = FALSE]
      }
      return(rpt)
    })

    gene_selected <- shiny::reactive({
      ngs <- inputData()
      shiny::req(ngs)
      i <- as.integer(genetable$rows_selected())
      if (is.null(i) || is.na(i) || length(i) == 0) i <- 1
      rpt <- geneDetails()
      if (is.null(rpt) || nrow(rpt) == 0) {
        return(list(gene = NA, probe = NA))
      }
      sel.gene <- rownames(rpt)[i]
      gene <- as.character(rpt$gene_name[i])
      probe <- rownames(ngs$genes)[match(gene, ngs$genes$gene_name)]
      return(list(gene = gene, probe = probe))
    })

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    WATERMARK <- FALSE

    # Top enriched gene sets

    enrichment_plot_top_enrich_gsets_server(
      "topEnriched",
      inputData = inputData,
      getFilteredGeneSetTable = getFilteredGeneSetTable,
      gs_contrast = shiny::reactive(input$gs_contrast),
      gseatable = gseatable,
      watermark = WATERMARK
    )

    # Frequency in top gene sets

    enrichment_plot_freq_top_gsets_server(
      "topEnrichedFreq",
      inputData = inputData,
      getFilteredGeneSetTable = getFilteredGeneSetTable,
      gs_contrast = shiny::reactive(input$gs_contrast),
      gseatable = gseatable,
      watermark = WATERMARK
    )

    # Volcano plot

    enrichment_plot_volcano_server(
      "subplot_volcano",
      inputData = inputData,
      gs_contrast = shiny::reactive(input$gs_contrast),
      selected_gxmethods = selected_gxmethods,
      gset_selected = gset_selected,
      gs_fdr = shiny::reactive(input$gs_fdr),
      gs_lfc = shiny::reactive(input$gs_lfc),
      subplot.MAR = subplot.MAR,
      watermark = WATERMARK
    )

    # Enrichment barplot

    enrichment_plot_barplot_server(
      "subplot_barplot",
      inputData = inputData,
      gset_selected = gset_selected,
      gs_contrast = shiny::reactive(input$gs_contrast),
      subplot.MAR = subplot.MAR,
      watermark = WATERMARK
    )

    # Expression geneplot

    enrichment_plot_geneplot_server(
      "subplot_geneplot",
      inputData = inputData,
      gs_contrast = shiny::reactive(input$gs_contrast),
      gene_selected = gene_selected,
      subplot.MAR = subplot.MAR,
      watermark = WATERMARK
    )

    # Enrichment vs. expression

    enrichment_plot_scatter_server(
      "subplot_scatter",
      inputData = inputData,
      gene_selected = gene_selected,
      gs_contrast = shiny::reactive(input$gs_contrast),
      subplot.MAR = subplot.MAR,
      gset_selected = gset_selected,
      watermark = WATERMARK
    )

    # Enrichment of geneset across multiple contrasts

    enrichment_plot_compare_server(
      "compare",
      inputData = inputData,
      gs_contrast = shiny::reactive(input$gs_contrast),
      gset_selected = gset_selected,
      selected_gsetmethods = selected_gsetmethods,
      watermark = WATERMARK
    )

    # Volcano plots for all contrasts

    enrichment_plot_volcanoall_server(
      "volcanoAll",
      inputData = inputData,
      gs_features = shiny::reactive(input$gs_features),
      gs_statmethod = shiny::reactive(input$gs_statmethod),
      gs_fdr = shiny::reactive(input$gs_fdr),
      gs_lfc = shiny::reactive(input$gs_lfc),
      calcGsetMeta = calcGsetMeta,
      watermark = WATERMARK
    )

    # Volcano plots for all methods

    enrichment_plot_volcanomethods_server(
      "volcanoMethods",
      inputData = inputData,
      gs_features = shiny::reactive(input$gs_features),
      gs_contrast = shiny::reactive(input$gs_contrast),
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
      geneDetails = geneDetails
    )

    # Gene set enrichment for all contrasts

    enrichment_table_gset_enrich_all_contrasts_server(
      "fctable",
      inputData = inputData,
      getFilteredGeneSetTable = getFilteredGeneSetTable
    )

    # Number of significant gene sets

    enrichment_table_n_sig_gsets_server(
      "FDRtable",
      inputData = inputData,
      gs_statmethod = shiny::reactive(input$gs_statmethod)
    )

    ## reactive values to return to parent environment
    outx <- list(selected_gsetmethods = selected_gsetmethods)
    return(outx)
  }) ## end of moduleServer
} ## end-of-Board
