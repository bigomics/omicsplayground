##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ExpressionBoard <- function(id, pgx, labeltype = shiny::reactive("feature"),
                            board_observers = board_observers ) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 780
    rowH <- 390 ## row height of panels
    imgH <- 340 ## height of images
    tabV <- "70vh" ## height of tables
    tabH <- 320 ## row height of panels

    gx_infotext <- tspan("The <strong>Differential Expression Analysis</strong> module compares expression between two conditions (i.e. tumor versus control),
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
     frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>", js = FALSE)

    GX.DEFAULTTEST <- "trend.limma"
    GX.DEFAULTTEST <- c("trend.limma", "edger.qlf", "deseq2.wald", "edger.lrt")

    #================================================================================
    #=============================== observers ======================================
    #================================================================================

    my_observers <- list()
    
    my_observers[[1]] <- shiny::observeEvent(input$gx_info,
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
    my_observers[[2]] <- shiny::observe({
      pgx <- pgx
      shiny::req(pgx$X)

      contr <- colnames(pgx$model.parameters$contr.matrix)
      shiny::updateSelectInput(session, "gx_contrast", choices = sort(contr))
      fam <- playbase::pgx.getFamilies(pgx, nmin = 10, extended = FALSE)
      if (length(fam) == 0) {
        fam <- "<all>"
      }
      names(fam) <- sub(".*:", "", fam)
      shiny::updateSelectInput(session, "gx_features", choices = fam)

      ## available statistical methods
      gx.methods <- colnames(pgx$gx.meta$meta[[1]]$fc) ## available
      sel1 <- c(intersect(GX.DEFAULTTEST, gx.methods), gx.methods)
      sel1 <- head(unique(sel1), 3) ## maximum three!!

      shiny::updateCheckboxGroupInput(session, "gx_statmethod",
        choices = sort(gx.methods),
        selected = sel1
      )

      shiny::updateCheckboxInput(session, "gx_grouped", value = (ncol(pgx$X) <= 8))
    })

    my_observers[[3]] <- observeEvent(
      {
        input$show_pv
      },
      {
        if (input$show_pv) {
          shinyalert::shinyalert(
            title = "",
            text = "WARNING: Nominal p-values are NOT corrected for multiple testing. We do not advice their use.",
            type = ""
          )
        }
      }
    )
    
    # observe functions to project DT from invalidating equal row_select
    gsettable_rows_selected <- reactiveVal()

    my_observers[[4]] <- observe({
      gsettable_rows_selected(gsettable$rows_selected())
    })

    genetable_rows_selected <- reactiveVal()

    my_observers[[5]] <- observe({
      genetable_rows_selected(genetable$rows_selected())
    })

    ## add to list global of observers. suspend by default.
    my_observers <- my_observers[!sapply(my_observers,is.null)]
    # lapply( my_observers, function(b) b$suspend() )
    board_observers[[id]] <- my_observers

    ## =========================================================================
    ## ============================= REACTIVES =================================
    ## =========================================================================


    selected_gxmethods <- shiny::reactive({
      req(pgx$X)
      gx.methods0 <- colnames(pgx$gx.meta$meta[[1]]$fc)
      test <- input$gx_statmethod
      test <- intersect(test, gx.methods0)
      test
    })

    shiny::observeEvent(input$show_pv, {
      shiny::updateSelectInput(
        session, "gx_fdr",
        label = if (input$show_pv) "P-value" else "FDR"
      )
    })

    pval_cap <- shiny::reactive({
      pval_cap <- input$pval_cap
      if (pval_cap == "Uncaped") {
        pval_cap <- 1e-999
      } else {
        pval_cap <- as.numeric(input$pval_cap)
      }
      return(pval_cap)
    })

    # functions #########
    comparison <- 1
    testmethods <- c("trend.limma")
    add.pq <- 0
    getDEGtable <- function(pgx, testmethods, comparison, add.pq, lfc, fdr) {
      shiny::req(pgx$X)

      if (is.null(testmethods)) {
        return(NULL)
      }
      if (is.null(comparison)) {
        return(NULL)
      }
      if (length(testmethods) == 0 || all(testmethods == "")) {
        return(NULL)
      }
      if (length(comparison) == 0 || all(comparison == "")) {
        return(NULL)
      }

      ## build meta table
      mx <- pgx$gx.meta$meta[[comparison]]
      if (is.null(mx)) {
        return(NULL)
      }

      ## get rid of AsIs so that behaves as data.frame/list
      ## For some reason, this seems to be really important.
      mm <- colnames(unclass(mx$p))
      testmethods <- intersect(mm, testmethods)
      mx.p <- unclass(mx$p[, testmethods, drop = FALSE])
      mx.q <- unclass(mx$q[, testmethods, drop = FALSE])
      mx.fc <- unclass(mx$fc[, testmethods, drop = FALSE])
      rownames(mx.p) <- rownames(mx)
      rownames(mx.q) <- rownames(mx)
      rownames(mx.fc) <- rownames(mx)

      mx.fc[which(is.infinite(mx.fc) | is.nan(mx.fc))] <- NA
      mx.p[which(is.infinite(mx.p) | is.nan(mx.p))] <- NA
      mx.q[which(is.infinite(mx.q) | is.nan(mx.q))] <- NA

      ## !!!!!!!!!!!!!!!!!!!!!!!! NEED RETHINK !!!!!!!!!!!!!!!!!!!!!!!!
      ## must recompute meta parameters (maxQ method)
      mx$meta.p <- apply(mx.p, 1, max, na.rm = TRUE)
      mx$meta.q <- apply(mx.q, 1, max, na.rm = TRUE)
      mx$meta.fx <- rowMeans(mx.fc, na.rm = TRUE)

      stars.fdr <- fdr
      is.sig <- 1 * (mx.q <= stars.fdr) * (abs(mx$meta.fx) >= lfc)
      stars <- sapply(rowSums(is.sig, na.rm = TRUE), playbase::star.symbols)

      ## recalculate group averages???
      y0 <- pgx$model.parameters$exp.matrix[, comparison]
      names(y0) <- rownames(pgx$model.parameters$exp.matrix)
      AveExpr1 <- rowMeans(pgx$X[rownames(mx), names(which(y0 > 0)), drop = FALSE], na.rm = TRUE)
      AveExpr0 <- rowMeans(pgx$X[rownames(mx), names(which(y0 < 0)), drop = FALSE], na.rm = TRUE)

      logFC <- mx$meta.fx
      ## [hack] adjust averages to match logFC...
      mean0 <- (AveExpr0 + AveExpr1) / 2
      AveExpr1 <- mean0 + logFC / 2
      AveExpr0 <- mean0 - logFC / 2

      if (all(c("map", "chr") %in% colnames(pgx$genes))) {
        colnames(pgx$genes)[which(colnames(pgx$genes) == "chr")] <- "chr0"
        colnames(pgx$genes)[which(colnames(pgx$genes) == "map")] <- "chr"
      }

      aa <- intersect(c("gene_name", "gene_title", "chr"), colnames(pgx$genes))
      gene.annot <- pgx$genes[rownames(mx), aa]
      gene.annot$chr <- sub("_.*", "", gene.annot$chr) ## strip any alt postfix
      res <- data.frame(gene.annot,
        logFC = logFC,
        AveExpr0,
        AveExpr1,
        stars = stars,
        meta.q = mx$meta.q,
        check.names = FALSE
      )

      colnames(mx.q) <- paste0(colnames(mx.q), ".q")
      colnames(mx.p) <- paste0(colnames(mx.p), ".p")

      if (add.pq) {
        mx.q <- mx.q[rownames(mx), , drop = FALSE]
        mx.p <- mx.p[rownames(mx), , drop = FALSE]
        res <- cbind(res, mx.q, meta.p = mx$meta.p, mx.p)
      }

      return(res)
    }

    fullDiffExprTable <- shiny::reactive({
      ## return the full DE table
      if (is.null(pgx)) {
        return(NULL)
      }
      comp <- input$gx_contrast
      tests <- input$gx_statmethod
      fdr <- as.numeric(input$gx_fdr)
      lfc <- as.numeric(input$gx_lfc)

      req(input$gx_contrast, input$gx_statmethod, input$gx_fdr, input$gx_lfc)

      if (is.null(comp)) {
        return(NULL)
      }
      if (is.null(tests)) {
        return(NULL)
      }

      res <- getDEGtable(pgx,
        testmethods = tests, comparison = comp,
        add.pq = TRUE, lfc = lfc, fdr = fdr
      )

      ## Filter on features/genesets
      psel <- rownames(res)
      gx_features <- 1
      gx_features <- input$gx_features
      if (gx_features != "<all>") {
        gset <- unlist(playdata::getGSETS(gx_features))
        psel <- playbase::filterProbes(pgx$genes, gset)
      }
      res <- res[which(rownames(res) %in% psel), , drop = FALSE]
      res <- res[order(-abs(res$logFC)), , drop = FALSE]


      # TODO the annot could be moved to geDEGtable function!!!
      annot <- pgx$genes[rownames(res),
        c("feature", "symbol", "human_ortholog"),
        drop = FALSE
      ]
      # res$gene_name <- NULL, better remove before showing the plot

      res <- cbind(annot, res)

      if (all(is.na(res$human_ortholog))) {
        res$human_ortholog <- NULL
      }
      return(res)
    })

    filteredDiffExprTable <- shiny::reactive({
      ## DE table filtered by FDR and gene family

      shiny::req(pgx$X, input$gx_features, input$gx_fdr, input$gx_lfc)

      comp <- 1
      test <- "trend.limma"
      comp <- input$gx_contrast
      tests <- input$gx_statmethod
      fdr <- as.numeric(input$gx_fdr)
      lfc <- as.numeric(input$gx_lfc)
      res <- fullDiffExprTable()

      if (is.null(res) || nrow(res) == 0) {
        return(NULL)
      }

      fx.col <- grep("mean.diff|logfc|foldchange|meta.fx", colnames(res), ignore.case = TRUE)[1]
      res <- res[order(-abs(res[, fx.col])), ]

      ## just show significant genes
      if (!is.null(input$gx_showall) && !input$gx_showall) {
        n <- length(tests)
        ## sel <- which(res$stars == playbase::star.symbols(n))
        if (input$show_pv) {
          sel <- which(res$meta.p <= fdr & abs(res$logFC) >= lfc)
        } else {
          sel <- which(res$meta.q <= fdr & abs(res$logFC) >= lfc)
        }
        res <- res[sel, , drop = FALSE]
      }

      if (nrow(res) == 0) {
        shiny::validate(shiny::need(nrow(res) > 0, tspan("No genes passed the statistical thresholds. Please update the thresholds on the settings sidebar.", js = FALSE)))
        return(NULL)
      }

      ## limit number of rows???
      return(res)
    })

    genes_in_sel_geneset <- shiny::reactive({
      req(pgx$X, pgx$name)

      if (!is.data.frame(gx_related_genesets()) && gx_related_genesets() == tspan("No geneset for selected gene.", js = FALSE)) {
        sel_gene <- filteredDiffExprTable()$symbol[genetable_rows_selected()]
        return(sel_gene)
      }
      sel_gset <- rownames(gx_related_genesets()[gsettable_rows_selected(), ])
      sel_genes <- pgx$GMT[, sel_gset]
      # return sel_genes that are not zero
      sel_genes <- sel_genes[which(sel_genes > 0)]

      # convert symbol to rownames (module is based on rownames)

      ortholog_genes <- ifelse(is.na(pgx$genes$human_ortholog), pgx$genes$symbol, pgx$genes$human_ortholog)

      matched_rownames <- rownames(pgx$genes[match(names(sel_genes), ortholog_genes), ])
      return(matched_rownames)
    })

    genes_selected <- shiny::reactive({
      shiny::req(input$gx_features)

      df1 <- filteredDiffExprTable()
      df2 <- gx_related_genesets()
      res <- fullDiffExprTable()
      features <- input$gx_features
      fam.genes <- res$symbol
      fdr <- as.numeric(input$gx_fdr) ## ps: can also be P-value if user checks box.
      lfc <- as.numeric(input$gx_lfc)
      if (features != "<all>") {
        gset <- playdata::getGSETS(features)
        fam.genes <- unique(unlist(gset))
      }
      ##      jj <- match(fam.genes, res$symbol)
      jj <- which(res$symbol %in% fam.genes)
      sel.features <- rownames(res)[setdiff(jj, NA)]

      features <- rownames(res)

      qval <- res[, grep("adj.P.Val|meta.q|qval|padj", colnames(res))[1]]
      qval <- pmax(qval, pval_cap())
      pval <- res[, grep("pvalue|meta.p|pval|p_value", colnames(res))[1]]
      pval <- pmax(pval, pval_cap())

      x <- res[, grep("logFC|meta.fx|fc", colnames(res))[1]]
      y <- -log10(qval + pval_cap())
      scaled.x <- scale(x, center = FALSE)
      scaled.y <- scale(y, center = FALSE)

      if (input$show_pv) {
        sig.genes <- features[which(pval <= fdr & abs(x) > lfc)]
      } else {
        sig.genes <- features[which(qval <= fdr & abs(x) > lfc)]
      }
      sel.features <- intersect(sig.genes, sel.features)

      impt <- function(g) {
        j <- match(g, features)
        x1 <- scaled.x[j]
        y1 <- scaled.y[j]
        x <- sign(x1) * (x1**2 + 0.25 * y1**2)
        names(x) <- g
        x
      }

      gene.selected <- !is.null(genetable_rows_selected()) && !is.null(df1)
      gset.selected <- !is.null(gsettable_rows_selected()) && !is.null(df2)

      if (gene.selected && !gset.selected) {
        ## only gene selected: color all genes with same name, label selected
        this.feature <- rownames(df1)[genetable_rows_selected()]
        this.symbol <- res[this.feature, "symbol"]
        sel.features <- rownames(res)[which(res$symbol == this.symbol)]
        ## lab.features <- this.feature
        lab.features <- sel.features
        lab.cex <- 1.3
      } else if (gene.selected && gset.selected) {
        ## gset selected: color all significant gene in geneset, label top 20
        sel.features <- genes_in_sel_geneset()
        lab.features <- c(
          head(sel.features[order(impt(sel.features))], 10),
          head(sel.features[order(-impt(sel.features))], 10)
        )
        lab.cex <- 1
      } else {
        ## no selection: color all significant features, label top 20
        lab.features <- c(
          head(sel.features[order(impt(sel.features))], 10),
          head(sel.features[order(-impt(sel.features))], 10)
        )
        lab.cex <- 1
      }

      res <- list(
        "sel.genes" = sel.features,
        "lab.genes" = lab.features,
        "fc.genes" = features
      )
      return(res)
    })

    # Plotting ###

    # tab differential expression > Plot ####
    expression_plot_volcano_server(
      id = "plots_volcano",
      comp1 = shiny::reactive(input$gx_contrast),
      fdr = shiny::reactive(input$gx_fdr),
      lfc = shiny::reactive(input$gx_lfc),
      show_pv = shiny::reactive(input$show_pv),
      res = fullDiffExprTable,
      genes_selected = genes_selected,
      labeltype = labeltype,
      watermark = WATERMARK,
      pval_cap = pval_cap,
      pgx = pgx
    )

    expression_plot_maplot_server(
      id = "plots_maplot",
      pgx = pgx,
      gx_fdr = reactive(input$gx_fdr),
      gx_contrast = reactive(input$gx_contrast),
      gx_lfc = reactive(input$gx_lfc),
      gx_features = reactive(input$gx_features),
      res = fullDiffExprTable,
      genes_selected = genes_selected,
      labeltype = labeltype,
      watermark = WATERMARK
    )

    expression_plot_barplot_server(
      id = "plots_barplot",
      comp = shiny::reactive(input$gx_contrast),
      pgx = pgx,
      sel = genetable_rows_selected,
      res = filteredDiffExprTable,
      watermark = WATERMARK
    )

    expression_plot_topfoldchange_server(
      id = "plots_topfoldchange",
      comp = shiny::reactive(input$gx_contrast),
      pgx = pgx,
      sel = genetable_rows_selected,
      res = filteredDiffExprTable,
      watermark = WATERMARK
    )

    # tab differential expression > Top genes ####

    getAllContrasts <- shiny::reactive({
      if (is.null(pgx)) {
        return(NULL)
      }

      comp <- names(pgx$gx.meta$meta)
      if (length(comp) == 0) {
        return(NULL)
      }

      tests <- colnames(pgx$gx.meta$meta[[1]]$p)
      tests <- input$gx_statmethod
      if (is.null(tests)) {
        return(NULL)
      }

      i <- 1
      F <- list()
      Q <- list()
      P <- list()
      shiny::withProgress(message = "computing contrasts ...", value = 0, {
        for (i in 1:length(comp)) {
          res <- getDEGtable(pgx,
            testmethods = tests, comparison = comp[i],
            add.pq = TRUE, lfc = 0, fdr = 1
          )
          fc.gene <- rownames(res)
          qv.col <- grep("qval|adj.p|padj|fdr|meta.q", colnames(res), ignore.case = TRUE)[1]
          pv.col <- grep("pval|pvalue|padj|meta.p", colnames(res), ignore.case = TRUE)[1]
          fx.col <- grep("mean.diff|logfc|foldchange|meta.fx", colnames(res), ignore.case = TRUE)[1]
          qval <- res[, qv.col]
          pval <- res[, pv.col]
          fx <- res[, fx.col]
          names(qval) <- names(pval) <- names(fx) <- fc.gene
          F[[i]] <- fx
          Q[[i]] <- qval
          P[[i]] <- pval
          if (!interactive()) shiny::incProgress(1 / length(comp))
        }
      })
      names(Q) <- names(P) <- names(F) <- comp
      F <- do.call(cbind, F)
      P <- do.call(cbind, P)
      Q <- do.call(cbind, Q)
      ct <- list(Q = Q, P = P, F = F)
      return(ct)
    })

    expression_plot_topgenes_server(
      id = "topgenes",
      comp = shiny::reactive(input$gx_contrast),
      pgx = pgx,
      res = filteredDiffExprTable,
      rows_current = genetable$rows_current,
      watermark = WATERMARK
    )

    # tab differential expression > Volcano All ####

    expression_plot_volcanoAll_server(
      id = "volcanoAll",
      pgx = pgx,
      getAllContrasts = getAllContrasts,
      fdr = shiny::reactive(input$gx_fdr),
      lfc = shiny::reactive(input$gx_lfc),
      show_pv = shiny::reactive(input$show_pv),
      genes_selected = genes_selected,
      labeltype = labeltype,
      pval_cap = pval_cap,
      watermark = WATERMARK
    )

    # tab differential expression > Volcano Methods ####

    expression_plot_volcanoMethods_server(
      id = "volcanoMethods",
      pgx = pgx,
      comp = shiny::reactive(input$gx_contrast),
      fdr = shiny::reactive(input$gx_fdr),
      lfc = shiny::reactive(input$gx_lfc),
      show_pv = shiny::reactive(input$show_pv),
      genes_selected = genes_selected,
      labeltype = labeltype,
      pval_cap = pval_cap,
      watermark = WATERMARK
    )

    # tab differential expression > Volcano Methods ####

    # rendering tables ####
    gx_related_genesets <- shiny::reactive({
      res <- filteredDiffExprTable()

      X <- pgx$X
      if (is.null(res) || nrow(res) == 0) {
        return(NULL)
      }
      contr <- input$gx_contrast
      if (is.null(contr)) {
        return(NULL)
      }
      ## get table
      sel.row <- genetable_rows_selected()
      if (is.null(sel.row)) {
        return(NULL)
      }
      gene1 <- res[sel.row, ]

      gmt_gene_mapping <- ifelse(is.na(gene1$human_ortholog), gene1$symbol, gene1$human_ortholog)

      j <- which(rownames(pgx$GMT) == gmt_gene_mapping)
      if (length(j) == 0) {
        return(NULL)
      } else {
        gset <- names(which(pgx$GMT[j, ] != 0))
        gset <- intersect(gset, rownames(pgx$gsetX))
      }

      if (length(gset) == 0) {
        return(tspan("No geneset for selected gene.", js = FALSE))
      }

      fx <- pgx$gset.meta$meta[[contr]]$meta.fx
      names(fx) <- rownames(pgx$gset.meta$meta[[contr]])
      fx <- round(fx[gset], digits = 4)

      # X <- playbase::rename_by(X, pgx$genes, "feature")
      rho <- cor(t(pgx$gsetX[gset, , drop = FALSE]), X[rownames(gene1), ])[, 1]
      rho <- round(rho, digits = 3)
      gset1 <- substring(gset, 1, 60)

      df <- data.frame(geneset = gset1, rho = rho, fx = fx, check.names = FALSE)
      rownames(df) <- gset
      df <- df[order(-abs(df$fx)), ]

      return(df)
    })

    genetable <- expression_table_genetable_server(
      id = "genetable",
      res = filteredDiffExprTable,
      organism = pgx$organism,
      show_pv = shiny::reactive(input$show_pv),
      height = c(tabH - 10, 700),
      scrollY = "200px",
      cont = shiny::reactive(input$gx_contrast)
    )

    gsettable <- expression_table_gsettable_server(
      id = "gsettable",
      gx_related_genesets = gx_related_genesets,
      height = c(tabH - 10, 700),
      width = c("100%", 800),
      scrollY = "200px",
      watermark = WATERMARK
    )

    expression_table_fctable_server(
      id = "fctable",
      pgx = pgx,
      res = filteredDiffExprTable,
      metaFC = metaFC,
      metaQ = metaQ,
      height = c(tabH, 700),
      scrollY = "200px",
      watermark = WATERMARK
    )

    expression_table_FDRtable_server(
      id = "FDRtable",
      pgx = pgx,
      methods = shiny::reactive(input$gx_statmethod),
      height = c(tabH, 700),
      scrollY = "200px",
      watermark = WATERMARK
    )

    # reactive values to return to parent environment  #########
    metaQ <- shiny::reactive({
      req(pgx$X)
      methods <- selected_gxmethods()
      metaQ <- sapply(pgx$gx.meta$meta, function(m) {
        apply(m$q[, methods, drop = FALSE], 1, max, na.rm = TRUE)
      })
      rownames(metaQ) <- rownames(pgx$gx.meta$meta[[1]])
      metaQ
    })

    metaFC <- shiny::reactive({
      req(pgx$X)
      methods <- selected_gxmethods()
      metaFC <- sapply(pgx$gx.meta$meta, function(m) m$meta.fx)
      rownames(metaFC) <- rownames(pgx$gx.meta$meta[[1]])
      metaFC
    })

    outx <- list(
      selected_gxmethods = selected_gxmethods
    )
    return(outx)
  }) ## end of moduleServer
} ## end-of-ExpressionBoard
