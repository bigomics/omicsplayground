##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ExpressionBoard <- function(id, pgx) {
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


    # observe functions ##########


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
      pgx <- pgx
      shiny::req(pgx)

      contr <- colnames(pgx$model.parameters$contr.matrix)
      shiny::updateSelectInput(session, "gx_contrast", choices = sort(contr))
      fam <- playbase::pgx.getFamilies(pgx, nmin = 10, extended = FALSE)
      names(fam) <- sub(".*:","",fam)
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


    # observe functions to project DT from invalidating equal row_select

    gsettable_rows_selected <- reactiveVal()

    observe({
      gsettable_rows_selected(gsettable$rows_selected())
    })

    genetable_rows_selected <- reactiveVal()

    observe({
      ## req(genetable$rows_selected())
      genetable_rows_selected(genetable$rows_selected())
    })

    # reactives ########


    selected_gxmethods <- shiny::reactive({
      req(pgx)
      gx.methods0 <- colnames(pgx$gx.meta$meta[[1]]$fc)
      test <- input$gx_statmethod
      test <- intersect(test, gx.methods0)
      test
    })


    # functions #########
    comparison <- 1
    testmethods <- c("trend.limma")
    add.pq <- 0
    getDEGtable <- function(pgx, testmethods, comparison, add.pq,
                            lfc, fdr) {
      ## if(is.null(pgx)) return(NULL)
      shiny::req(pgx)

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

      ## build meta table
      mx <- pgx$gx.meta$meta[[comparison]]
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
      stars <- sapply(rowSums(is.sig, na.rm = TRUE), playbase::star.symbols)

      ## recalculate group averages???
      y0 <- pgx$model.parameters$exp.matrix[, comparison]
      names(y0) <- rownames(pgx$model.parameters$exp.matrix)
      AveExpr1 <- rowMeans(pgx$X[rownames(mx), names(which(y0 > 0)), drop = FALSE])
      AveExpr0 <- rowMeans(pgx$X[rownames(mx), names(which(y0 < 0)), drop = FALSE])

      ## logFC <- unclass(pgx$gx.meta$meta[[comparison]][,"fc"])[,"trend.limma"]
      ## logFC <- pgx$gx.meta$meta[[comparison]][,"meta.fx"]
      logFC <- mx$meta.fx
      ## logFC <- (AveExpr1 - AveExpr0)  ## override ??? yes: see "contrast in R" Rose Maier 2015...
      ## [hack] adjust averages to match logFC...
      mean0 <- (AveExpr0 + AveExpr1) / 2
      AveExpr1 <- mean0 + logFC / 2
      AveExpr0 <- mean0 - logFC / 2

      ## gene.annot = mx[,grep("^gene|^chr",colnames(mx)),drop=FALSE]
      aa <- intersect(c("gene_name", "gene_title", "chr"), colnames(pgx$genes))
      gene.annot <- pgx$genes[rownames(mx), aa]
      gene.annot$chr <- sub("_.*", "", gene.annot$chr) ## strip any alt postfix
      res <- data.frame(gene.annot,
        logFC = logFC,
        stars = stars, meta.q = mx$meta.q,
        AveExpr0, AveExpr1, check.names = FALSE
      )
      rownames(res) <- rownames(mx)

      if (add.pq) {
        ## add extra columns
        ## res <- cbind( res, q=mx$q, p=mx$p)
        colnames(mx.q) <- paste0("q.", colnames(mx.q))
        res <- cbind(res, mx.q[rownames(mx), , drop = FALSE])
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
        ## gset <- GSETS[[gx_features]]
        gset <- unlist(getGSETS(gx_features))
        psel <- playbase::filterProbes(pgx$genes, gset)
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
      ## if(is.null(pgx)) return(NULL)
      shiny::req(pgx, input$gx_features, input$gx_fdr, input$gx_lfc)

      comp <- 1
      test <- "trend.limma"
      comp <- input$gx_contrast
      tests <- input$gx_statmethod
      fdr <- as.numeric(input$gx_fdr)
      lfc <- as.numeric(input$gx_lfc)

      ## res = getDEGtable(pgx, testmethods="trend.limma", comparison=1,add.pq=FALSE)
      ## res = getDEGtable(pgx, testmethods=tests, comparison=comp,
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
        sel <- which(res$stars == playbase::star.symbols(n))
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

    # Plotting ###

    # tab differential expression > Plot ####
    expression_plot_volcano_server(
      id = "plots_volcano",
      comp1 = shiny::reactive(input$gx_contrast),
      fdr = shiny::reactive(input$gx_fdr),
      lfc = shiny::reactive(input$gx_lfc),
      features = shiny::reactive(input$gx_features),
      res = fullDiffExprTable,
      sel1 = genetable_rows_selected,
      df1 = filteredDiffExprTable,
      sel2 = gsettable_rows_selected,
      df2 = gx_related_genesets
    )

    expression_plot_maplot_server(
      id = "plots_maplot",
      pgx = pgx,
      gx_fdr = reactive(input$gx_fdr),
      gx_contrast = reactive(input$gx_contrast),
      gx_lfc = reactive(input$gx_lfc),
      gx_features = reactive(input$gx_features),
      res = fullDiffExprTable,
      sel1 = genetable_rows_selected,
      df1 = filteredDiffExprTable,
      sel2 = gsettable_rows_selected,
      df2 = gx_related_genesets,
      fam.genes = res$gene_name,
      watermark = FALSE
    )

    expression_plot_barplot_server(
      id = "plots_barplot",
      comp = shiny::reactive(input$gx_contrast),
      pgx = pgx,
      sel = genetable_rows_selected,
      res = filteredDiffExprTable,
      watermark = FALSE
    )

    expression_plot_topfoldchange_server(
      id = "plots_topfoldchange",
      comp = shiny::reactive(input$gx_contrast),
      pgx = pgx,
      sel = genetable_rows_selected,
      res = filteredDiffExprTable,
      watermark = FALSE
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
      ## if(is.null(input$gx_features)) return(NULL)

      ## fdr=1;lfc=0
      ## fdr = as.numeric(input$gx_fdr)
      ## lfc = as.numeric(input$gx_lfc)
      tests <- colnames(pgx$gx.meta$meta[[1]]$p)
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
          res <- getDEGtable(pgx,
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
      comp = shiny::reactive(input$gx_contrast),
      pgx = pgx,
      res = filteredDiffExprTable,
      ii = genetable$rows_current,
      watermark = FALSE
    )

    # tab differential expression > Volcano All ####

    expression_plot_volcanoAll_server(
      id = "volcanoAll",
      pgx = pgx,
      getAllContrasts = getAllContrasts,
      features = shiny::reactive(input$gx_features),
      fdr = shiny::reactive(input$gx_fdr),
      lfc = shiny::reactive(input$gx_lfc),
      watermark = FALSE
    )

    # tab differential expression > Volcano Methods ####

    expression_plot_volcanoMethods_server(
      id = "volcanoMethods",
      pgx = pgx,
      comp = shiny::reactive(input$gx_contrast),
      features = shiny::reactive(input$gx_features),
      fdr = shiny::reactive(input$gx_fdr),
      lfc = shiny::reactive(input$gx_lfc),
      watermark = FALSE
    )

    # tab differential expression > Volcano Methods ####

    # rendering tables ####
    gx_related_genesets <- shiny::reactive({
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
      sel.row <- genetable_rows_selected()
      if (is.null(sel.row)) return(NULL)

      gene0 <- rownames(res)[sel.row]
      gene1 <- toupper(sub(".*:", "", gene0)) ## always uppercase...

      j <- which(toupper(rownames(pgx$GMT)) == gene1)
      gset <- names(which(pgx$GMT[j, ] != 0))
      gset <- intersect(gset, rownames(pgx$gsetX))
      if (length(gset) == 0) {
        return(NULL)
      }

      fx <- pgx$gset.meta$meta[[contr]]$meta.fx
      names(fx) <- rownames(pgx$gset.meta$meta[[contr]])
      fx <- round(fx[gset], digits = 4)

      rho <- cor(t(pgx$gsetX[gset, ]), pgx$X[gene0, ])[, 1]
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
      height = c(tabH - 10, 700),
      scrollY = "200px"
    )

    gsettable <- expression_table_gsettable_server(
      id = "gsettable",
      gx_related_genesets = gx_related_genesets,
      height = c(tabH - 10, 700),
      width = c("100%", 800),
      scrollY = "200px",
      watermark = FALSE
    )

    expression_table_fctable_server(
      id = "fctable",
      pgx = pgx,
      res = filteredDiffExprTable,
      metaFC = metaFC,
      metaQ = metaQ,
      height = c(tabH, 700),
      scrollY = "200px",
      watermark = FALSE
    )

    expression_table_FDRtable_server(
      id = "FDRtable",
      pgx = pgx,
      methods = shiny::reactive(input$gx_statmethod),
      height = c(tabH, 700),
      scrollY = "200px",
      watermark = FALSE
    )

    # reactive values to return to parent environment  #########
    metaQ <- shiny::reactive({
      req(pgx)
      methods <- selected_gxmethods()
      metaQ <- sapply(pgx$gx.meta$meta, function(m)
        apply(m$q[, methods, drop = FALSE], 1, max, na.rm = TRUE))
      rownames(metaQ) <- rownames(pgx$gx.meta$meta[[1]])
      metaQ
    })

    metaFC <- shiny::reactive({
      req(pgx)
      methods <- selected_gxmethods()
      ## metaFC <- sapply(pgx$gx.meta$meta, function(m) rowMeans(m$fc[,methods,drop=FALSE]))
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
