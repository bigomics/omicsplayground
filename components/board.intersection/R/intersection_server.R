##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

IntersectionBoard <- function(id, pgx, selected_gxmethods, selected_gsetmethods) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 800 # row height of panel

    infotext <-
      "The <strong>Intersection analysis module</strong> enables users to compare multiple contrasts by intersecting the genes of profiles. The main goal is to identify contrasts showing similar profiles.

<br><br>For the selected contrasts, the platform provides volcano plots and pairwise correlation plots between the profiles in the <strong>Pairs</strong> panel. Simultaneously, a Venn diagram with the number of intersecting genes between the profiles is plotted in <strong>Venn diagram</strong> panel. Details of intersecting genes are also reported in an interactive table. A more detailed scatter plot of two profiles is possible under the <strong>Two-pairs</strong> panel. Users can check the pairwise correlations of the contrasts under the <b>Contrast heatmap</b> panel. Alternatively, the <strong>Connectivity Map (CMap)</strong> shows the similarity of the contrasts profiles as a t-SNE plot.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=5' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>

"

    ## delayed input
    input_comparisons <- shiny::reactive({
      input$comparisons
    }) %>% shiny::debounce(500)

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Intersection Analysis Board</strong>"),
        shiny::HTML(infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    ## update choices upon change of data set
    shiny::observe({
      ## req(pgx)
      if (is.null(pgx)) {
        return(NULL)
      }
      comparisons <- colnames(pgx$model.parameters$contr.matrix)
      comparisons <- sort(comparisons)
      shiny::updateSelectInput(session, "comparisons",
        choices = comparisons,
        selected = head(comparisons, 3)
      )
    })

    ## update choices upon change of feature level
    ## observeEvent( input$level, {
    shiny::observe({
      ## shiny::req(pgx,input$level)
      if (is.null(pgx)) {
        return(NULL)
      }
      shiny::req(input$level)
      ## flt.choices = names(pgx$families)
      if (input$level == "geneset") {
        ft <- names(COLLECTIONS)
        nn <- sapply(COLLECTIONS, function(x) sum(x %in% rownames(pgx$gsetX)))
        ft <- ft[nn >= 10]
      } else {
        ## gene level
        ft <- playbase::pgx.getFamilies(pgx, nmin = 10, extended = FALSE)
      }
      ft <- sort(ft)
      ## if(input$level=="gene") ft = sort(c("<custom>",ft))
      ## ft = sort(c("<custom>",ft))
      names(ft) <- sub(".*:","",ft)
      shiny::updateSelectInput(session, "filter", choices = ft, selected = "<all>")
    })

    ## shiny::observe({
    ##     splom.sel <- plotly::event_data("plotly_selected", source="splom")
    ##     sel.keys <- as.character(splom.sel$key)
    ##     if(0 && length(sel.keys)>0) {
    ##         shiny::updateSelectInput(session, "filter", selected="<custom>")
    ##         sel.keys = paste(sel.keys, collapse=" ")
    ##         shiny::updateTextAreaInput(session, "customlist", value=sel.keys)
    ##     }
    ## })


    ## ================================================================================
    ## ========================= REACTIVE FUNCTIONS ===================================
    ## ================================================================================

    getFoldChangeMatrix <- shiny::reactive({
      ##
      ## Get full foldchange matrix from pgx object.
      ##
      ##
      ##
      fc0 <- NULL
      qv0 <- NULL
      alertDataLoaded(session, pgx)
      shiny::req(pgx)

      sel <- names(pgx$gset.meta$meta)
      ## sel = input_comparisons()
      ## sel = intersect(sel, names(pgx$gset.meta$meta))
      ## if(length(sel)==0) return(NULL)

      if (input$level == "geneset") {
        gsetmethods <- c("gsva", "camera", "fgsea")
        gsetmethods <- selected_gsetmethods()
        if (length(gsetmethods) < 1 || gsetmethods[1] == "") {
          return(NULL)
        }

        ## fc0 = sapply(pgx$gset.meta$meta[sel], function(x)
        ##    rowMeans(unclass(x$fc)[,gsetmethods,drop=FALSE]))
        fc0 <- sapply(pgx$gset.meta$meta[sel], function(x) x$meta.fx)
        rownames(fc0) <- rownames(pgx$gset.meta$meta[[1]])
        qv0 <- sapply(pgx$gset.meta$meta[sel], function(x) {
          apply(unclass(x$q)[, gsetmethods, drop = FALSE], 1, max)
        })
        rownames(qv0) <- rownames(pgx$gset.meta$meta[[1]])

        ## apply user selected filter
        gsets <- rownames(fc0)
        if (input$filter == "<custom>") {
          gsets <- strsplit(input$customlist, split = "[, ;]")[[1]]
          if (length(gsets) > 0) {
            gsets <- intersect(rownames(pgx$gsetX), gsets)
          }
        } else if (input$filter != "<all>") {
          gsets <- unique(unlist(COLLECTIONS[input$filter]))
        }
        gsets <- intersect(gsets, rownames(fc0))
        fc1 <- fc0[gsets, , drop = FALSE]
        qv1 <- qv0[gsets, , drop = FALSE]
      } else {
        ## Gene
        ##
        gxmethods <- "trend.limma"
        gxmethods <- c("trend.limma", "edger.qlf", "deseq2.wald")
        gxmethods <- selected_gxmethods() ## reactive object from EXPRESSION section

        mq1 <- pgx$gx.meta$meta[[1]]$meta.q

        if (length(gxmethods) < 1 || gxmethods[1] == "") {
          return(NULL)
        }

        fc0 <- sapply(pgx$gx.meta$meta[sel], function(x) x$meta.fx)
        rownames(fc0) <- rownames(pgx$gx.meta$meta[[1]])
        qv0 <- sapply(pgx$gx.meta$meta[sel], function(x) {
          apply(unclass(x$q)[, gxmethods, drop = FALSE], 1, max)
        })
        rownames(qv0) <- rownames(pgx$gx.meta$meta[[1]])
        dim(fc0)
        dim(qv0)

        ## filter with active filter
        sel.probes <- rownames(fc0) ## default to all probes
        if (input$filter == "<custom>") {
          genes <- strsplit(input$customlist, split = "[, ;]")[[1]]
          if (length(genes) > 0) {
            sel.probes <- playbase::filterProbes(pgx$genes, genes)
          }
        } else if (input$filter != "<all>") {
          ## gset <- GSETS[[input$filter]]
          gset.genes <- unlist(getGSETS(input$filter))
          sel.probes <- playbase::filterProbes(pgx$genes, gset.genes)
        }
        sel.probes <- intersect(sel.probes, rownames(fc0))
        fc1 <- fc0[sel.probes, , drop = FALSE]
        qv1 <- qv0[sel.probes, , drop = FALSE]
      }
      fc1 <- fc1[, !duplicated(colnames(fc1)), drop = FALSE]
      qv1 <- qv1[, !duplicated(colnames(qv1)), drop = FALSE]

      res <- list(fc = fc1, qv = qv1, fc.full = fc0, qv.full = qv0)
      return(res)
    })


    getActiveFoldChangeMatrix <- shiny::reactive({
      res <- getFoldChangeMatrix()
      ## if(is.null(res)) return(NULL)
      shiny::req(res)

      ## match with selected/active contrasts
      ## comp = head(colnames(res$fc),3)
      comp <- input_comparisons()
      kk <- match(comp, colnames(res$fc))
      if (length(kk) == 0) {
        return(NULL)
      }
      if (length(kk) == 1) kk <- c(kk, kk)
      res$fc <- res$fc[, kk, drop = FALSE]
      res$qv <- res$qv[, kk, drop = FALSE]
      res$fc.full <- res$fc.full[, kk, drop = FALSE]
      res$qv.full <- res$qv.full[, kk, drop = FALSE]

      return(res)
    })

    getCurrentSig <- shiny::reactive({
      ## Switch between FC profile or NMF vectors
      ##
      ##
      shiny::req(pgx)
      progress <- shiny::Progress$new()
      on.exit(progress$close())

      ## ------------ UMAP clustering (genes) -----------------
      progress$inc(0.33, "calculating UMAP for genes...")
      if ("cluster.genes" %in% names(pgx)) {
        pos <- pgx$cluster.genes$pos[["umap2d"]]
      } else {
        X1 <- pgx$X
        X1 <- (X1 - rowMeans(X1)) / mean(apply(X1, 1, sd, na.rm = TRUE))
        pos <- playbase::pgx.clusterBigMatrix(
          t(X1),
          methods = "umap", dims = 2, reduce.sd = -1
        )[[1]]
        pos <- playbase::pos.compact(pos)
      }

      ## ------------ UMAP clustering (genesets) -----------------
      progress$inc(0.33, "calculating UMAP for genesets...")
      if ("cluster.gsets" %in% names(pgx)) {
        gsea.pos <- pgx$cluster.gsets$pos[["umap2d"]]
      } else {
        X2 <- pgx$gsetX
        X2 <- (X2 - rowMeans(X2)) / mean(apply(X2, 1, sd, na.rm = TRUE))
        gsea.pos <- playbase::pgx.clusterBigMatrix(
          t(X2),
          methods = "umap", dims = 2, reduce.sd = -1
        )[[1]]
        gsea.pos <- playbase::pos.compact(gsea.pos)
        dim(gsea.pos)
      }

      ## ------------ get signature matrices -----------------
      F <- playbase::pgx.getMetaMatrix(pgx, level = "gene")
      G <- playbase::pgx.getMetaMatrix(pgx, level = "geneset")
      ## f.score <- F$fc * -log10(F$qv)
      ## g.score <- G$fc * -log10(G$qv)
      f.score <- F$fc * (1 - F$qv)**4 ## q-weighted FC
      g.score <- G$fc * (1 - G$qv)**4

      ii <- intersect(rownames(pos), rownames(f.score))
      sig <- f.score[ii, , drop = FALSE]
      pos <- pos[ii, ]
      ii <- order(-rowMeans(sig))
      sig <- sig[ii, , drop = FALSE]
      pos <- pos[ii, ]

      ii <- intersect(rownames(gsea.pos), rownames(g.score))
      gsea <- g.score[ii, , drop = FALSE]
      gsea.pos <- gsea.pos[ii, ]
      ii <- order(-rowMeans(gsea))
      gsea <- gsea[ii, , drop = FALSE]
      gsea.pos <- gsea.pos[ii, ]

      out <- list(sig = sig, pos = pos, gsea = gsea, gsea.pos = gsea.pos)

      progress$close()
      out
    })

    ## ================================================================================
    ## Gene/geneset UMAP plots
    ## ================================================================================

    ctGeneUMAP.RENDER <- shiny::reactive({
      out <- getCurrentSig()
      pos <- out$pos

      W <- out$sig
      sel0 <- input_comparisons()
      shiny::req(sel0)
      if (!all(sel0 %in% colnames(W))) {
        return(NULL)
      }
      W <- W[, sel0, drop = FALSE]
      dim(W)

      sel1 <- sel2 <- 1
      sel1 <- ctGeneTable_module$rows_selected()
      sel2 <- ctGseaTable_module$rows_selected()

      is.single <- (NCOL(W) == 1)
      hilight <- NULL
      hilight2 <- NULL
      cexlab <- 1.2
      cex <- ifelse(ncol(W) > 9, 0.5, 0.8)
      gset <- NULL

      if (length(sel1) > 0) {
        df <- getGeneTable()
        sel.gene <- rownames(df)[1]
        sel.gene <- rownames(df)[sel1]
        hilight2 <- sel.gene
        hilight <- c(hilight, sel.gene)
        cexlab <- 1.8
      }
      if (length(sel2) > 0) {
        gse <- out$gsea
        gset <- rownames(gse)[sel2]
        gset.genes <- unlist(getGSETS(gset))
        hilight <- c(hilight, gset.genes)
        hilight2 <- c(hilight2, hilight)
      }

      no.hilight <- (length(hilight) == 0)
      is.nmf <- all(grepl("NMF", colnames(W)))
      ntop <- as.integer(input$ntop_genes)
      lab.pos <- NULL

      ## Spectral NMF vectors
      nc <- ceiling(sqrt(ncol(W)))
      nr <- ceiling(ncol(W) / nc)
      par(mfrow = c(nr, nc))
      par(mar = c(2.7, 2.8, 0.7, 0.2), mgp = c(1.4, 0.5, 0), cex.axis = 0.9, cex.lab = 0.9)
      i <- 1
      for (i in 1:ncol(W)) {
        jj <- match(rownames(pos), rownames(W))
        fc <- W[jj, i]
        fc[is.na(fc)] <- 0

        if (no.hilight) {
          ## no selected gene or geneset
          sorted.fc <- sort(fc)
          hilight <- c(head(names(sorted.fc), ntop / 2), tail(names(sorted.fc), ntop / 2))
          hilight2 <- hilight
        }
        sorted.h2 <- hilight2[order(-fc[hilight2]**2)]
        ## hilight2 = c(head(sorted.h2,ntop),tail(sorted.h2,ntop))
        hilight2 <- head(sorted.h2, ntop)
        opacity <- ifelse(length(hilight) > 0 || length(hilight2) > 0, 0.15, 1)
        zsym <- ifelse(is.nmf, FALSE, TRUE)
        if (min(fc, na.rm = TRUE) >= 0) zsym <- FALSE

        plt <- playbase::pgx.scatterPlotXY.BASE(
          pos,
          var = fc,
          lab.pos = lab.pos,
          softmax = 1,
          zsym = zsym,
          cex.lab = cexlab,
          opacity = opacity,
          cex = cex,
          cex.legend = 0.9,
          hilight.cex = 1.3 * cex,
          hilight.col = NULL,
          hilight.lwd = 0.5,
          hilight = hilight,
          hilight2 = hilight2,
          set.par = FALSE,
          lab.xpd = FALSE, ## clip repel
          dlim = c(0.05, 0.05),
          bty = "n",
          xlab = "UMAP-y  (genes)",
          ylab = "UMAP-y  (genes)",
          legend.pos = "bottomleft"
        )
        title(colnames(W)[i], cex.main = 1.1, line = -0.5)
        if (!is.null(gset)) {
          tt <- substring(tolower(gset), 1, 80)
          title(tt, cex.main = 0.8, line = -1.3, font.main = 1)
          lab.pos <- plt$lab.pos
        }
      }
    })

    ctGeneUMAP_info <- "<b>CORSA module analysis.</b> Functional analysis of NMF modules."

    ctGeneUMAP_opts <- shiny::tagList(
      shiny::radioButtons(ns("ntop_genes"), "Gene labels:",
        choices = c(0, 10, 30, 100), selected = 10, inline = TRUE
      )
    )

    shiny::callModule(
      plotModule,
      id = "ctGeneUMAP", ## ns=ns,
      title = "GENE SIGNATURE", label = "a",
      func = ctGeneUMAP.RENDER,
      func2 = ctGeneUMAP.RENDER,
      download.fmt = c("png", "pdf"),
      options = ctGeneUMAP_opts,
      info.text = ctGeneUMAP_info,
      height = c(480, 780), width = c("auto", 1200),
      pdf.height = 8, pdf.width = 12,
      res = c(75, 95),
      add.watermark = WATERMARK
    )

    ## -----------------------------------------------------------------
    ## ----------------------  NMF gsea --------------------------------
    ## -----------------------------------------------------------------

    ctGseaUMAP.RENDER <- shiny::reactive({
      out <- getCurrentSig()
      gse <- out$gsea

      sel0 <- 1:4
      sel0 <- input_comparisons()
      shiny::req(sel0)
      if (!all(sel0 %in% colnames(gse))) {
        return(NULL)
      }
      gse <- gse[, sel0, drop = FALSE]
      gse0 <- gse

      sel2 <- ctGseaTable_module$rows_selected()
      gset <- NULL
      if (length(sel2) > 0) {
        gset <- rownames(gse)[sel2]
      }

      ii <- grep("HALLMARK", rownames(gse))
      ii <- ctGseaTable_module$rows_all()
      shiny::req(ii)
      gse <- gse[ii, , drop = FALSE]

      ## multiple group layout
      gse.scores <- lapply(1:ncol(gse), function(i) gse[, i])
      names(gse.scores) <- colnames(gse)

      ## layout
      ngse <- length(gse.scores)
      nc <- ceiling(sqrt(ngse))
      nr <- ceiling(ngse / nc)
      cex1 <- 1.15
      if (nc == 2 && nr > 1) cex1 <- 0.9
      if (nr == 1) cex1 <- 0.8
      par(mfrow = c(nr, nc))

      gstype <- input$gstype
      if (gstype == "bar") {
        par(oma = c(0, 0, 0, 0))
        par(mar = c(2.8, 2, 1.8, 3), mgp = c(1.6, 0.6, 0))
        ntop <- round(28 / nr)
        ## maximum 24 terms
        gse.scores <- lapply(gse.scores, function(x) {
          head(sort(x, decreasing = TRUE), ntop)
        })
        xlim <- c(0, max(abs(unlist(gse.scores))))

        i <- 1
        for (i in 1:ngse) {
          playbase::gsea.barplot(
            gse.scores[[i]],
            names = names(gse.scores[[i]]),
            n = ntop, xlim = xlim,
            main = "", cex.text = cex1
          )
          title(names(gse.scores)[i], cex.main = 1.1, line = +0.3)
        }
      }

      if (gstype == "umap") {
        par(mar = c(2.7, 2.8, 0.7, 0.2), mgp = c(1.4, 0.5, 0), cex.axis = 0.9, cex.lab = 0.9)
        pos <- out$gsea.pos
        ntop <- as.integer(input$ntop_gsets)
        cex <- ifelse(ngse > 9, 0.5, 0.8)

        i <- 1
        for (i in 1:length(gse.scores)) {
          k <- names(gse.scores)[i]
          if (ncol(gse) == 1) {
            var <- gse0[, 1]
          } else {
            var <- gse0[, k]
          }
          if (!is.null(gset)) {
            hmarks <- gset
          } else {
            var1 <- gse.scores[[i]]
            ss <- names(sort(var))
            ss <- intersect(ss, names(var1))
            hmarks <- c(head(ss, ntop / 2), tail(ss, ntop / 2))
          }
          var <- var[match(rownames(pos), names(var))]
          opacity <- ifelse(length(hmarks) > 0, 0.15, 1)
          zsym <- ifelse(min(var, na.rm = TRUE) >= 0, FALSE, TRUE)

          playbase::pgx.scatterPlotXY(
            pos,
            var = var,
            zsym = zsym, set.par = FALSE, softmax = 1,
            cex = cex, cex.legend = 0.9, cex.lab = 1.2, bty = "n",
            plotlib = "base", col = "grey70", dlim = c(0.2, 0.08),
            hilight = hmarks, hilight.col = NULL, opacity = opacity,
            xlab = "UMAP-y  (genesets)", ylab = "UMAP-y  (genesets)",
            hilight.lwd = 0.5, hilight.cex = 1.3
          )
          title(k, cex.main = 1.1, line = -0.5)
        }
      }
    })

    ctGseaUMAP.info <- "<b>Functional signature.</b> Functional enrichment signatures for selected comparisons are shown as barplots or UMAP plots."

    ctGseaUMAP.opts <- shiny::tagList(
      shiny::radioButtons(ns("gstype"), "Geneset plot type:",
        choices = c("bar", "umap"),
        selected = "bar", inline = TRUE
      ),
      shiny::radioButtons(ns("ntop_gsets"), "Geneset labels:",
        choices = c(0, 6, 10, 30, 100), selected = 6, inline = TRUE
      )
    )

    shiny::callModule(
      plotModule,
      id = "ctGseaUMAP", ## ns=ns,
      title = "FUNCTIONAL SIGNATURE", label = "b",
      func = ctGseaUMAP.RENDER,
      func2 = ctGseaUMAP.RENDER,
      download.fmt = c("png", "pdf"),
      options = ctGseaUMAP.opts,
      info.text = ctGseaUMAP.info,
      height = c(480, 780), width = c("auto", 1200),
      pdf.height = 8, pdf.width = 12,
      res = c(75, 95),
      add.watermark = WATERMARK
    )

    ## -------------------------------------------
    ## --------------gene table ------------------
    ## -------------------------------------------

    getGeneTable <- shiny::reactive({
      out <- getCurrentSig()

      W <- out$sig
      sel0 <- 1:ncol(W)
      sel0 <- input_comparisons()
      shiny::req(sel0)
      if (length(sel0) == 0) {
        return(NULL)
      }
      if (!all(sel0 %in% colnames(W))) {
        return(NULL)
      }

      ## only genes
      W <- W[rownames(W) %in% rownames(pgx$X), , drop = FALSE]
      W <- W[, sel0, drop = FALSE]

      tt <- NA
      tt <- GENE.TITLE[rownames(W)]
      tt <- substring(tt, 1, 80)
      df <- data.frame(gene = rownames(W), title = tt, W, check.names = FALSE)
      sel1 <- ctGseaTable_module$rows_selected()
      if (length(sel1) > 0) {
        gset <- rownames(out$gsea)[sel1]
        gset.genes <- unlist(getGSETS(gset))
        gg <- intersect(rownames(df), gset.genes)
        df <- df[gg, , drop = FALSE]
      }
      df
    })

    ctGeneTable.RENDER <- shiny::reactive({
      df <- getGeneTable()
      shiny::req(df)
      numeric.cols <- colnames(df)[3:ncol(df)]

      DT::datatable(
        df,
        rownames = FALSE, ## escape = c(-1,-2),
        ## filter = 'top',
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip", ## buttons = c('copy','csv','pdf'),
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE, ## scrollY = TRUE,
          ## scrollY = 170,
          scrollY = "70vh",
          scroller = TRUE, deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    ctGeneTable.info <- "Genes."
    ctGeneTable.caption <- "<b>Module gene table.</b> <b>(a)</b> ..."
    ctGeneTable.opts <- shiny::tagList(
      ## selectInput(ns('wgcna_ctGeneTable_selmodule'),'module:', choices=NULL)
    )

    ctGeneTable_module <- shiny::callModule(
      tableModule,
      id = "ctGeneTable",
      title = "GENE TABLE", label = "I",
      ## caption = wgcna_ctGeneTable_caption,
      func = ctGeneTable.RENDER, ## ns=ns,
      options = ctGeneTable.opts,
      info.text = ctGeneTable.info,
      height = c(225, 750), width = c("100%", 1500)
    )

    ## -------------------------------------------
    ## -------------- GSEA table --------------
    ## -------------------------------------------

    ctGseaTable.RENDER <- shiny::reactive({
      out <- getCurrentSig()
      df <- out$gsea
      sel0 <- input_comparisons()
      shiny::req(sel0)
      if (!all(sel0 %in% colnames(df))) {
        return(NULL)
      }
      df <- df[, sel0, drop = FALSE]
      numeric.cols <- colnames(df)
      gset <- substring(rownames(df), 1, 80)
      df <- data.frame(geneset = gset, df, check.names = FALSE)

      DT::datatable(
        df,
        rownames = FALSE, ## escape = c(-1,-2),
        ## filter = 'top',
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip", ## buttons = c('copy','csv','pdf'),
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE, ## scrollY = TRUE,
          ## scrollY = 170,
          scrollY = "70vh",
          scroller = TRUE, deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    ctGseaTable.info <- "NMF GSEA."
    ctGseaTable.caption <- "<b>Module enrichment table.</b> <b>(a)</b> ..."
    ctGseaTable.opts <- shiny::tagList(
      ## selectInput(ns('wgcna_ctGseaTable_selmodule'),'module:', choices=NULL)
    )

    ctGseaTable_module <- shiny::callModule(
      tableModule,
      id = "ctGseaTable",
      title = "ENRICHMENT TABLE", label = "II",
      ## caption = ctGseaTable_caption,
      func = ctGseaTable.RENDER, ## ns=ns,
      options = ctGseaTable.opts,
      info.text = ctGseaTable.info,
      height = c(225, 750), width = c("100%", 1500)
    )

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    ## first tab ---------------------------------------

    WATERMARK <- FALSE

    intersection_plot_venn_diagram_server(
      "venndiagram",
      pgx           = pgx,
      level               = input$level,
      input_comparisons   = input_comparisons,
      getFoldChangeMatrix = getFoldChangeMatrix,
      watermark           = WATERMARK
    )

    intersection_scatterplot_pairs_server(
      "scatterplot",
      getActiveFoldChangeMatrix = getActiveFoldChangeMatrix,
      level                     = input$level,
      pgx                 = pgx,
      watermark                 = WATERMARK
    )

    ## second tab ---------------------------------------

    foldchange_heatmap_server(
      "FoldchangeHeatmap",
      getFoldChangeMatrix       = getFoldChangeMatrix,
      getActiveFoldChangeMatrix = getActiveFoldChangeMatrix,
      pgx                 = pgx,
      level                     = input$level,
      watermark                 = WATERMARK
    )

    contrast_correlation_server(
      "ctcorrplot",
      getFoldChangeMatrix = getFoldChangeMatrix,
      pgx           = pgx,
      input_comparisons   = input_comparisons
    )
  })
} ## end-of-Board
