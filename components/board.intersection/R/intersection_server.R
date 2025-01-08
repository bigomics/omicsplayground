##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

IntersectionBoard <- function(
    id,
    pgx,
    selected_gxmethods = reactive(colnames(pgx$gx.meta$meta[[1]]$fc)),
    selected_gsetmethods = reactive(colnames(pgx$gset.meta$meta[[1]]$fc)),
    board_observers = NULL
    ) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 800 # row height of panel

    infotext <-
      tspan("The <strong>Intersection analysis module</strong> enables users to compare multiple contrasts by intersecting the genes of profiles. The main goal is to identify contrasts showing similar profiles.

<br><br>For the selected contrasts, the platform provides volcano plots and pairwise correlation plots between the profiles in the <strong>Pairs</strong> panel. Simultaneously, a Venn diagram with the number of intersecting genes between the profiles is plotted in <strong>Venn diagram</strong> panel. Details of intersecting genes are also reported in an interactive table. A more detailed scatter plot of two profiles is possible under the <strong>Two-pairs</strong> panel. Users can check the pairwise correlations of the contrasts under the <b>Contrast heatmap</b> panel. Alternatively, the <strong>Connectivity Map (CMap)</strong> shows the similarity of the contrasts profiles as a t-SNE plot.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=5' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>

", js = FALSE)

    ## delayed input
    input_comparisons <- shiny::reactive({
      input$comparisons
    }) %>% shiny::debounce(500)

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    my_observers <- list()
    
    my_observers[[1]] <- shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Intersection Analysis Board</strong>"),
        shiny::HTML(infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    ## update choices upon change of data set
    my_observers[[2]] <- shiny::observe({
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

    my_observers[[3]] <- shiny::observeEvent(pgx$X, {
      choices <- c("gene", "geneset")
      choices_names <- c(tspan("gene", js = FALSE), tspan("geneset", js = FALSE))
      names(choices) <- choices_names
      shiny::updateRadioButtons(session, "level", choices = choices)
      shiny::updateTextAreaInput(session, "customlist", placeholder = tspan("Paste your custom gene list", js = FALSE))
    })

    ## update choices upon change of feature level
    ## observeEvent( input$level, {
    my_observers[[4]] <- shiny::observe({
      if (is.null(pgx)) {
        return(NULL)
      }
      shiny::req(input$level)
      if (input$level == "geneset") {
        gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
        ft <- names(gset_collections)
        nn <- sapply(gset_collections, function(x) sum(x %in% rownames(pgx$gsetX)))
        ft <- ft[nn >= 10]
      } else {
        ## gene level
        ft <- playbase::pgx.getFamilies(pgx, nmin = 10, extended = FALSE)
      }
      ft <- sort(ft)
      names(ft) <- sub(".*:", "", ft)
      shiny::updateSelectInput(session, "filter", choices = ft, selected = "<all>")
    })

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Pairwise scatter" = list(disable = NULL),
      ## "Signature clustering" = list(disable = c("comparisons"))
      "Signature clustering" = list(disable = NULL)
    )
    my_observers[[5]] <- shiny::observeEvent(input$tabs1, {
      bigdash::update_tab_elements(input$tabs1, tab_elements)
    })

    ## add to list global of observers. suspend by default.
    my_observers <- my_observers[!sapply(my_observers,is.null)]
    # lapply( my_observers, function(b) b$suspend() )
    if(!is.null(board_observers)) board_observers[[id]] <- my_observers
    
    ## ================================================================================
    ## ========================= REACTIVE FUNCTIONS ===================================
    ## ================================================================================

    getFoldChangeMatrix <- shiny::reactive({
      ## Get full foldchange matrix from pgx object.
      fc0 <- NULL
      qv0 <- NULL
      alertDataLoaded(session, pgx)
      shiny::req(pgx$X)

      sel <- names(pgx$gset.meta$meta)
      if (input$level == "geneset") {
        gsetmethods <- c("gsva", "camera", "fgsea")
        gsetmethods <- selected_gsetmethods()
        if (length(gsetmethods) < 1 || gsetmethods[1] == "") {
          return(NULL)
        }

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
          gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
          gsets <- unique(unlist(gset_collections[input$filter]))
        }
        gsets <- intersect(gsets, rownames(fc0))
        fc1 <- fc0[gsets, , drop = FALSE]
        qv1 <- qv0[gsets, , drop = FALSE]
      } else {
        ## Gene
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

        ## filter with active filter
        sel.probes <- rownames(fc0) ## default to all probes
        if (input$filter == "<custom>") {
          genes <- strsplit(input$customlist, split = "[, ;]")[[1]]
          if (length(genes) > 0) {
            sel.probes <- playbase::filterProbes(pgx$genes, genes)
          }
        } else if (input$filter != "<all>") {
          gset.genes <- unlist(playdata::getGSETS(input$filter))
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
      shiny::req(res)

      ## match with selected/active contrasts
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
      shiny::req(pgx$X)
      progress <- shiny::Progress$new()
      on.exit(progress$close())

      ## ------------ UMAP clustering (genes) -----------------
      progress$inc(0.33, tspan("calculating UMAP for genes...", js = FALSE))
      if ("cluster.genes" %in% names(pgx)) {
        pos <- pgx$cluster.genes$pos[["umap2d"]]
      } else {
        X1 <- pgx$X
        X1 <- (X1 - rowMeans(X1, na.rm = TRUE)) / mean(apply(X1, 1, sd, na.rm = TRUE))
        pos <- playbase::pgx.clusterBigMatrix(
          t(X1),
          methods = "umap", dims = 2, reduce.sd = -1
        )[[1]]
        pos <- playbase::pos.compact(pos)
      }

      ## ------------ UMAP clustering (genesets) -----------------
      progress$inc(0.33, tspan("calculating UMAP for genesets...", js = FALSE))
      if ("cluster.gsets" %in% names(pgx)) {
        gsea.pos <- pgx$cluster.gsets$pos[["umap2d"]]
      } else {
        X2 <- pgx$gsetX
        X2 <- (X2 - rowMeans(X2, na.rm = TRUE)) / mean(apply(X2, 1, sd, na.rm = TRUE))
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
      f.score <- F$fc * (1 - F$qv)**4 ## q-weighted FC
      g.score <- G$fc * (1 - G$qv)**4

      ii <- intersect(rownames(pos), rownames(f.score))
      sig <- f.score[ii, , drop = FALSE]
      pos <- pos[ii, ]
      ii <- order(-rowMeans(sig, na.rm = TRUE))
      sig <- sig[ii, , drop = FALSE]
      pos <- pos[ii, ]

      ii <- intersect(rownames(gsea.pos), rownames(g.score))
      gsea <- g.score[ii, , drop = FALSE]
      gsea.pos <- gsea.pos[ii, ]
      ii <- order(-rowMeans(gsea, na.rm = TRUE))
      gsea <- gsea[ii, , drop = FALSE]
      gsea.pos <- gsea.pos[ii, ]

      out <- list(sig = sig, pos = pos, gsea = gsea, gsea.pos = gsea.pos)

      progress$close()
      out
    })

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
      tt <- playdata::GENE_TITLE[rownames(W)]
      tt <- substring(tt, 1, 80)
      df <- data.frame(gene = rownames(W), title = tt, W, check.names = FALSE)
      sel1 <- ctGseaTable_module$rows_selected()
      if (length(sel1) > 0) {
        gset <- rownames(out$gsea)[sel1]
        gset.genes <- unlist(playdata::getGSETS(gset))
        gg <- intersect(rownames(df), gset.genes)
        df <- df[gg, , drop = FALSE]
      }
      df
    })



    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    ## first tab ---------------------------------------

    intersection_plot_venn_diagram_server(
      "venndiagram",
      pgx = pgx,
      level = shiny::reactive(input$level),
      input_comparisons = input_comparisons,
      getFoldChangeMatrix = getFoldChangeMatrix,
      watermark = WATERMARK
    )

    intersection_scatterplot_pairs_server(
      "scatterplot",
      getActiveFoldChangeMatrix = getActiveFoldChangeMatrix,
      level = shiny::reactive(input$level),
      pgx = pgx,
      watermark = WATERMARK
    )

    ## second tab ---------------------------------------

    foldchange_heatmap_server(
      "FoldchangeHeatmap",
      getFoldChangeMatrix = getFoldChangeMatrix,
      getActiveFoldChangeMatrix = getActiveFoldChangeMatrix,
      pgx = pgx,
      input_comparisons = input_comparisons,
      level = input$level,
      watermark = WATERMARK
    )

    contrast_correlation_server(
      "ctcorrplot",
      getFoldChangeMatrix = getFoldChangeMatrix,
      pgx = pgx,
      input_comparisons = input_comparisons
    )
  })
} ## end-of-Board
