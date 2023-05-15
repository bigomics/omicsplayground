##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

CompareBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 770 # row height of panel
    tabH <- "70vh"

    infotext <-
      "The <strong>Compare Datasets</strong> module enables users to compare their dataset to other datasets.
         This module allows side-by-side comparison of volcano, scatter or gene t-SNE plots.
         It provides pairwise correlation plots and/or enrichment plots with signatures from other data sets.
        <br><br><br><br>
        <center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=5'
        frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>"

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Compare Experiments</strong>"),
        shiny::HTML(infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    shiny::observe({

      comparisons1 <- names(pgx$gx.meta$meta)
      sel1 <- comparisons1[1]
      shiny::updateSelectInput(session, "contrast1", choices = comparisons1, selected = sel1)

      pgx.files <- sort(dir(file.path(OPG, "data"), pattern = "pgx$"))
      shiny::updateSelectInput(session, "dataset2", choices = c("<this>", pgx.files))
    })

    shiny::observe({
      df <- getOmicsScoreTable()
      if (is.null(df)) {
        return(NULL)
      }
      ntop <- as.integer(input$ntop)
      higenes <- rownames(df)[order(df$score, decreasing = TRUE)]
      higenes <- head(higenes, ntop)
      higenes <- paste(higenes, collapse = " ")
      shiny::updateTextAreaInput(session, "genelist", value = higenes)
    })

    ## ================================================================================
    ## ========================= REACTIVE FUNCTIONS ===================================
    ## ================================================================================

    cum_fc <- shiny::reactive({
      pgx1 <- pgx
      pgx2 <- dataset2()

      ct1 <- head(names(pgx1$gx.meta$meta), 2)
      ct2 <- head(names(pgx2$gx.meta$meta), 2)
      ct1 <- input.contrast1()
      ct2 <- input.contrast2()
      shiny::req(ct1)
      shiny::req(ct2)
      if (!all(ct1 %in% names(pgx1$gx.meta$meta))) {
        return(NULL)
      }
      if (!all(ct2 %in% names(pgx2$gx.meta$meta))) {
        return(NULL)
      }

      F1 <- playbase::pgx.getMetaMatrix(pgx1)$fc[, ct1, drop = FALSE]
      F2 <- playbase::pgx.getMetaMatrix(pgx2)$fc[, ct2, drop = FALSE]

      gg <- intersect(toupper(rownames(F1)), toupper(rownames(F2)))
      g1 <- rownames(F1)[match(gg, toupper(rownames(F1)))]
      g2 <- rownames(F2)[match(gg, toupper(rownames(F2)))]
      F1 <- F1[g1, , drop = FALSE]
      F2 <- F2[g2, , drop = FALSE]
      colnames(F1) <- paste0("1:", colnames(F1))
      colnames(F2) <- paste0("2:", colnames(F2))

      return(cbind(F1, F2))
    })

    dataset2 <- shiny::reactive({
      shiny::req(input$dataset2)
      if (input$dataset2 == "<this>") {
        pgx <- pgx
      } else {
        ##load(file.path(OPG, "data", input$dataset2))
        file2 <- file.path(OPG, "data", input$dataset2)
        pgx <- local(get(load(file2, verbose=0)))
      }
      comparisons2 <- names(pgx$gx.meta$meta)
      sel2 <- tail(head(comparisons2, 2), 1)
      shiny::updateSelectInput(session, "contrast2", choices = comparisons2, selected = sel2)
      pgx
    })

    getOmicsScoreTable <- shiny::reactive({
      pgx1 <- pgx
      pgx2 <- dataset2()
      shiny::req(pgx1)
      shiny::req(pgx2)

      ct1 <- head(names(pgx1$gx.meta$meta), 2)
      ct2 <- head(names(pgx2$gx.meta$meta), 2)
      ct1 <- input.contrast1()
      ct2 <- input.contrast2()
      shiny::req(ct1)
      shiny::req(ct2)
      if (!all(ct1 %in% names(pgx1$gx.meta$meta))) {
        return(NULL)
      }
      if (!all(ct2 %in% names(pgx2$gx.meta$meta))) {
        return(NULL)
      }

      F1 <- playbase::pgx.getMetaMatrix(pgx1)$fc[, ct1, drop = FALSE]
      F2 <- playbase::pgx.getMetaMatrix(pgx2)$fc[, ct2, drop = FALSE]

      gg <- intersect(rownames(pgx1$X), rownames(pgx2$X))
      F1 <- F1[match(gg, rownames(F1)), , drop = FALSE]
      F2 <- F2[match(gg, rownames(F2)), , drop = FALSE]
      rownames(F1) <- gg
      rownames(F2) <- gg
      colnames(F1) <- paste0("1:", colnames(F1))
      colnames(F2) <- paste0("2:", colnames(F2))
      rho <- 1

      kk <- intersect(colnames(pgx1$X), colnames(pgx2$X))
      if (length(kk) >= 10) {
        X1 <- scale(t(pgx1$X[gg, kk]))
        X2 <- scale(t(pgx2$X[gg, kk]))
        rho <- colSums(X1 * X2) / (nrow(X1) - 1)
      }

      fc1 <- sqrt(rowMeans(F1**2))
      fc2 <- sqrt(rowMeans(F2**2))
      score <- rho * fc1 * fc2

      title <- pgx1$genes[gg, "gene_title"]
      title <- substring(title, 1, 60)

      df <- data.frame(title, score, rho, F1, F2, check.names = FALSE)
      df <- df[order(-df$score), ]
      df
    })

    hilightgenes <- shiny::reactive({
      genes <- as.character(input$genelist)
      genes <- strsplit(genes, split = "[\t, \n]")[[1]]
      gsub("[ ]", "", genes)
    })

    input.contrast1 <- shiny::reactive({
      input$contrast1
    }) %>% shiny::debounce(2500)

    input.contrast2 <- shiny::reactive({
      input$contrast2
    }) %>% shiny::debounce(2500)

    ## ============================================================================
    ## ScatterPlot 1
    ## ============================================================================


    createPlot <- function(pgx, pgx1, pgx2, ct, type, cex.lab, higenes, ntop) {
      p <- NULL
      genes1 <- rownames(pgx$X)
      higenes1 <- genes1[match(toupper(higenes), toupper(genes1))]

      # if(type %in% c('UMAP1','UMAP2')) {
      if (type %in% c("UMAP1", "UMAP2")) {
        if (type == "UMAP1") {
          pos <- pgx1$cluster.genes$pos[["umap2d"]]
        } else if (type == "UMAP2") {
          pos <- pgx2$cluster.genes$pos[["umap2d"]]
        }
        dim(pos)
        gg <- intersect(toupper(rownames(pgx$X)), toupper(rownames(pos)))
        jj <- match(gg, toupper(rownames(pos)))
        pos <- pos[jj, ]
        ii <- match(toupper(rownames(pos)), toupper(rownames(pgx$X)))
        rownames(pos) <- rownames(pgx$X)[ii]

        p <- playbase::pgx.plotGeneUMAP(
          pgx,
          contrast = ct, pos = pos,
          cex = 0.9, cex.lab = cex.lab,
          hilight = higenes1, ntop = ntop,
          zfix = TRUE, par.sq = TRUE,
          plotlib = "base"
        )
      } else if (type == "heatmap") {
        gg <- intersect(toupper(higenes), toupper(rownames(pgx$X)))
        if (length(gg) > 1) {
          jj <- match(gg, toupper(rownames(pgx$X)))
          X1 <- pgx$X[jj, , drop = FALSE]
          Y1 <- pgx$samples
          playbase::gx.splitmap(X1,
            nmax = 40, col.annot = Y1,
            softmax = TRUE, show_legend = FALSE
          )
        }
      } else {
        genes1 <- rownames(pgx$X)
        gg <- intersect(toupper(higenes), toupper(genes1))
        higenes1 <- genes1[match(gg, toupper(genes1))]
        p <- playbase::pgx.plotContrast(
          pgx,
          contrast = ct, hilight = higenes1,
          ntop = ntop, cex.lab = cex.lab, ## dlim=0.06,
          par.sq = TRUE, type = type, plotlib = "base"
        )
      }
      p <- grDevices::recordPlot()
      p
    }


    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    # Dataset 1

    compare_plot_compare1_server(
      "dt1",
      pgx = pgx,
      input.contrast1 = input.contrast1,
      hilightgenes = hilightgenes,
      createPlot = createPlot,
      plottype = shiny::reactive(input$plottype),
      dataset2 = dataset2,
      watermark = WATERMARK
    )

    # Dataset 2

    compare_plot_compare2_server(
      "dt2",
      pgx = pgx,
      input.contrast2 = input.contrast2,
      hilightgenes = hilightgenes,
      createPlot = createPlot,
      plottype = shiny::reactive(input$plottype),
      dataset2 = dataset2,
      watermark = WATERMARK
    )

    # FC Correlation

    compare_plot_fc_correlation_server(
      "fcfcplot",
      pgx = pgx,
      dataset2 = dataset2,
      hilightgenes = hilightgenes,
      input.contrast1 = input.contrast1,
      input.contrast2 = input.contrast2,
      watermark = WATERMARK
    )

    # Cumulative FC

    compare_plot_cum_fc1_server(
      "cumfcplot1",
      pgx = pgx,
      dataset2 = dataset2,
      cum_fc = cum_fc,
      input.contrast1 = input.contrast1,
      input.contrast2 = input.contrast2,
      watermark = WATERMARK
    )

    compare_plot_cum_fc2_server(
      "cumfcplot2",
      pgx = pgx,
      dataset2 = dataset2,
      cum_fc = cum_fc,
      input.contrast1 = input.contrast1,
      input.contrast2 = input.contrast2,
      watermark = WATERMARK
    )

    # Correlation score

    score_table <- compare_table_corr_score_server(
      "score_table",
      getOmicsScoreTable = getOmicsScoreTable,
      watermark = WATERMARK
    )

    # Expression

    compare_plot_expression_server(
      "multibarplot",
      pgx = pgx,
      dataset2 = dataset2,
      input.contrast1 = input.contrast1,
      input.contrast2 = input.contrast2,
      hilightgenes = hilightgenes,
      getOmicsScoreTable = getOmicsScoreTable,
      score_table = score_table,
      watermark = WATERMARK
    )

    # Gene correlation

    compare_plot_gene_corr_server(
      "genecorr",
      pgx = pgx,
      dataset2 = dataset2,
      input.contrast1 = input.contrast1,
      input.contrast2 = input.contrast2,
      hilightgenes = hilightgenes,
      getOmicsScoreTable = getOmicsScoreTable,
      score_table = score_table,
      contrast1 = shiny::reactive(input$contrast1),
      watermark = WATERMARK
    )
  })
} ## end-of-Board
