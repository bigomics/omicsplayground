##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

CompareBoard <- function(id, pgx, pgx_dir = reactive(file.path(OPG, "data", "mini-example"))) {
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

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Compare expression" = list(
        enable = NULL,
        disable = c("pcor_ntop")
      ),
      "Foldchange" = list(
        enable = NULL,
        disable = c("plottype")
      ),
      "Gene Correlation" = list(
        enable = NULL,
        disable = c("plottype", "ntop", "hilighttype", "genelist")
      )
    )

    shiny::observeEvent(input$tabs1, {
      bigdash::update_tab_elements(input$tabs1, tab_elements)
    })

    score_table <- reactiveVal(NULL)

    shiny::observe({
      shiny::req(score_table_temp$rows_all())
      score_table(score_table_temp$rows_all())
    })

    shiny::observe({
      comparisons1 <- names(pgx$gx.meta$meta)
      sel1 <- comparisons1[1]
      shiny::updateSelectInput(session, "contrast1", choices = comparisons1, selected = sel1)
      pgx.files <- sort(dir(pgx_dir(), pattern = "pgx$"))
      shiny::updateSelectInput(session, "dataset2", choices = c("<this>", pgx.files))
    })
    hilightgenes <- reactiveVal(NULL)

    shiny::observe({
      shiny::req(contrast1())
      shiny::req(contrast2())

      df <- getOmicsScoreTable()
      req(df)
      ntop <- as.integer(input$ntop)
      higenes <- rownames(df)[order(df$score, decreasing = TRUE)]
      higenes <- head(higenes, ntop)
      higenes <- paste(higenes, collapse = " ")
      if (input$hilighttype == "custom") {
        higenes <- input$genelist
      }
      hilightgenes({
        genes <- strsplit(higenes, split = "[\t, \n]")[[1]]
        gsub("[ ]", "", genes)
      })
      shiny::updateTextAreaInput(session, "genelist", value = higenes)
    })

    ## allow trigger on explicit compare button
    contrast1 <- shiny::reactiveVal()
    contrast2 <- shiny::reactiveVal()
    shiny::observeEvent(input$compare_button, {
      contrast1(input$contrast1)
      contrast2(input$contrast2)
    })

    ## ================================================================================
    ## ========================= REACTIVE FUNCTIONS ===================================
    ## ================================================================================

    # Retreive the 2nd dataset
    dataset2 <- shiny::reactive({
      shiny::req(input$dataset2)
      if (input$dataset2 == "<this>") {
        pgx <- pgx
      } else {
        file2 <- file.path(pgx_dir(), input$dataset2)
        pgx <- playbase::pgx.load(file2)
        pgx <- playbase::pgx.initialize(pgx)
      }
      comparisons2 <- names(pgx$gx.meta$meta)
      sel2 <- tail(head(comparisons2, 2), 1)
      shiny::updateSelectInput(session, "contrast2", choices = comparisons2, selected = sel2)
      return(pgx)
    })

    # Cummulative FC
    cum_fc <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(dataset2())
      shiny::req(contrast1())
      shiny::req(contrast2())
      pgx1 <- pgx
      pgx2 <- dataset2()
      org1 <- playbase::pgx.getOrganism(pgx1)
      org2 <- playbase::pgx.getOrganism(pgx2)

      ct1 <- head(names(pgx1$gx.meta$meta), 2)
      ct2 <- head(names(pgx2$gx.meta$meta), 2)
      ct1 <- contrast1()
      ct2 <- contrast2()

      shiny::validate(shiny::need(all(ct1 %in% names(pgx1$gx.meta$meta)), "Warning: No common contrasts."))
      shiny::validate(shiny::need(all(ct2 %in% names(pgx2$gx.meta$meta)), "Warning: No common contrasts."))

      F1 <- playbase::pgx.getMetaMatrix(pgx1)$fc[, ct1, drop = FALSE]
      F2 <- playbase::pgx.getMetaMatrix(pgx2)$fc[, ct2, drop = FALSE]
      gg <- intersect(toupper(rownames(F1)), toupper(rownames(F2)))

      if (is.null(pgx1$version) && is.null(pgx2$version)) {
        gg <- intersect(toupper(rownames(F1)), toupper(rownames(F2)))
        g1 <- rownames(F1)[match(gg, toupper(rownames(F1)))]
        g2 <- rownames(F2)[match(gg, toupper(rownames(F2)))]
      } else if (org1 == org2) {
        # For same org. we ensure compare on symbol
        F1 <- playbase::rename_by(F1, pgx1$genes, "symbol")
        F2 <- playbase::rename_by(F2, pgx2$genes, "symbol")
        gg <- intersect(rownames(F1), rownames(F2))
        g1 <- g2 <- gg
      } else if (org1 != org2) {
        # For different org. we ensure compare on human_ortholog
        # If it is not present, use gene_name
        target_col1 <- target_col2 <- "human_ortholog"
        if (!target_col1 %in% colnames(pgx1$genes)) target_col1 <- "gene_name"
        if (!target_col2 %in% colnames(pgx2$genes)) target_col2 <- "gene_name"
        F1 <- playbase::rename_by(F1, pgx1$genes, target_col1)
        F2 <- playbase::rename_by(F2, pgx2$genes, target_col2)
        gg <- intersect(rownames(F1), rownames(F2))
        g1 <- g2 <- gg
      }

      F1 <- F1[g1, , drop = FALSE]
      F2 <- F2[g2, , drop = FALSE]

      # TODO: implement average or sum
      F1 <- F1[!duplicated(rownames(F1)), , drop = FALSE]
      F2 <- F2[!duplicated(rownames(F2)), , drop = FALSE]

      colnames(F1) <- paste0("1:", colnames(F1))
      colnames(F2) <- paste0("2:", colnames(F2))

      return(cbind(F1, F2))
    })

    getOmicsScoreTable <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(dataset2())
      shiny::req(contrast1()) ## trigger on button
      shiny::req(contrast2())
      shiny::req(
        contrast2() %in% colnames(playbase::pgx.getMetaMatrix(dataset2())$fc)
      )
      pgx1 <- pgx
      pgx2 <- dataset2()
      org1 <- playbase::pgx.getOrganism(pgx1)
      org2 <- playbase::pgx.getOrganism(pgx2)

      F1_2 <- cum_fc()
      F1 <- F1_2[, 1, drop = FALSE]
      F2 <- F1_2[, 2, drop = FALSE]
      rho <- 1
      gg <- rownames(F1_2)

      if (is.null(pgx1$version) && is.null(pgx2$version)) {
        target_col1 <- target_col2 <- "gene_name"
      } else if (org1 == org2) {
        target_col1 <- target_col2 <- "symbol"
      } else if (org1 != org2) {
        target_col1 <- target_col2 <- "human_ortholog"
        if (!target_col1 %in% colnames(pgx1$genes)) target_col1 <- "gene_name"
        if (!target_col2 %in% colnames(pgx2$genes)) target_col2 <- "gene_name"
      }
      kk <- intersect(colnames(pgx1$X), colnames(pgx2$X))
      if (length(kk) >= 10) {
        X1 <- playbase::rename_by(pgx1$X, pgx1$genes, target_col1)
        X2 <- playbase::rename_by(pgx2$X, pgx2$genes, target_col2)
        X1 <- scale(t(X1[gg, kk]))
        X2 <- scale(t(X2[gg, kk]))
        rho <- colSums(X1 * X2) / (nrow(X1) - 1)
      }
      fc1 <- sqrt(rowMeans(F1**2))
      fc2 <- sqrt(rowMeans(F2**2))
      score <- rho * fc1 * fc2

      # TODO: get teh gene_title
      # mii <- which(pgx1$genes[ , target_col1] %in% gg)
      # ii <- ii[!duplicated(ii)]
      title <- pgx1$genes[gg, "gene_title"]
      title <- substring(title, 1, 60)

      df <- data.frame(title, score, rho, F1, F2, check.names = FALSE)
      df <- df[order(-df$score), ]

      return(df)
    })

    ## ============================================================================
    ## ScatterPlot 1
    ## ============================================================================


    createPlot <- function(pgx, pgx1, pgx2, ct, type, cex.lab, higenes, ntop, get_data = FALSE) {
      p <- NULL
      genes1 <- rownames(pgx$X)
      higenes1 <- genes1[match(toupper(higenes), toupper(genes1))]
      if (type %in% c("UMAP1", "UMAP2")) {
        if (type == "UMAP1") {
          pos <- pgx1$cluster.genes$pos[["umap2d"]]
        } else if (type == "UMAP2") {
          pos <- pgx2$cluster.genes$pos[["umap2d"]]
        }

        gg <- intersect(toupper(rownames(pgx$X)), toupper(rownames(pos)))
        jj <- match(gg, toupper(rownames(pos)))
        pos <- pos[jj, ]
        ii <- match(toupper(rownames(pos)), toupper(rownames(pgx$X)))
        rownames(pos) <- rownames(pgx$X)[ii]
        if (get_data) {
          return(
            playbase::pgx.plotGeneUMAP(
              pgx,
              contrast = ct, pos = pos,
              cex = 0.9, cex.lab = cex.lab,
              hilight = higenes1, ntop = ntop,
              zfix = TRUE, par.sq = TRUE,
              plotlib = "base", data = TRUE
            )
          )
        }
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
          if (get_data) {
            return(
              playbase::gx.splitmap(X1,
                nmax = 40, col.annot = Y1,
                softmax = TRUE, show_legend = FALSE,
                data = TRUE
              )
            )
          }
          playbase::gx.splitmap(X1,
            nmax = 40, col.annot = Y1,
            softmax = TRUE, show_legend = FALSE
          )
        }
      } else {
        genes1 <- rownames(pgx$X)
        gg <- intersect(toupper(higenes), toupper(genes1))
        higenes1 <- genes1[match(gg, toupper(genes1))]
        if (get_data) {
          return(
            playbase::pgx.plotContrast(
              pgx,
              contrast = ct, hilight = higenes1,
              ntop = ntop, cex.lab = cex.lab, #
              par.sq = TRUE, type = type, plotlib = "base",
              data = TRUE
            )
          )
        }
        p <- playbase::pgx.plotContrast(
          pgx,
          contrast = ct, hilight = higenes1,
          ntop = ntop, cex.lab = cex.lab, #
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
      "dataset1",
      pgx = pgx,
      ## input.contrast1 = shiny::reactive(input$contrast1),
      input.contrast1 = contrast1,
      ## input.contrast2 = shiny::reactive(input$contrast2),
      hilightgenes = hilightgenes,
      createPlot = createPlot,
      plottype = shiny::reactive(input$plottype),
      dataset2 = dataset2,
      ## compute = shiny::reactive(input$compare_button),
      watermark = WATERMARK
    )

    # Dataset 2
    compare_plot_compare2_server(
      "dataset2",
      pgx = pgx,
      input.contrast2 = shiny::reactive(input$contrast2),
      input.contrast1 = shiny::reactive(input$contrast1),
      hilightgenes = hilightgenes,
      createPlot = createPlot,
      plottype = shiny::reactive(input$plottype),
      dataset2 = dataset2,
      compute = shiny::reactive(input$compare_button),
      watermark = WATERMARK
    )

    # FC Correlation

    compare_plot_fc_correlation_server(
      "fcfcplot",
      cum_fc = cum_fc,
      hilightgenes = hilightgenes,
      input.contrast1 = reactive(input$contrast1),
      input.contrast2 = reactive(input$contrast2),
      compute = shiny::reactive(input$compare_button),
      watermark = WATERMARK
    )

    # Cumulative FC

    compare_plot_cum_fc1_server(
      "cumfcplot1",
      pgx = pgx,
      dataset2 = dataset2,
      cum_fc = cum_fc,
      ## compute = shiny::reactive(input$compare_button),
      watermark = WATERMARK
    )

    compare_plot_cum_fc2_server(
      "cumfcplot2",
      pgx = pgx,
      dataset2 = dataset2,
      cum_fc = cum_fc,
      compute = shiny::reactive(input$compare_button),
      watermark = WATERMARK
    )

    # Correlation score

    score_table_temp <- compare_table_corr_score_server(
      "score_table",
      getOmicsScoreTable = getOmicsScoreTable,
      watermark = WATERMARK
    )

    # Expression

    compare_plot_expression_server(
      "multibarplot",
      pgx = pgx,
      dataset2 = dataset2,
      input.contrast1 = shiny::reactive(input$contrast1),
      input.contrast2 = shiny::reactive(input$contrast2),
      hilightgenes = hilightgenes,
      getOmicsScoreTable = getOmicsScoreTable,
      score_table = score_table,
      compute = shiny::reactive(input$compare_button),
      watermark = WATERMARK
    )

    # Gene correlation

    compare_plot_gene_corr_server(
      "genecorr",
      pgx = pgx,
      dataset2 = dataset2,
      input.contrast1 = shiny::reactive(input$contrast1),
      input.contrast2 = shiny::reactive(input$contrast2),
      hilightgenes = hilightgenes,
      getOmicsScoreTable = getOmicsScoreTable,
      score_table = score_table,
      compute = shiny::reactive(input$compare_button),
      watermark = WATERMARK
    )
  })
} ## end-of-Board
