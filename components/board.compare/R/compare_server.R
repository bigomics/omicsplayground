##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

CompareBoard <- function(id, pgx, pgx_dir = reactive(file.path(OPG, "data", "mini-example")), labeltype = shiny::reactive("feature")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 770 # row height of panel
    tabH <- "70vh"

    infotext <-
      tspan("The <strong>Compare Datasets</strong> module enables users to compare their dataset to other datasets. This module allows side-by-side comparison of volcano, scatter or gene t-SNE plots. It provides pairwise correlation plots and/or enrichment plots with signatures from other data sets. <br><br><br><br>
        <center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=5' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>",
        js = FALSE
      )

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

    score_table_rows <- reactive({
      score_table$rownames_all()
    })

    ## upon new pgx upload
    shiny::observeEvent(
      {
        list(pgx$X)
      },
      {
        comparisons1 <- names(pgx$gx.meta$meta)
        sel1 <- comparisons1[1]
        shiny::updateSelectInput(session, "contrast1",
          choices = comparisons1,
          selected = sel1
        )
        contrast1(sel1)

        pgx.files <- sort(dir(pgx_dir(), pattern = "pgx$"))
        shiny::updateSelectInput(session, "dataset2",
          choices = c("<this>", pgx.files),
          sel = "<this>"
        )

        sel2 <- rev(head(comparisons1, 2))[1]
        shiny::updateSelectInput(session, "contrast2",
          choices = comparisons1,
          selected = sel2
        )
        contrast2(sel2)

        shiny::updateTextAreaInput(
          session, "genelist",
          placeholder = tspan("Paste your custom gene list", js = FALSE)
        )
      }
    )

    ## keep a list of highlighted/selected features
    hilightgenes <- reactiveVal(NULL)
    shiny::observe({
      shiny::req(contrast1())
      shiny::req(contrast2())

      df <- getScoreTable()
      req(df)
      ntop <- as.integer(input$ntop)
      higenes <- NULL
      if (input$hilighttype == "top scoring") {
        higenes <- rownames(df)[order(df$score, decreasing = TRUE)]
        higenes <- head(higenes, ntop)
      }
      if (input$hilighttype == "custom") {
        higenes <- input$genelist
        higenes <- strsplit(higenes, split = "[\t, \n]")[[1]]
        higenes <- trimws(gsub("[ ]", "", higenes))
      }
      hilightgenes(higenes)
    })


    ## allow trigger on explicit compare button
    contrast1 <- shiny::reactiveVal()
    contrast2 <- shiny::reactiveVal()
    shiny::observeEvent(
      {
        list(pgx$X, input$compare_button)
      },
      {
        shiny::req(pgx$X)
        shiny::req(dataset2()$X)

        ## check if contrast1 is NULL or invalid (e.g. from previous
        ## used dataset), if invalid, force update
        ct1 <- input$contrast1
        all.ct1 <- playbase::pgx.getContrasts(pgx)
        valid.ct1 <- !is.null(ct1) && all(ct1 %in% all.ct1)
        if (!valid.ct1) {
          ct1 <- all.ct1[1]
          shiny::updateSelectInput(session, "contrast1",
            choices = all.ct1,
            selected = ct1
          )
        }
        contrast1(ct1)

        ## check if contrast2 is NULL or invalid (e.g. from previous
        ## used dataset), if invalid, force update
        ct2 <- input$contrast2
        pgx2 <- dataset2()
        all.ct2 <- playbase::pgx.getContrasts(pgx2)
        valid.ct2 <- !is.null(ct2) && all(ct2 %in% all.ct2)
        if (!valid.ct2) {
          ct2 <- tail(head(all.ct2, 2), 1)
          shiny::updateSelectInput(session, "contrast2",
            choices = all.ct2,
            selected = ct2
          )
        }
        contrast2(ct2)
      },
      ignoreInit = FALSE,
      ignoreNULL = FALSE
    )


    # Retrieve the 2nd dataset
    dataset2 <- shiny::eventReactive(
      {
        list(pgx$X, input$dataset2)
      },
      {
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
        shiny::updateSelectInput(
          session, "contrast2",
          choices = comparisons2,
          selected = sel2
        )
        contrast2(sel2)

        return(pgx)
      }
    )

    ## ============================================================================
    ## ========================= REACTIVE FUNCTIONS ===============================
    ## ============================================================================

    # Cummulative FC
    getMatrices <- shiny::reactive({
      pgx1 <- pgx
      pgx2 <- dataset2()
      shiny::req(pgx1$X)
      shiny::req(pgx2$X)

      ct1 <- contrast1()
      ct2 <- contrast2()
      valid.ct1 <- !is.null(ct1) && all(ct1 %in% playbase::pgx.getContrasts(pgx1))
      valid.ct2 <- !is.null(ct2) && all(ct2 %in% playbase::pgx.getContrasts(pgx2))
      if (!valid.ct1 || !valid.ct2) {
        ## if NULL refresh
        dbg("[getMatrices] WARNING contrast RV not valid! force update! ")
        contrast1(input$contrast1)
        contrast2(input$contrast2)
        shinyjs::click("compare_button")
        shinyjs::click(ns("compare_button"))
        return(NULL)
      }

      org1 <- tolower(playbase::pgx.getOrganism(pgx1))
      org2 <- tolower(playbase::pgx.getOrganism(pgx2))

      F1 <- playbase::pgx.getMetaMatrix(pgx1)$fc[, ct1, drop = FALSE]
      F2 <- playbase::pgx.getMetaMatrix(pgx2)$fc[, ct2, drop = FALSE]

      target_col <- "rownames"
      if (is.null(pgx1$version) && is.null(pgx2$version)) {
        # For old versions
        target_col <- "rownames"
        rownames(F1) <- toupper(rownames(F1))
        rownames(F2) <- toupper(rownames(F2))
      } else if (pgx1$name == pgx2$name && nrow(pgx1$X) == nrow(pgx2$X)) {
        # For same dataset. we compare on rownames
        target_col <- "rownames"
      } else if (org1 == org2) {
        # For same org. we ensure compare on symbol
        target_col <- "symbol"
      } else if (org1 != org2) {
        # For different org. we ensure compare on human_ortholog
        # If it is not present, use gene_name
        target_col <- "human_ortholog"
      }
      target_col

      if (target_col != "rownames") {
        F1 <- playbase::rename_by(F1, pgx1$genes, target_col)
        F2 <- playbase::rename_by(F2, pgx2$genes, target_col)
      }

      F1 <- playbase::rowmean(F1, rownames(F1))
      F2 <- playbase::rowmean(F2, rownames(F2))
      gg <- intersect(rownames(F1), rownames(F2))
      F1 <- F1[gg, , drop = FALSE]
      F2 <- F2[gg, , drop = FALSE]
      colnames(F1) <- paste0("1:", colnames(F1))
      colnames(F2) <- paste0("2:", colnames(F2))

      ## ---------- X matrices ---------------
      X1 <- pgx1$X
      X2 <- pgx2$X
      if (target_col != "rownames") {
        X1 <- playbase::rename_by(X1, pgx1$genes, target_col)
        X2 <- playbase::rename_by(X2, pgx2$genes, target_col)
      }

      ## collapse duplicates
      X1 <- playbase::rowmean(X1, rownames(X1))
      X2 <- playbase::rowmean(X2, rownames(X2))
      X1 <- X1[match(rownames(F1), rownames(X1)), ]
      X2 <- X2[match(rownames(F2), rownames(X2)), ]

      ## compute correlation if we have enough common samples
      rho <- NULL
      kk <- intersect(colnames(pgx1$X), colnames(pgx2$X))
      if (length(kk) >= 10) {
        scaled_X1 <- scale(Matrix::t(X1[, kk]))
        scaled_X2 <- scale(Matrix::t(X2[, kk]))
        rho <- colSums(scaled_X1 * scaled_X2, na.rm = TRUE) / (nrow(X1) - 1)
      }

      ## get gene_title
      if (target_col == "rownames") {
        match.target <- rownames(pgx1$genes)
      } else {
        match.target <- pgx1$genes[, target_col]
      }
      gene_title <- pgx1$genes[match(rownames(F1), match.target), "gene_title"]

      pos1 <- pos2 <- NULL
      if (0) {
        ## we also need collapsed/aligned UMAP positions
        pos1 <- pgx1$cluster.genes$pos[["umap2d"]]
        pos2 <- pgx2$cluster.genes$pos[["umap2d"]]
        if (target_col != "rownames") {
          pos1 <- playbase::rename_by(pos1, pgx1$genes, target_col)
          pos2 <- playbase::rename_by(pos2, pgx2$genes, target_col)
        }
        pos1 <- playbase::rowmean(pos1, rownames(pos1))
        pos2 <- playbase::rowmean(pos2, rownames(pos2))
        pos1 <- pos1[match(rownames(F1), rownames(pos1)), ]
        pos2 <- pos1[match(rownames(F2), rownames(pos2)), ]
      }

      list(
        F1 = F1,
        F2 = F2,
        X1 = X1,
        X2 = X2,
        pos1 = pos1,
        pos2 = pos2,
        rho = rho,
        gene_title = gene_title,
        target_col = target_col
      )
    })

    getScoreTable <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(dataset2())
      shiny::req(contrast1()) ## trigger on button
      shiny::req(contrast2())

      valid.ct1 <- all(contrast1() %in% playbase::pgx.getContrasts(pgx))
      valid.ct2 <- all(contrast2() %in% playbase::pgx.getContrasts(dataset2()))
      shiny::req(valid.ct1 && valid.ct2)

      pgx1 <- pgx
      pgx2 <- dataset2()

      res <- getMatrices()
      F1 <- res$F1
      F2 <- res$F2
      rho <- res$rho

      ## compute score
      fc1 <- sqrt(rowMeans(F1**2, na.rm = TRUE))
      fc2 <- sqrt(rowMeans(F2**2, na.rm = TRUE))
      if (!is.null(rho)) {
        score <- rho * fc1 * fc2
      } else {
        score <- fc1 * fc2
      }
      if (is.null(rho)) rho <- rep(NA, nrow(F1))

      # get gene_title
      title <- res$gene_title

      df <- data.frame(title, score, rho, F1, F2, check.names = FALSE)
      rownames(df) <- rownames(F1)

      return(df)
    })

    ## ============================================================================
    ## ScatterPlot 1
    ## ============================================================================

    createPlot <- function(pgx, pgx1, pgx2, ct, target_col, type, cex.lab,
                           higenes, ntop, get_data = FALSE, labeltype = shiny::reactive("feature")) {
      p <- NULL
      ## map hilighted genes to pgx probes
      label <- playbase::map_probes(pgx$genes, higenes, ignore.case = TRUE)
      if (type %in% c("UMAP1", "UMAP2")) {
        if (type == "UMAP1") {
          pos <- pgx1$cluster.genes$pos[["umap2d"]]
        } else if (type == "UMAP2") {
          pos <- pgx2$cluster.genes$pos[["umap2d"]]
        }
        posx <- playbase::rename_by(pos, pgx1$genes, new_id = target_col)
        mapped.pos <- playbase::rename_by2(posx, pgx$genes,
          new_id = "rownames",
          unique = FALSE
        )
        mapped.pos <- mapped.pos[match(rownames(pgx$X), rownames(mapped.pos)), ]
        rownames(mapped.pos) <- rownames(pgx$X)

        p <- playbase::pgx.plotGeneUMAP(
          pgx,
          contrast = ct,
          pos = mapped.pos,
          cex = 0.9,
          cex.lab = cex.lab,
          hilight = NULL,
          label = label,
          ntop = ntop,
          zfix = TRUE,
          par.sq = TRUE,
          plotlib = "base",
          data = get_data,
          labeltype = labeltype()
        )
        ## } else if (type == "heatmap") {
        ##   gg <- intersect(toupper(higenes), toupper(rownames(pgx$X)))
        ##     if (length(gg) > 1) {
        ##       jj <- match(gg, toupper(rownames(pgx$X)))
        ##       X1 <- pgx$X[jj, , drop = FALSE]
        ##       Y1 <- pgx$samples
        ##       if (get_data) {
        ##         return(
        ##           playbase::gx.splitmap(X1,
        ##             nmax = 40, col.annot = Y1,
        ##             softmax = TRUE, show_legend = FALSE,
        ##             data = TRUE
        ##           )
        ##         )
        ##       }
        ##       playbase::gx.splitmap(X1,
        ##         nmax = 40, col.annot = Y1,
        ##         softmax = TRUE, show_legend = FALSE
        ##       )
        ##     }
      } else {
        p <- playbase::pgx.plotContrast(
          pgx,
          contrast = ct,
          hilight = hilight,
          ntop = ntop,
          cex.lab = cex.lab, #
          par.sq = TRUE,
          type = type,
          plotlib = "base",
          data = get_data
        )
      }
      if (get_data) {
        return(p)
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
      contrast1 = contrast1,
      hilightgenes = hilightgenes,
      createPlot = createPlot,
      plottype = shiny::reactive(input$plottype),
      dataset2 = dataset2,
      getMatrices = getMatrices,
      watermark = WATERMARK,
      labeltype = labeltype
    )

    # Dataset 2
    compare_plot_compare2_server(
      "dataset2",
      pgx = pgx,
      contrast2 = contrast2,
      ## contrast1 = contrast1,
      hilightgenes = hilightgenes,
      createPlot = createPlot,
      plottype = shiny::reactive(input$plottype),
      dataset2 = dataset2,
      getMatrices = getMatrices,
      watermark = WATERMARK,
      labeltype = labeltype
    )

    # FC Correlation
    compare_plot_fcfc_server(
      "fcfcplot",
      getMatrices = getMatrices,
      hilightgenes = hilightgenes,
      watermark = WATERMARK,
      labeltype = labeltype,
      pgx = pgx
    )

    # Cumulative FC
    compare_plot_cum_fc1_server(
      "cumfcplot1",
      pgx = pgx,
      labeltype = labeltype,
      getMatrices = getMatrices,
      watermark = WATERMARK
    )

    compare_plot_cum_fc2_server(
      "cumfcplot2",
      pgx = pgx,
      labeltype = labeltype,
      getMatrices = getMatrices,
      watermark = WATERMARK
    )

    # Correlation score

    score_table <- compare_table_corr_score_server(
      "score_table",
      getScoreTable = getScoreTable,
      watermark = WATERMARK
    )

    # Expression

    compare_plot_expression_server(
      "multibarplot",
      pgx = pgx,
      dataset2 = dataset2,
      contrast1 = contrast1,
      contrast2 = contrast2,
      hilightgenes = hilightgenes,
      getMatrices = getMatrices,
      getScoreTable = getScoreTable,
      selected = score_table_rows,
      watermark = WATERMARK
    )

    # Gene correlation

    compare_plot_genecorr_server(
      "genecorr",
      pgx = pgx,
      dataset2 = dataset2,
      contrast1 = contrast1,
      contrast2 = contrast2,
      getMatrices = getMatrices,
      getScoreTable = getScoreTable,
      selected = score_table_rows,
      watermark = WATERMARK
    )
  })
} ## end-of-Board
