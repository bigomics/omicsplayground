##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

if (0) {
  pgx <- playdata::GEIGER_PGX
  rawX <- pgx$X
  normX <- pgx$X
  X <- normX
  Y <- pgx$samples[, c(2, 4, 5)]
}

PCAexplorer_server <- function(id, pgx) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      pgx <- playdata::GEIGER_PGX
      rX <- shiny::reactive(pgx$X)
      rY <- shiny::reactive({
        Y <- pgx$samples
        ngrp <- apply(Y, 2, function(y) length(table(y[!is.na(y)])))
        sel.valid <- ngrp > 1 & ngrp < colSums(!is.na(Y))
        Y[, sel.valid, drop = FALSE]
      })

      get_pcaX <- shiny::reactive({
        rawX <- rX()
        Y <- rY()

        m <- "CPM+quantile"
        normX <- playbase::normalizeExpression(rawX, method = m,
          ref = NULL, prior = 0)
        
        cX <- normX - rowMeans(normX, na.rm = TRUE)  ## important
        sel <- which(rowMeans(is.na(cX)) == 0)  ## only complete rows
        nv <- min(10, dim(cX) - 1)
        res <- irlba::irlba(cX[sel, ], nv = nv)
        rownames(res$u) <- rownames(cX)[sel]
        rownames(res$v) <- colnames(cX)
        colnames(res$u) <- paste0("PC", 1:ncol(res$u))
        colnames(res$v) <- paste0("PC", 1:ncol(res$v))

        H <- playbase::expandPhenoMatrix(Y, drop.ref = FALSE)
        res$R <- stats::cor(H, res$v)

        res$X <- normX
        res$Y <- Y
        res$H <- H

        res
      })

      ## ------------------------------------------------------------------
      ## Sidebar visibility — CSS class only, and only when board changes.
      ##
      ## Do NOT use shinyjs::toggle()/show()/hide(): jQuery .show()/.hide()
      ## re-apply inline display on every call and re-init ionRangeSliders,
      ## which flickers the settings panel on every plot-driven flush.
      ## toggleClass with a cached board is a no-op when nothing changed.
      ## Nested arrow opts are handled client-side via conditionalPanel.
      ## ------------------------------------------------------------------
      board <- shiny::reactive({
        b <- input$board
        if (is.null(b) || !nzchar(as.character(b))) "Biplot" else as.character(b)
      })

      shiny::observeEvent(board(), {
        b <- board()
        need <- list(
          sb_colorby     = b %in% c("Biplot", "PCA pairs"),
          sb_arrowtype   = identical(b, "Biplot"),
          sb_show_arrows = b %in% c("Biplot", "Feature PCA"),
          sb_arrow_opts  = b %in% c("Biplot", "Feature PCA"),
          sb_min_fc      = identical(b, "Feature PCA"),
          sb_ellipse     = b %in% c("Biplot", "PCA pairs")
        )
        for (nm in names(need)) {
          ## condition = TRUE  → add hide class when widget is NOT needed
          shinyjs::toggleClass(
            id = nm,
            class = "pca-sb-hide",
            condition = !isTRUE(need[[nm]])
          )
        }
      }, ignoreNULL = TRUE)

      ## seed selectInput choices from data (widgets stay mounted)
      shiny::observeEvent(rY(), {
        cols <- colnames(rY())
        shiny::updateSelectInput(session, "colorby", choices = cols)
      }, ignoreNULL = TRUE)

      ## Stretch min_fc slider once when PCA is ready (not on every plot tick)
      shiny::observeEvent(get_pcaX(), {
        res <- get_pcaX()
        normX <- res$X
        shiny::req(res, normX)
        H <- if (!is.null(res$H)) res$H else {
          playbase::expandPhenoMatrix(rY(), drop.ref = FALSE)
        }
        genes <- intersect(rownames(res$u), rownames(normX))
        if (length(genes) < 2L) return()
        max_fc <- pcaexplorer.feature_max_fc(normX[genes, , drop = FALSE], H)
        max_fc <- max_fc[is.finite(max_fc)]
        if (!length(max_fc)) return()
        hi <- as.numeric(stats::quantile(max_fc, 0.99, na.rm = TRUE))
        if (!is.finite(hi) || hi <= 0) hi <- max(max_fc, na.rm = TRUE)
        if (!is.finite(hi) || hi <= 0) hi <- 1
        hi <- round(hi, 2)
        step <- max(0.01, round(hi / 40, 2))
        shiny::updateSliderInput(
          session, "min_fc",
          min = 0, max = hi,
          step = step,
          value = 0
        )
      }, ignoreNULL = TRUE, once = TRUE)

      ## shared sidebar reactives for board modules
      colorby <- shiny::reactive(input$colorby)
      numarrows <- shiny::reactive({
        n <- suppressWarnings(as.integer(input$numarrows))
        if (!is.finite(n) || n < 1L) 10L else n
      })
      show_arrows <- shiny::reactive(isTRUE(input$show_arrows))
      arrowtype <- shiny::reactive(input$arrowtype)
      ellipse <- shiny::reactive(isTRUE(input$ellipse))
      min_fc <- shiny::reactive({
        v <- suppressWarnings(as.numeric(input$min_fc))
        if (!is.finite(v) || v < 0) 0 else v
      })
      var_scale <- shiny::reactive({
        v <- suppressWarnings(as.numeric(input$var_scale))
        if (!is.finite(v) || v <= 0) 0.65 else v
      })

      PCAexplorer_biplot_module(
        "biplot",
        get_pcaX,
        colorby = colorby,
        numarrows = numarrows,
        show_arrows = show_arrows,
        arrowtype = arrowtype,
        ellipse = ellipse,
        var_scale = var_scale
      )

      PCAexplorer_pairs_module(
        "pairs",
        get_pcaX,
        colorby = colorby,
        ellipse = ellipse
      )

      PCAexplorer_loadings_module(
        "loadings",
        get_pcaX
      )

      PCAexplorer_feature_module(
        "feature",
        get_pcaX,
        rY,
        numarrows = numarrows,
        show_arrows = show_arrows,
        min_fc = min_fc,
        var_scale = var_scale
      )


    } ## end-of-server
  )
}
