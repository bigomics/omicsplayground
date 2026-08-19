##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Server for the "PCA explorer" Qsee board (formerly the standalone
## app_pcaexplorer app, now merged in as an extra board). Thin wrapper:
## computation lives in app_qsee/R/pcaexplorer-compute.R, plotting in
## app_qsee/R/pcaexplorer-plots.R.
##

qsee_pcaexplorer_server <- function(id, rX, rY) {
  shiny::moduleServer(id, function(input, output, session) {

    OmicsBoard("board", pgx = NULL, title = "PCA explorer", infotext = NULL)
    is_visible <- bigdash::bd_is_visible(
      input, purge = qsee_purge_enabled(), label = "qsee_pcaexplorer_server"
    )
    redraw_tick <- bigdash::bd_redraw_tick(session = session)

    output$ui_output <- shiny::renderUI({
      qsee_pcaexplorer_ui_output(session$ns)
    })

    ## Qsee's own rY() does not drop degenerate phenotype columns
    ## (constant, or all-unique like a sample-id column); without that,
    ## qsee_pcaexplorer_plot_ftest_vs_pcdim() errors on an all-unique
    ## column ("length of 'dimnames' [2] not equal to array extent").
    rYf <- shiny::reactive({
      Y <- rY()
      shiny::req(Y)
      qsee_pcaexplorer_filter_pheno(Y)
    })

    shiny::observeEvent(rYf(), {
      cols <- colnames(rYf())
      shiny::updateSelectInput(session, "colorby", choices = cols, selected = cols[1])
    }, ignoreNULL = TRUE)

    ## Lazy: the plot/table outputs that read this are suspended while the
    ## tab is hidden. The one non-output consumer is the min_fc observer
    ## below, which is guarded on is_visible(). See qsee_visibility.R.
    get_pcaX <- shiny::reactive({
      rawX <- rX()
      Y <- rYf()
      shiny::req(rawX, Y)
      message("[qsee_pcaexplorer_server] computing...")
      qsee_pcaexplorer_compute(rawX, Y)
    })

    ## ------------------------------------------------------------------
    ## Sidebar visibility -- CSS class only, toggled on internal-tab change.
    ##
    ## Do NOT use shinyjs::toggle()/show()/hide(): jQuery .show()/.hide()
    ## re-apply inline display on every call and re-init ionRangeSliders,
    ## which flickers the settings panel. toggleClass is a no-op when
    ## nothing changed. shinyjs functions auto-namespace to this module's
    ## session, so plain (unprefixed) element ids below is correct.
    ## ------------------------------------------------------------------
    board_tab <- shiny::reactive({
      b <- input$board_tab
      if (is.null(b) || !nzchar(as.character(b))) "Biplot" else as.character(b)
    })

    shiny::observeEvent(board_tab(), {
      b <- board_tab()
      need <- list(
        sb_colorby     = b %in% c("Biplot", "PCA pairs"),
        sb_arrowtype   = identical(b, "Biplot"),
        sb_show_arrows = b %in% c("Biplot", "Feature PCA"),
        sb_arrow_opts  = b %in% c("Biplot", "Feature PCA"),
        sb_min_fc      = identical(b, "Feature PCA"),
        sb_ellipse     = b %in% c("Biplot", "PCA pairs")
      )
      for (nm in names(need)) {
        shinyjs::toggleClass(
          id = nm,
          class = "pca-sb-hide",
          condition = !isTRUE(need[[nm]])
        )
      }
    }, ignoreNULL = TRUE)

    ## Stretch min_fc slider once when PCA is first ready.
    ##
    ## This is the only eager consumer of get_pcaX(), so it must not pull
    ## before the board is visible -- without req(is_visible()) the PCA would
    ## be computed at startup even if this tab is never opened.
    ##
    ## Not observeEvent(..., once = TRUE): bindEvent() destroys the observer
    ## via on.exit(), so it self-destructs even when the handler req()s out,
    ## and the slider would never get stretched. Latch on a plain flag
    ## instead -- it is only read after an invalidation, so it needs no
    ## reactivity of its own.
    minfc_stretched <- FALSE
    shiny::observe({
      shiny::req(!minfc_stretched, is_visible())
      rng <- qsee_pcaexplorer_minfc_range(get_pcaX())
      shiny::req(rng)
      minfc_stretched <<- TRUE
      shiny::updateSliderInput(
        session, "min_fc",
        min = rng$min, max = rng$max, step = rng$step, value = rng$value
      )
    })

    ## shared sidebar reactives
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
    n_arrows <- shiny::reactive({
      if (isTRUE(show_arrows())) as.integer(numarrows()) else 0L
    })

    ## ------------------------------------------------------------------
    ## Biplot panel
    ## ------------------------------------------------------------------
    output$biplot_2d <- plotly::renderPlotly({
      redraw_tick()
      res <- get_pcaX()
      ph <- colorby()
      shiny::req(res, ph)
      qsee_pcaexplorer_plot_biplot_2d(
        res, res$Y, xpc = 1, ypc = 2,
        colorby_var = ph, arrows = arrowtype(),
        arrow.col = "black",
        numarrows = n_arrows(),
        var.scale = var_scale(),
        legend = TRUE,
        add.ellipse = ellipse()
      )
    })

    output$biplot_3d <- plotly::renderPlotly({
      redraw_tick()
      res <- get_pcaX()
      ph <- colorby()
      shiny::req(res, ph)
      qsee_pcaexplorer_plot_biplot_3d(
        res, res$Y, xpc = 1, ypc = 2, zpc = 3,
        colorby_var = ph, arrows = arrowtype(),
        arrow.col = "black",
        numarrows = n_arrows(),
        var.scale = var_scale(),
        legend = TRUE,
        add.ellipse = ellipse()
      )
    })

    output$biplot_table <- shiny::renderDataTable({
      res <- get_pcaX()
      shiny::req(res)
      dt <- round(res$v, 4)

      DT::datatable(
        data = dt,
        class = "compact hover",
        rownames = TRUE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(
          mode = "single",
          target = "row",
          selected = 1
        ),
        options = list(
          dom = "lfrtipB",
          buttons = c('copy', 'csv', 'excel'),
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = "700px",
          scrollResize = TRUE,
          deferRender = TRUE,
          autoWidth = TRUE
        )
      )
    })

    render.variance_proportion <- function() {
      res <- get_pcaX()
      shiny::req(res)
      qsee_pcaexplorer_plot_variance_proportion(res)
    }
    render.variance_cumulative <- function() {
      res <- get_pcaX()
      shiny::req(res)
      qsee_pcaexplorer_plot_variance_cumulative(res)
    }

    PlotModuleServer(
      "variance_proportion",
      plotlib = "base",
      func = render.variance_proportion,
      add.watermark = FALSE
    )
    PlotModuleServer(
      "variance_cumulative",
      plotlib = "base",
      func = render.variance_cumulative,
      add.watermark = FALSE
    )

    ## ------------------------------------------------------------------
    ## PCA pairs panel
    ## ------------------------------------------------------------------
    render.scatterpairs <- function() {
      res <- get_pcaX()
      ph <- colorby()
      shiny::req(res, ph)
      qsee_pcaexplorer_plot_scatterpairs(
        res, colorby_var = ph,
        add.ellipse = isTRUE(ellipse())
      )
    }
    render.pc_boxplots <- function() {
      res <- get_pcaX()
      ph <- colorby()
      shiny::req(res, ph)
      qsee_pcaexplorer_plot_pc_boxplots(res, colorby_var = ph, npc = 4)
    }
    render.ftest_vs_pcdim <- function() {
      res <- get_pcaX()
      shiny::req(res)
      graphics::par(mar = c(4, 4, 3, 1))
      qsee_pcaexplorer_plot_ftest_vs_pcdim(res)
    }

    PlotModuleServer(
      "scatterpairs",
      plotlib = "base",
      func = render.scatterpairs,
      add.watermark = FALSE
    )
    PlotModuleServer(
      "pc_boxplots",
      plotlib = "base",
      func = render.pc_boxplots,
      add.watermark = FALSE
    )
    PlotModuleServer(
      "ftest_vs_pcdim",
      plotlib = "base",
      func = render.ftest_vs_pcdim,
      add.watermark = FALSE
    )

    ## ------------------------------------------------------------------
    ## Loadings panel
    ## ------------------------------------------------------------------
    output$loadings <- shiny::renderPlot({
      res <- get_pcaX()
      shiny::req(res)
      qsee_pcaexplorer_plot_loadings_bars(res, npc = 4, ntop = 40)
    })

    ## Full feature × PC loading matrix (res$u)
    output$loadings_table <- shiny::renderDataTable({
      res <- get_pcaX()
      shiny::req(res, res$u)
      U <- as.matrix(res$u)
      if (ncol(U) >= 1L) {
        U <- U[order(-abs(U[, 1])), , drop = FALSE]
      }
      dt <- round(U, 4)

      DT::datatable(
        data = dt,
        class = "compact hover",
        rownames = TRUE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(
          mode = "single",
          target = "row",
          selected = 1
        ),
        options = list(
          dom = "lfrtipB",
          buttons = c('copy', 'csv', 'excel'),
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = "700px",
          scrollResize = TRUE,
          deferRender = TRUE,
          autoWidth = TRUE
        )
      )
    })

    render.loadings_heatmap <- function() {
      res <- get_pcaX()
      shiny::req(res)
      qsee_pcaexplorer_plot_loadings_heatmap(res, npc = 4, ntop = 15)
    }

    PlotModuleServer(
      "loadings_heatmap",
      plotlib = "base",
      func = render.loadings_heatmap,
      add.watermark = FALSE
    )

    ## ------------------------------------------------------------------
    ## Feature PCA panel
    ## ------------------------------------------------------------------
    output$feature_pca <- plotly::renderPlotly({
      redraw_tick()
      res <- get_pcaX()
      normX <- res$X
      shiny::req(res, normX)
      qsee_pcaexplorer_plot_feature_pca(
        res, normX = normX, Y = rYf(),
        numarrows = n_arrows(),
        min_fc = min_fc(),
        var.scale = var_scale()
      )
    })

    output$feature_pca_3d <- plotly::renderPlotly({
      redraw_tick()
      res <- get_pcaX()
      normX <- res$X
      shiny::req(res, normX)
      qsee_pcaexplorer_plot_feature_pca_3d(
        res, normX = normX, Y = rYf(),
        numarrows = n_arrows(),
        min_fc = min_fc(),
        var.scale = var_scale()
      )
    })

    output$feature_pca_all <- shiny::renderPlot({
      res <- get_pcaX()
      normX <- res$X
      shiny::req(res, normX)
      qsee_pcaexplorer_plot_feature_pca_all(
        res, normX = normX, Y = rYf(),
        min_fc = min_fc(),
        var.scale = var_scale(),
        show_arrows = isTRUE(show_arrows())
      )
    })

    render.trait_correlation <- function() {
      res <- get_pcaX()
      shiny::req(res)
      qsee_pcaexplorer_plot_trait_correlation(res, min_fc = min_fc())
    }
    render.loadings2 <- function() {
      res <- get_pcaX()
      shiny::req(res)
      qsee_pcaexplorer_plot_loadings_bars_compact(
        res, npc = 2, ntop = 40,
        normX = res$X, Y = rYf(),
        min_fc = min_fc()
      )
    }

    PlotModuleServer(
      "trait_correlation",
      plotlib = "base",
      func = render.trait_correlation,
      add.watermark = FALSE
    )
    PlotModuleServer(
      "loadings2",
      plotlib = "base",
      func = render.loadings2,
      add.watermark = FALSE
    )
  })
}
