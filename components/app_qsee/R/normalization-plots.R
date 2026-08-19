##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Shared `qsee_plotly_*` helper functions used by all QSEE `_plotly`
## implementations in `plotly-plots.R`. (The original base-R plot functions
## have been removed as they are no longer used.)
##

## ---------------------------------------------------------------------------
## Shared helpers (moved from plotly-plots.R)
## ---------------------------------------------------------------------------

#' Active omicsplots palette.
#' @noRd
qsee_plotly_pal <- function() {
  omicsplots::get_theme_palette()
}

#' Named color vector for a set of discrete levels, from the active theme.
#' @noRd
qsee_plotly_group_colors <- function(levels) {
  levels <- unique(as.character(levels))
  pal <- qsee_plotly_pal()
  stats::setNames(pal@categorical_unlimited(length(levels)), levels)
}

#' Placeholder plot carrying a message (empty / not-computable panels).
#' @noRd
qsee_plotly_empty <- function(msg = "No data") {
  p <- plotly::plot_ly(type = "scatter", mode = "markers")
  p <- plotly::layout(
    p,
    xaxis = list(visible = FALSE),
    yaxis = list(visible = FALSE),
    annotations = list(list(
      text = msg, x = 0.5, y = 0.5, xref = "paper", yref = "paper",
      showarrow = FALSE, font = list(size = 14)
    ))
  )
  omicsplots::plotly_theme(p)
}

#' Drop mixed/empty names that htmlwidgets:::shouldEval rejects.
#'
#' plotly::subplot() can leave `$x$attrs` as a mix of named and unnamed
#' entries (typical when an empty panel is merged with data panels).
#' htmlwidgets then errors with: 'options' must be a fully named list,
#' or have no names (NULL).
#' @noRd
qsee_plotly_fix_names <- function(x) {
  if (!is.list(x) || is.environment(x) || is.function(x)) {
    return(x)
  }
  if (inherits(x, c("data.frame", "POSIXt", "Date", "factor", "formula"))) {
    return(x)
  }
  if (typeof(x) %in% c("closure", "builtin", "special", "externalptr")) {
    return(x)
  }
  nm <- names(x)
  if (!is.null(nm) && length(x) && any(is.na(nm) | !nzchar(nm))) {
    names(x) <- NULL
  }
  for (i in seq_along(x)) {
    xi <- x[[i]]
    if (is.list(xi) && !is.environment(xi) && !is.function(xi)) {
      x[[i]] <- qsee_plotly_fix_names(xi)
    }
  }
  x
}

#' Sanitize a plotly widget for htmlwidgets serialization.
#' @noRd
qsee_plotly_sanitize <- function(p) {
  if (!inherits(p, "htmlwidget") || is.null(p$x)) {
    return(p)
  }
  p$x <- qsee_plotly_fix_names(p$x)
  p
}

#' Fetch an omicsplots internal, or NULL when unavailable.
#' @noRd
qsee_plotly_internal <- function(name) {
  tryCatch(
    utils::getFromNamespace(name, "omicsplots"),
    error = function(e) NULL
  )
}

#' Assemble a named list of plotly panels into a subplot grid.
#'
#' Uses the same assembler as `omicsplots`' own facetting so panel titles,
#' axis theming and legend de-duplication match the rest of the package.
#' Unlike `facet=` on the omicsplots plot functions (which sorts panels
#' alphabetically) this preserves the order of `panels`, which matters here
#' because panels are methods in a meaningful order.
#'
#' @param panels named list of plotly widgets
#' @param n_cols max number of columns
#' @param x_title,y_title shared axis labels (for homogeneous grids)
#' @param margin `plotly::subplot()` margins, `c(left, right, top, bottom)`.
#'   Bump the bottom entry for panels with rotated tick labels.
#' @param axis_titles keep each panel's own axis titles.  Needed for
#'   heterogeneous grids where a shared label makes no sense; the omicsplots
#'   assembler always drops them.
#' @noRd
qsee_plotly_grid <- function(panels, n_cols = 3L,
                             x_title = NULL, y_title = NULL,
                             share_x = FALSE, share_y = FALSE,
                             margin = c(0.02, 0.02, 0.05, 0.05),
                             axis_titles = FALSE) {
  panels <- Filter(Negate(is.null), panels)
  if (!length(panels)) {
    return(qsee_plotly_empty())
  }
  if (length(panels) == 1L) {
    return(qsee_plotly_sanitize(panels[[1]]))
  }

  build_light <- qsee_plotly_internal(".plotly_build_light")
  dedupe <- qsee_plotly_internal(".dedupe_legend_traces")
  assemble <- qsee_plotly_internal(".assemble_subplot")

  panels <- lapply(panels, function(p) {
    p$jsHooks <- NULL
    p$x$dynamic_sample <- NULL
    if (!is.null(build_light)) build_light(p) else plotly::plotly_build(p)
  })
  if (!is.null(dedupe)) panels <- dedupe(panels)

  if (!axis_titles && !is.null(assemble)) {
    return(qsee_plotly_sanitize(assemble(
      panels,
      n_cols = n_cols, share_x = share_x, share_y = share_y,
      margin = margin, x_title = x_title, y_title = y_title,
      build = "standard"
    )))
  }

  ## Own assembly: keeps the per-panel axis titles that .assemble_subplot
  ## strips (also the fallback if the omicsplots internals ever move).
  n_cols <- min(length(panels), n_cols)
  p <- plotly::subplot(
    panels,
    nrows = ceiling(length(panels) / n_cols),
    shareX = share_x, shareY = share_y,
    titleX = axis_titles, titleY = axis_titles, margin = margin
  )
  for (i in seq_along(panels)) {
    ax_x <- if (i == 1) "xaxis" else paste0("xaxis", i)
    ax_y <- if (i == 1) "yaxis" else paste0("yaxis", i)
    dx <- p$x$layout[[ax_x]]$domain
    dy <- p$x$layout[[ax_y]]$domain
    if (is.null(dx)) dx <- c(0, 1)
    if (is.null(dy)) dy <- c(0, 1)
    p <- plotly::add_annotations(
      p,
      text = paste0("<b>", names(panels)[i], "</b>"),
      x = mean(dx), y = dy[2], xref = "paper", yref = "paper",
      xanchor = "center", showarrow = FALSE, yshift = 10,
      font = list(size = 13)
    )
  }
  qsee_plotly_sanitize(
    plotly::layout(p, margin = list(l = 60, r = 20, t = 40, b = 60))
  )
}

#' Add a text-label layer to a plotly panel.
#'
#' `label_names=` on `pgx.plot_scatter()` emits layout annotations, which
#' `plotly::subplot()` collapses into a single vectorised annotation that
#' plotly.js renders as one stray "new text" label.  A text *trace* survives
#' subplotting, so labels are added that way instead.
#' @noRd
qsee_plotly_add_labels <- function(p, x, y, labels, nmax = 60) {
  if (is.null(labels) || !length(labels) || length(labels) > nmax) {
    return(p)
  }
  plotly::add_trace(
    p,
    x = x, y = y, type = "scatter", mode = "text", inherit = FALSE,
    text = labels, textposition = "top center",
    textfont = list(size = 10, color = qsee_plotly_pal()@chrome$foreground),
    showlegend = FALSE, hoverinfo = "skip"
  )
}

#' Bar plot that preserves the order of the categories.
#'
#' `omicsplots::pgx.plot_barplot()` re-sorts categories by value, which is
#' the right default for ranked gene/pathway summaries but wrong for the
#' QC plots here where the x axis is samples, principal components or
#' intensity bins.  This keeps the incoming order and otherwise mirrors the
#' omicsplots look (theme palette + `plotly_theme()`).
#'
#' @noRd
qsee_plotly_barplot <- function(x, y, fill = NULL, stacked = TRUE,
                                xlab = NULL, ylab = NULL, color = NULL,
                                showlegend = TRUE, hline = NULL, ylim = NULL) {
  x <- as.character(x)
  lev <- unique(x)
  pal <- qsee_plotly_pal()
  line <- list(color = pal@chrome$foreground, width = 0.3)
  p <- plotly::plot_ly()

  if (is.null(fill)) {
    col <- if (is.null(color)) pal@single else unname(color)[1]
    p <- plotly::add_trace(
      p,
      x = x, y = y, type = "bar",
      marker = list(color = col, line = line),
      hoverinfo = "text", textposition = "none",
      text = paste0("<b>", x, "</b><br>", signif(y, 4)),
      showlegend = FALSE
    )
  } else {
    fill <- as.character(fill)
    glev <- unique(fill)
    cmap <- if (is.null(color)) qsee_plotly_group_colors(glev) else color
    for (g in glev) {
      i <- which(fill == g)
      p <- plotly::add_trace(
        p,
        x = x[i], y = y[i], type = "bar", name = g, legendgroup = g,
        marker = list(color = unname(cmap[[g]]), line = line),
        hoverinfo = "text", textposition = "none",
        text = paste0("<b>", x[i], "</b><br>", g, ": ", signif(y[i], 4)),
        showlegend = showlegend
      )
    }
  }

  yaxis <- list()
  if (!is.null(ylab)) yaxis$title <- list(text = ylab)
  if (!is.null(ylim)) yaxis$range <- ylim
  xaxis <- list(type = "category", categoryorder = "array", categoryarray = lev)
  if (!is.null(xlab)) xaxis$title <- list(text = xlab)
  p <- plotly::layout(
    p,
    barmode = if (stacked) "stack" else "group",
    xaxis = xaxis,
    yaxis = yaxis
  )

  ## Reference lines are drawn as traces (not layout shapes) so they survive
  ## plotly::subplot() when the panel is put in a grid.
  for (h in hline) {
    p <- plotly::add_trace(
      p,
      x = c(lev[1], lev[length(lev)]), y = c(h, h),
      type = "scatter", mode = "lines", inherit = FALSE,
      line = list(color = "#bf616a", dash = "dot", width = 1),
      showlegend = FALSE, hoverinfo = "skip"
    )
  }
  omicsplots::plotly_theme(p)
}

#' Add a dashed vertical / horizontal reference line as a trace (unused in current plots but kept for completeness).
#' @noRd
qsee_plotly_refline <- function(p, v = NULL, h = NULL, xrange = NULL, yrange = NULL,
                                color = "#bf616a") {
  ln <- list(color = color, dash = "dash", width = 1)
  for (x0 in v) {
    p <- plotly::add_trace(
      p,
      x = c(x0, x0), y = yrange, type = "scatter", mode = "lines",
      inherit = FALSE, line = ln, showlegend = FALSE, hoverinfo = "skip"
    )
  }
  for (y0 in h) {
    p <- plotly::add_trace(
      p,
      x = xrange, y = c(y0, y0), type = "scatter", mode = "lines",
      inherit = FALSE, line = ln, showlegend = FALSE, hoverinfo = "skip"
    )
  }
  p
}

#' Grid of iheatmapr panels.
#'
#' iheatmapr widgets cannot be assembled with `plotly::subplot()`, so
#' multi-method heatmaps are laid out as `n` separate outputs in a CSS grid.
#' Unused slots simply stay blank.
#'
#' @param ns module namespace function
#' @param prefix output-id prefix, must match the server side
#' @param n number of slots to reserve
#' @noRd
qsee_plotly_hm_grid_ui <- function(ns, prefix, n = 6L, ncol = 3L,
                                   height = "330px") {
  cells <- lapply(seq_len(n), function(i) {
    shiny::div(
      style = "min-width: 0;",
      shiny::div(
        style = paste(
          "text-align:center; font-weight:600; font-size:13px;",
          "padding-bottom:2px; height:18px; overflow:hidden;"
        ),
        shiny::textOutput(ns(paste0(prefix, "_title_", i)), inline = TRUE)
      ),
      iheatmapr::iheatmaprOutput(ns(paste0(prefix, "_hm_", i)), height = height)
    )
  })
  shiny::div(
    style = sprintf(
      "display:grid; grid-template-columns: repeat(%d, minmax(0, 1fr)); gap:10px;",
      ncol
    ),
    cells
  )
}

#' Render a named list of iheatmapr widgets into a fixed set of grid slots.
#' @noRd
qsee_plotly_hm_grid_server <- function(output, prefix, panels, n = 6L) {
  for (i in seq_len(n)) {
    local({
      ii <- i
      output[[paste0(prefix, "_title_", ii)]] <- shiny::renderText({
        p <- panels()
        if (is.null(p) || length(p) < ii) "" else names(p)[ii]
      })
      output[[paste0(prefix, "_hm_", ii)]] <- htmlwidgets::shinyRenderWidget(
        quote({
          p <- panels()
          shiny::req(!is.null(p), length(p) >= ii)
          w <- p[[ii]]
          ## pgx.plot_heatmap() normally returns an already-converted
          ## "iheatmapr" widget, but accept a raw S4 Iheatmap too.
          if (methods::is(w, "Iheatmap")) iheatmapr::to_widget(w) else w
        }),
        iheatmapr::iheatmaprOutput,
        environment(),
        quoted = TRUE
      )
    })
  }
}

#' Shorthand for a full-height plotly panel inside a card.
#' @noRd
qsee_plotly_card <- function(ns, id, height = "720px") {
  bslib::card(
    full_screen = TRUE,
    plotly::plotlyOutput(ns(id), height = height)
  )
}

#' Render a single iheatmapr widget output, accepting either an already
#' converted htmlwidget or a raw S4 Iheatmap.
#' @noRd
qsee_plotly_hm_server <- function(output, id, widget_expr_fn) {
  output[[id]] <- htmlwidgets::shinyRenderWidget(
    quote({
      w <- widget_expr_fn()
      if (methods::is(w, "Iheatmap")) iheatmapr::to_widget(w) else w
    }),
    iheatmapr::iheatmaprOutput,
    environment(),
    quoted = TRUE
  )
}

#' Standard PCA-style scatter panel used by several boards.
#'
#' Mirrors `plot(pos, pch = 19, col = cc); text(pos, labels, pos = 3)`.
#' @noRd
qsee_plotly_pca_panel <- function(pos, color, xlab = "PC1", ylab = "PC2",
                                  label_nmax = 60, show_labels = TRUE) {
  labels <- rownames(pos)
  if (is.null(labels)) labels <- paste0("sample", seq_len(nrow(pos)))
  p <- omicsplots::pgx.plot_scatter(
    x = pos[, 1], y = pos[, 2],
    labels = if (show_labels) labels else NULL,
    color = color,
    xlab = xlab, ylab = ylab
  )
  if (show_labels) {
    p <- qsee_plotly_add_labels(p, pos[, 1], pos[, 2], labels, nmax = label_nmax)
  }
  p
}

#' Quadrant guide reproducing the legend panel of bc.CovariateAnalysisPlot().
#' @noRd
qsee_plotly_quadrant_legend <- function() {
  labs <- data.frame(
    x = c(0.2, 0.75, 0.2, 0.75),
    y = c(0.80, 0.80, 0.20, 0.20),
    t = c(
      "strong batch-effect", "strong model-parameter",
      "nuisance parameter", "weak model-parameter"
    ),
    stringsAsFactors = FALSE
  )
  fg <- qsee_plotly_pal()@chrome$foreground
  ln <- list(color = fg, dash = "dash", width = 1)
  p <- plotly::plot_ly()
  p <- plotly::add_trace(p,
    x = c(0.5, 0.5), y = c(0, 1), type = "scatter",
    mode = "lines", line = ln, showlegend = FALSE, hoverinfo = "skip"
  )
  p <- plotly::add_trace(p,
    x = c(0, 1), y = c(0.5, 0.5), type = "scatter",
    mode = "lines", line = ln, showlegend = FALSE, hoverinfo = "skip"
  )
  p <- plotly::add_trace(
    p,
    x = labs$x, y = labs$y, type = "scatter", mode = "text",
    text = labs$t, textposition = "middle center",
    textfont = list(size = 12, color = fg),
    showlegend = FALSE, hoverinfo = "skip"
  )
  p <- plotly::layout(
    p,
    xaxis = list(
      range = c(0, 1), showticklabels = FALSE, zeroline = FALSE,
      title = list(text = "correlation with phenotype")
    ),
    yaxis = list(
      range = c(0, 1), showticklabels = FALSE, zeroline = FALSE,
      title = list(text = "correlation with PC")
    )
  )
  omicsplots::plotly_theme(p)
}

## ---------------------------------------------------------------------------
## Normalization board (plotly versions)
## ---------------------------------------------------------------------------

#' Plotly version of [qsee_normalization_plot_boxplots()].
#'
#' @param nmax max number of samples (boxes) to draw
#' @param feature_nmax max number of features fed to each box.  plotly ships
#'   every observation to the browser, so the full matrix would produce a
#'   multi-MB widget; a random subset of features gives the same box.
qsee_normalization_plot_boxplots_plotly <- function(res, Y, ph, nmax = 40,
                                                    feature_nmax = 1000) {
  normX <- res$normX
  y <- factor(Y[, ph])
  ## One feature subset shared by all panels, so the boxes stay comparable.
  n1 <- nrow(normX[[1]])
  ii <- if (n1 > feature_nmax) sort(sample(n1, feature_nmax)) else seq_len(n1)
  panels <- list()
  for (name in names(normX)) {
    X <- normX[[name]]
    if (is.null(X) || !length(X)) {
      panels[[name]] <- qsee_plotly_empty("No data")
      next
    }
    jj <- seq_len(ncol(X))
    if (ncol(X) > nmax) jj <- sort(sample(ncol(X), nmax))
    X <- X[ii, jj, drop = FALSE]
    yj <- as.character(y)[jj]
    cmap <- qsee_plotly_group_colors(levels(y))
    ## pgx.plot_boxplot() maps group -> color; one box per sample, colored
    ## by that sample's phenotype (same as col=cc in the base version).
    sample_cols <- stats::setNames(unname(cmap[yj]), colnames(X))
    v <- as.vector(X)
    g <- rep(colnames(X), each = nrow(X))
    keep <- is.finite(v)
    bp <- omicsplots::pgx.plot_boxplot(
      values = v[keep],
      groups = factor(g[keep], levels = colnames(X)),
      type = "box", color = sample_cols, ylab = "log-expression"
    )
    ## pgx.plot_boxplot() hardcodes boxpoints = "all", which hides the box
    ## itself for matrices this size. Show outliers only, like boxplot().
    ## The legend would be one entry per sample; the base plot has none.
    panels[[name]] <- plotly::layout(
      plotly::style(bp, boxpoints = "outliers"),
      showlegend = FALSE,
      xaxis = list(categoryorder = "array", categoryarray = colnames(X))
    )
  }
  qsee_plotly_grid(panels,
    n_cols = 3, y_title = "log-expression",
    margin = c(0.02, 0.02, 0.10, 0.15)
  )
}

#' Plotly version of [qsee_normalization_plot_histograms()].
#'
#' The base version overlays per-sample histogram outlines; the plotly
#' equivalent in omicsplots is the per-sample density chart.
qsee_normalization_plot_histograms_plotly <- function(res, Y, ph, nmax = 40) {
  normX <- res$normX
  y <- factor(Y[, ph])
  cmap <- qsee_plotly_group_colors(levels(y))
  panels <- list()
  for (name in names(normX)) {
    X <- normX[[name]]
    if (is.null(X) || !length(X)) {
      panels[[name]] <- qsee_plotly_empty("No data")
      next
    }
    jj <- seq_len(ncol(X))
    if (ncol(X) > nmax) jj <- sort(sample(ncol(X), nmax))
    X <- X[, jj, drop = FALSE]
    p <- omicsplots::pgx.plot_density(
      X,
      groups = as.character(y)[jj], color = cmap,
      n_breaks = 80, xlab = "signal", showlegend = FALSE
    )
    ## Only show tooltip for the exact trace/group under the cursor
    ## (default hovermode = "x" shows all groups at that x-position)
    panels[[name]] <- plotly::layout(p, hovermode = "closest")
  }
  qsee_plotly_grid(panels,
    n_cols = 3,
    margin = c(0.04, 0.04, 0.07, 0.07),
    x_title = "signal", y_title = "Density"
  )
}

#' Plotly version of [qsee_normalization_plot_pca()].
qsee_normalization_plot_pca_plotly <- function(res, Y, ph, show_labels = TRUE) {
  pcaX <- Filter(Negate(is.null), res$pcaX)
  y <- factor(Y[, ph])
  panels <- lapply(names(pcaX), function(name) {
    qsee_plotly_pca_panel(pcaX[[name]], color = y, show_labels = show_labels)
  })
  names(panels) <- names(pcaX)
  qsee_plotly_grid(panels,
    n_cols = 3,
    margin = c(0.04, 0.04, 0.07, 0.07),
    x_title = "PC1", y_title = "PC2"
  )
}
