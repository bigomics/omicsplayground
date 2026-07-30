## Biplot grid for the BC panel: one biplot per method (sample scores as
## points, arrows overlaid) on the two selected PCs. Same method list /
## grid as plot_all_methods. Two arrow sources:
##  - "loadings": top feature loadings (res$loadings, same PCA as scores)
##  - "pheno": each annotation's correlation with the two PCs, precomputed
##    in playbase (res$pheno.cor; eigencorplot-style, plain cor())
pcaexplorer.plot_biplot <- function(res, samples, xpc=1, ypc=2, colorby_var=NULL,
                                    main = "PCA biplot", cex = 2, cex.labels = 0.85,
                                    arrows = "loadings", arrow.col = "#B3444488",
                                    numarrows = 10, arrow.lwd = 2,
                                    var.scale = 0.5,
                                    legend = TRUE, add.ellipse=TRUE ) {

  if(0) {

    X <- playbase::logCPM(playbase::COUNTS)
    X <- X - rowMeans(X,na.rm=TRUE)
    Y <- playbase::SAMPLES
    res <- svd(X)
    rownames(res$v) <- rownames(Y)
    rownames(res$u) <- rownames(X)
    H <- playbase::expandPhenoMatrix(Y, drop.ref=FALSE)
    res$R <- cor(H, res$v)
    xpc=1;ypc=2
  }
  
  ## arrow endpoints for a method: named 2-col matrix in loading- or
  ## correlation-units (scaled to the score cloud later). Both come
  ## precomputed from playbase::compare_batchcorrection_methods.
  arrow_ends <- function(res, scores, arrows="loadings") {
    if (arrows == "loadings") {
      L <- res$u
      if (is.null(L) || max(xpc, ypc) > ncol(L)) {
        return(NULL)
      }
      A <- L[, c(xpc, ypc), drop = FALSE]
    } else {
      Rc <- res$R
      if (is.null(Rc) || max(xpc, ypc) > ncol(Rc)) {
        return(NULL)
      }
      A <- Rc[, c(xpc, ypc), drop = FALSE]
      if (!nrow(A)) {
        return(NULL)
      }
    }
    A <- A[stats::complete.cases(A), , drop = FALSE]
    A[head(order(-(A[, 1]^2 + A[, 2]^2)), numarrows), , drop = FALSE]
  }
  
  P <- res$v
  if (is.null(P) || max(xpc, ypc) > ncol(P)) {
    plot.new()
    text(0.45, 0.5, "method failed")
    title("PCA biplot", cex.main = 1.3, adj = 0)
    return(invisible())
  }
  varx <- round(100 * res$d**2 / sum(res$d**2),1)
  scores <- P[, c(xpc, ypc), drop = FALSE]
  varx <- varx[c(xpc, ypc)]
  
  smp <- samples[rownames(scores), , drop = FALSE]
  group <- factor(smp[, colorby_var])
  color <- 1 + as.integer(group)
  names(color) <- group
  cex1 <- cut(nrow(scores), c(0, 40, 100, 250, 1000, 999999),
    c(1, 0.85, 0.7, 0.55, 0.4))
  cex1 <- 2 * cex * as.numeric(as.character(cex1))
  if (is.na(cex1)) cex1 <- 1
  
  lim <- c(-1, 1) * max(abs(scores)) * 1.2
  plot(scores,
    col = color, pch = 20, cex = cex1, las = 1, xlim = lim, ylim = lim,
    xlab = paste0("PC", xpc, " (", varx[1],"%)"),
    ylab = paste0("PC", ypc, " (", varx[2],"%)"),
  )
  title(main, cex.main = 1.3)
  abline(h = 0, v = 0, lty = 3, col = "grey70")

  plotrix::thigmophobe.labels(
    scores[,1], scores[,2], rownames(scores),
    cex = cex.labels, col = color, offset = 0.4, xpd = NA
  )

  if(numarrows > 0) {
    A <- arrow_ends(res, scores, arrows = arrows)
    if (!is.null(A) && nrow(A) && max(abs(A)) > 0) {
      ## scale arrows to fill the score cloud (loadings/correlations sit on a
      ## much smaller scale than the scores)
      s <- var.scale * max(abs(scores)) / max(abs(A))
      ax <- A[, 1] * s
      ay <- A[, 2] * s
      arrows(0, 0, ax, ay, length = 0.15, col = arrow.col, lwd = arrow.lwd)
      ## repel labels off each other (base-graphics; places each away from
      ## its nearest neighbour) rather than a fixed left/right offset
      plotrix::thigmophobe.labels(ax, ay, rownames(A),
        cex = cex.labels, col = arrow.col, offset = 0.4, xpd = NA
      )
    }
  }
  
  if(add.ellipse) {
    g = group[1]
    for(g in unique(group)) {
      ii <- which(group == g)
      mat <- scores[ii,,drop=FALSE]
      if(nrow(mat) < 2) next
      
      # 3. Calculate ellipse boundary from the matrix covariance
      # (t = radius multiplier; t = 2 captures roughly 95% of data for bivariate normal)
      ellipse_points <- ellipse::ellipse(
        cov(mat), centre = colMeans(mat), t = 2)
      
      # 4. Layer the ellipse onto the plot
      ##lines(ellipse_points, col = color[g], lwd = 2)
      c1 <- adjustcolor( color[g], alpha.f = 0.12)
      polygon(ellipse_points, col = c1, border = color[g], lwd = 1.4)
    }    
  }

  if(legend) {
    cc <- color[levels(group)]
    legend("bottomright", levels(group), fill=cc, cex=0.85, y.intersp=0.85)
  }
  
}



## 2D plotly counterpart of pcaexplorer.plot_biplot: same score/arrow logic
## (top loadings or top pheno-correlations, scaled to fill the score cloud),
## rendered as an interactive plotly scatter. Optional group ellipses as
## closed line traces (same covariance ellipse as the base plot).
pcaexplorer.plot_biplot_2d <- function(res, samples, xpc=1, ypc=2, colorby_var=NULL,
                                       main = "PCA biplot", arrows = "loadings",
                                       arrow.col = "black", numarrows = 10,
                                       var.scale = 0.65, marker.size = 10,
                                       legend = TRUE, add.ellipse = TRUE) {

  arrow_ends <- function(res, arrows = "loadings") {
    A <- if (arrows == "loadings") res$u else res$R
    if (is.null(A) || max(xpc, ypc) > ncol(A)) {
      return(NULL)
    }
    A <- A[, c(xpc, ypc), drop = FALSE]
    A <- A[stats::complete.cases(A), , drop = FALSE]
    if (!nrow(A)) {
      return(NULL)
    }
    A[head(order(-(A[, 1]^2 + A[, 2]^2)), numarrows), , drop = FALSE]
  }

  P <- res$v
  if (is.null(P) || max(xpc, ypc) > ncol(P)) {
    return(plotly::plotly_empty(type = "scatter") %>%
      plotly::layout(title = "method failed"))
  }

  varx <- round(100 * res$d**2 / sum(res$d**2), 1)
  scores <- P[, c(xpc, ypc), drop = FALSE]
  varx <- varx[c(xpc, ypc)]

  smp <- samples[rownames(scores), , drop = FALSE]
  group <- factor(smp[, colorby_var])
  clrs <- grDevices::palette()[1 + seq_along(levels(group))]
  names(clrs) <- levels(group)

  df <- data.frame(
    x = scores[, 1], y = scores[, 2],
    sample = rownames(scores), group = group
  )

  lim <- c(-1, 1) * max(abs(scores)) * 1.2
  df$color <- clrs[as.character(df$group)]

  fig <- plotly::plot_ly(df) %>%
    plotly::add_markers(
      x = ~x, y = ~y,
      type = "scatter", mode = "markers",
      color = ~group, colors = clrs,
      marker = list(size = marker.size, line = list(color = "grey20", width = 0.5)),
      text = ~sample, hoverinfo = "text",
      showlegend = legend
    )

  ## sample labels: same color as their group points (one trace per group
  ## so textfont color is not forced to a single value)
  for (g in levels(group)) {
    dfg <- df[df$group == g, , drop = FALSE]
    if (!nrow(dfg)) next
    fig <- fig %>% plotly::add_text(
      data = dfg, x = ~x, y = ~y, text = ~sample,
      type = "scatter", mode = "text",
      textposition = "top center",
      textfont = list(size = 14, color = clrs[[g]]),
      showlegend = FALSE, hoverinfo = "skip", inherit = FALSE
    )
  }

  if (isTRUE(add.ellipse)) {
    for (g in levels(group)) {
      ii <- which(group == g)
      mat <- scores[ii, , drop = FALSE]
      if (nrow(mat) < 2) next
      ep <- tryCatch(
        ellipse::ellipse(cov(mat), centre = colMeans(mat), t = 2),
        error = function(e) NULL
      )
      if (is.null(ep)) next
      fig <- fig %>% plotly::add_paths(
        x = ep[, 1], y = ep[, 2],
        type = "scatter", mode = "lines",
        line = list(color = clrs[[g]], width = 1.2),
        fill = "toself",
        fillcolor = grDevices::adjustcolor(clrs[[g]], alpha.f = 0.12),
        name = g, legendgroup = g,
        showlegend = FALSE, hoverinfo = "skip", inherit = FALSE
      )
    }
  }

  ## Arrows via annotations (empty text so the label is not stuck at the
  ## tail/origin) + separate tip labels (showarrow=FALSE) so names sit at
  ## the arrow heads.
  arrow_ann <- list()
  if (numarrows > 0) {
    A <- arrow_ends(res, arrows = arrows)
    if (!is.null(A) && nrow(A) && max(abs(A)) > 0) {
      s <- var.scale * max(abs(scores)) / max(abs(A))
      ax <- A[, 1] * s
      ay <- A[, 2] * s
      r <- sqrt(ax^2 + ay^2)
      r[r == 0] <- 1
      nudge <- 0.05 * max(abs(scores))
      lx <- ax + nudge * ax / r
      ly <- ay + nudge * ay / r
      arrow_ann <- c(
        ## shaft + arrowhead from origin to tip
        lapply(seq_len(nrow(A)), function(i) {
          list(
            x = ax[i], y = ay[i], ax = 0, ay = 0,
            xref = "x", yref = "y", axref = "x", ayref = "y",
            showarrow = TRUE, text = "",
            arrowhead = 3, arrowsize = 1.2, arrowwidth = 1.8,
            arrowcolor = arrow.col
          )
        }),
        ## variable label just beyond the tip
        lapply(seq_len(nrow(A)), function(i) {
          list(
            x = lx[i], y = ly[i],
            xref = "x", yref = "y",
            showarrow = FALSE,
            text = rownames(A)[i],
            xanchor = if (ax[i] >= 0) "left" else "right",
            yanchor = if (ay[i] >= 0) "bottom" else "top",
            font = list(color = arrow.col, size = 15),
            bgcolor = "rgba(255,255,255,0.65)",
            borderpad = 2
          )
        })
      )
    }
  }

  ## simple_white-like theme (named templates are a no-op on plotly 4.12's
  ## bundled plotly.js). White panel, light grid, axis lines + full box.
  axis_sw <- function(title_txt, ...) {
    utils::modifyList(
      list(
        title = list(text = title_txt, standoff = 8),
        showgrid = TRUE, gridcolor = "rgb(230,230,230)", gridwidth = 1,
        zeroline = TRUE, zerolinecolor = "rgb(200,200,200)", zerolinewidth = 1,
        showline = TRUE, linecolor = "black", linewidth = 1,
        mirror = TRUE,  ## draw opposite spine -> full box with both axes
        ticks = "outside", tickcolor = "black", ticklen = 4,
        automargin = TRUE
      ),
      list(...)
    )
  }
  fig %>% plotly::layout(
    title = list(text = main, font = list(size = 20, color = "black")),
    font = list(family = "Arial, Helvetica, sans-serif", size = 16, color = "black"),
    paper_bgcolor = "white",
    plot_bgcolor = "white",
    xaxis = axis_sw(
      paste0("PC", xpc, " (", varx[1], "%)"),
      range = lim, scaleanchor = "y", scaleratio = 1
    ),
    yaxis = axis_sw(
      paste0("PC", ypc, " (", varx[2], "%)"),
      range = lim
    ),
    annotations = arrow_ann,
    showlegend = legend,
    ## horizontal legend under the plot; entrywidth forces wrapping into
    ## multiple rows when there are many groups
    legend = list(
      orientation = "h",
      x = 0.5, xanchor = "center",
      y = -0.12, yanchor = "top",
      entrywidth = 80,
      entrywidthmode = "pixels",
      traceorder = "normal",
      itemsizing = "constant",
      bgcolor = "rgba(255,255,255,0.7)"
    ),
    margin = list(b = 100, l = 60, r = 20, t = 50)
  )
}


## 3D counterpart of pcaexplorer.plot_biplot: same score/arrow logic (top
## loadings or top pheno-correlations, scaled to fill the score cloud), but
## rendered as an interactive plotly scatter3d instead of base graphics.
## Optional group ellipsoids as scatter3d wireframes (t=2 covariance).
pcaexplorer.plot_biplot_3d <- function(res, samples, xpc=1, ypc=2, zpc=3, colorby_var=NULL,
                                       main = "PCA biplot", arrows = "loadings",
                                       arrow.col = "black", numarrows = 10,
                                       var.scale = 0.5, marker.size = 6, legend = TRUE,
                                       add.ellipse = TRUE) {

  arrow_ends <- function(res, arrows="loadings") {
    A <- if (arrows == "loadings") res$u else res$R
    if (is.null(A) || max(xpc, ypc, zpc) > ncol(A)) {
      return(NULL)
    }
    A <- A[, c(xpc, ypc, zpc), drop = FALSE]
    A <- A[stats::complete.cases(A), , drop = FALSE]
    if (!nrow(A)) {
      return(NULL)
    }
    A[head(order(-(A[, 1]^2 + A[, 2]^2 + A[, 3]^2)), numarrows), , drop = FALSE]
  }

  ## Covariance wireframe for scatter3d lines (NA-separated segments).
  ## Full rank-3 -> ellipsoid (meridians+parallels). Rank-2 (e.g. 3 coplanar
  ## points) -> thin ellipsoid: missing axis inflated to thin.frac of the
  ## geometric-mean valid radius (default 20%). Rank-1 -> line segment. t=2 as 2D.
  ellipsoid_wireframe <- function(mu, S, t = 2, n_u = 24, n_v = 12,
                                  thin.frac = 0.10) {
    eig <- tryCatch(eigen(S, symmetric = TRUE), error = function(e) NULL)
    if (is.null(eig) || any(!is.finite(eig$values))) {
      return(NULL)
    }
    vals <- pmax(as.numeric(eig$values), 0)
    pos <- which(vals > max(vals, 0) * 1e-8 + 1e-14)
    if (!length(pos)) {
      return(NULL)
    }
    ## Rank-2: inflate the zero axis so we draw a thin 3D ellipsoid instead
    ## of a flat circle. Radius_thin = thin.frac * geom.mean(valid radii).
    if (length(pos) == 2L && length(vals) >= 3L) {
      r_pos <- sqrt(vals[pos])
      r_thin <- thin.frac * exp(mean(log(pmax(r_pos, .Machine$double.eps))))
      zero <- setdiff(seq_along(vals), pos)[1L]
      vals[zero] <- r_thin^2
      pos <- c(pos, zero)
    }
    na3 <- c(NA_real_, NA_real_, NA_real_)
    if (length(pos) >= 3) {
      R <- eig$vectors[, pos[1:3], drop = FALSE] %*%
        diag(t * sqrt(vals[pos[1:3]]), 3)
      u <- seq(0, 2 * pi, length.out = n_u)
      v <- seq(0, pi, length.out = n_v)
      segs <- vector("list", length(v) + length(u))
      k <- 0L
      for (vj in v) {
        xyz <- cbind(cos(u) * sin(vj), sin(u) * sin(vj), rep(cos(vj), length(u)))
        pts <- sweep(xyz %*% t(R), 2, mu, "+")
        k <- k + 1L
        segs[[k]] <- rbind(pts, na3)
      }
      for (ui in u) {
        xyz <- cbind(rep(cos(ui), length(v)) * sin(v),
                     rep(sin(ui), length(v)) * sin(v),
                     cos(v))
        pts <- sweep(xyz %*% t(R), 2, mu, "+")
        k <- k + 1L
        segs[[k]] <- rbind(pts, na3)
      }
      return(do.call(rbind, segs))
    }
    ## rank-1: segment along the single axis
    v1 <- eig$vectors[, pos[1]]
    r <- t * sqrt(vals[pos[1]])
    rbind(mu - r * v1, mu + r * v1, na3)
  }

  P <- res$v
  if (is.null(P) || max(xpc, ypc, zpc) > ncol(P)) {
    return(plotly::plotly_empty(type = "scatter3d") %>%
      plotly::layout(title = "method failed"))
  }

  varx <- round(100 * res$d**2 / sum(res$d**2), 1)
  scores <- P[, c(xpc, ypc, zpc), drop = FALSE]
  varx <- varx[c(xpc, ypc, zpc)]

  smp <- samples[rownames(scores), , drop = FALSE]
  group <- factor(smp[, colorby_var])
  clrs <- grDevices::palette()[1 + seq_along(levels(group))]
  names(clrs) <- levels(group)

  df <- data.frame(
    x = scores[, 1], y = scores[, 2], z = scores[, 3],
    sample = rownames(scores), group = group
  )

  fig <- plotly::plot_ly(df, mode = "markers+text") %>%
    plotly::add_markers(
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d",
      color = ~group, colors = clrs,
      marker = list(size = marker.size, line = list(color = "grey20", width = 0.3)),
      text = ~sample, hoverinfo = "text",
      textposition = "top center"
    )

  ## covariance ellipsoid / planar-ellipse wireframes (one lines trace / group)
  if (isTRUE(add.ellipse)) {
    for (g in levels(group)) {
      ii <- which(group == g)
      mat <- scores[ii, , drop = FALSE]
      if (nrow(mat) < 2) next
      S <- tryCatch(stats::cov(mat), error = function(e) NULL)
      if (is.null(S) || anyNA(S)) next
      wf <- ellipsoid_wireframe(colMeans(mat), S, t = 2)
      if (is.null(wf) || !nrow(wf)) next
      ## rgba line color is more reliable than opacity alone for scatter3d
      line_col <- grDevices::adjustcolor(clrs[[g]], alpha.f = 0.8)
      fig <- fig %>% plotly::add_trace(
        x = wf[, 1], y = wf[, 2], z = wf[, 3],
        type = "scatter3d", mode = "lines",
        line = list(color = line_col, width = 2),
        opacity = 0.60,
        name = g, legendgroup = g,
        showlegend = FALSE, hoverinfo = "skip", inherit = FALSE
      )
    }
  }

  if (numarrows > 0) {
    A <- arrow_ends(res, arrows = arrows)
    if (!is.null(A) && nrow(A) && max(abs(A)) > 0) {
      ## scale arrows to fill the score cloud, same as the 2D biplot
      s <- var.scale * max(abs(scores)) / max(abs(A))
      ax <- A[, 1] * s
      ay <- A[, 2] * s
      az <- A[, 3] * s
      for (i in seq_len(nrow(A))) {
        fig <- fig %>% plotly::add_trace(
          x = c(0, ax[i]), y = c(0, ay[i]), z = c(0, az[i]),
          type = "scatter3d", mode = "lines",
          line = list(color = arrow.col, width = 4),
          showlegend = FALSE, hoverinfo = "skip", inherit = FALSE
        )
      }
      fig <- fig %>% plotly::add_text(
        x = ax, y = ay, z = az, text = rownames(A),
        type = "scatter3d", mode = "text",
        textfont = list(color = arrow.col, size = 14),
        showlegend = FALSE, hoverinfo = "skip", inherit = FALSE
      )
    }
  }

  fig %>% plotly::layout(
    title = main,
    scene = list(
      xaxis = list(title = paste0("PC", xpc, " (", varx[1], "%)")),
      yaxis = list(title = paste0("PC", ypc, " (", varx[2], "%)")),
      zaxis = list(title = paste0("PC", zpc, " (", varx[3], "%)")),
      aspectmode = "cube"
    ),
    showlegend = legend
  )
}


pcaexplorer.plot_scatterpairs <- function(res, colorby_var, add.ellipse) {
  group <- factor(res$Y[rownames(res$v), colorby_var])
  cc <- 1 + as.integer(group)
  V <- res$v[, 1:min(4L, ncol(res$v)), drop = FALSE]
  colnames(V) <- paste0("PC", seq_len(ncol(V)))
  npc <- ncol(V)

  ## mean silhouette width for every unordered PC pair, using euclidean
  ## distance in that 2D PC plane and the colorby labels as clusters.
  ## Singleton groups (size < 2) are dropped; score is computed on the
  ## remaining groups that each have size >= 2 (needs >= 2 such groups).
  sil_pair <- function(i, j) {
    ok <- !is.na(group)
    g0 <- droplevels(group[ok])
    xy0 <- V[ok, c(i, j), drop = FALSE]
    ## keep only groups large enough for silhouette
    tab <- table(g0)
    keep_lv <- names(tab)[tab >= 2L]
    if (length(keep_lv) < 2L) {
      return(NA_real_)
    }
    keep <- g0 %in% keep_lv
    g_ok <- droplevels(g0[keep])
    xy <- xy0[keep, , drop = FALSE]
    if (nrow(xy) < 3L) {
      return(NA_real_)
    }
    d <- stats::dist(xy, method = "euclidean")
    sil <- tryCatch(
      cluster::silhouette(as.integer(g_ok), d),
      error = function(e) NULL
    )
    if (is.null(sil)) {
      return(NA_real_)
    }
    mean(sil[, "sil_width"], na.rm = TRUE)
  }
  sil_mat <- matrix(NA_real_, npc, npc,
    dimnames = list(colnames(V), colnames(V)))
  if (npc >= 2L) {
    for (i in seq_len(npc - 1L)) {
      for (j in (i + 1L):npc) {
        s <- sil_pair(i, j)
        sil_mat[i, j] <- sil_mat[j, i] <- s
      }
    }
  }

  ## map a panel coordinate vector back to its PC index (pairs passes
  ## as.vector(V[, k]), so match on values)
  pc_index <- function(z) {
    for (k in seq_len(npc)) {
      if (isTRUE(all.equal(as.numeric(z), as.numeric(V[, k]),
            tolerance = 1e-10, check.attributes = FALSE))) {
        return(k)
      }
    }
    NA_integer_
  }

  par(cex = 1.4)
  ## upper triangle: scatter (+ optional ellipses)
  panel_upper <- function(x, y, col, ...) {
    points(x, y, col = col, ...)
    if (isTRUE(add.ellipse)) {
      pts <- cbind(x, y)
      for (c1 in unique(col)) {
        ii <- which(col == c1)
        if (length(ii) < 2) next
        mat <- pts[ii, , drop = FALSE]
        ep <- tryCatch(
          ellipse::ellipse(stats::cov(mat), centre = colMeans(mat), t = 2),
          error = function(e) NULL
        )
        if (is.null(ep)) next
        polygon(ep, col = grDevices::adjustcolor(c1, alpha.f = 0.15),
          border = c1, lwd = 1.5)
      }
    }
  }
  ## lower triangle: mean silhouette for that PC pair + phenotype
  ## pairs() lower.panel(x, y) gets x = V[, j], y = V[, i] with i > j
  panel_lower <- function(x, y, ...) {
    ix <- pc_index(x)
    iy <- pc_index(y)
    s <- if (!is.na(ix) && !is.na(iy)) sil_mat[ix, iy] else NA_real_
    ## soft fill by silhouette (−1..1 → blue–grey–red)
    if (is.finite(s)) {
      t01 <- max(0, min(1, (s + 1) / 2))
      bg <- grDevices::colorRampPalette(
        c("#3B4CC0", "#DDDDDD", "#B40426"))(101)[round(t01 * 100) + 1]
      bg <- adjustcolor(bg, alpha.f=0.5)
      usr <- graphics::par("usr")
      graphics::rect(usr[1], usr[3], usr[2], usr[4], col = bg, border = NA)
    }
    xl <- mean(range(x, na.rm = TRUE))
    yl <- mean(range(y, na.rm = TRUE))
    dy <- diff(range(y, na.rm = TRUE))
    ## if (!is.na(ix) && !is.na(iy)) {
    ##   graphics::text(xl, yl + 0.18 * dy,
    ##     labels = sprintf("PC%d vs PC%d", min(ix, iy), max(ix, iy)),
    ##     cex = 0.95, col = "grey25")
    ## }
    lbl <- if (is.finite(s)) sprintf("sil = %.2f", s) else "sil = NA"
    #graphics::text(xl, yl - 0.05 * dy, labels = lbl, cex = 1.4, font = 2)
    graphics::text(xl, yl, labels = lbl, cex = 1.8, font = 1)    
  }

  graphics::pairs(
    V,
    upper.panel = panel_upper,
    lower.panel = panel_lower,
    col = cc, pch = 19, cex = 2,
    oma = c(3, 2, 6, 2)
  )
  graphics::title(
    main = paste0("PC scatterpair matrix  ·  silhouette by ", colorby_var),
    line = 2.4, cex.main = 1.25
  )
}


## Boxplots of sample PC scores stratified by a phenotype (colorby).
## Complements the pairs matrix (bivariate) and the F-test panel
## (which phenotype separates on which PC).
pcaexplorer.plot_pc_boxplots <- function(res, colorby_var, npc = 4) {
  Y <- res$Y
  V <- res$v
  if (is.null(V) || is.null(Y) || is.null(colorby_var) ||
        !colorby_var %in% colnames(Y)) {
    plot.new()
    text(0.5, 0.5, "no data")
    return(invisible())
  }
  npc <- min(as.integer(npc), ncol(V))
  group <- factor(Y[rownames(V), colorby_var])
  clrs <- grDevices::palette()[1 + seq_along(levels(group))]
  names(clrs) <- levels(group)

  ## 2x2 for 4 PCs; fall back to a single row otherwise
  if (npc <= 2) {
    nr <- 1; nc <- npc
  } else if (npc <= 4) {
    nr <- 2; nc <- 2
  } else {
    nr <- 2; nc <- ceiling(npc / 2)
  }

  bm <- 2.5 + 2.5*min((max(nchar(as.character(group)))/12),1)
  las <- 2
  if(length(levels(group)) <= 3) {
    bm <- 2
    las <- 1
  } 
  par(mfrow = c(nr, nc), mar = c(bm, 3.5, 2.2, 0.3),
    mgp = c(2.3, 0.7, 0), cex = 0.95)
  
  varx <- 100 * res$d**2 / sum(res$d**2)
  for (k in seq_len(npc)) {
    sc <- V[, k]
    boxplot(
      sc ~ group,
      col = clrs[levels(group)],
      las = las, outline = FALSE,
      xlab = "", ylab = "score",
      main = sprintf("PC%d (%.1f%%)", k, varx[k]),
      cex.main = 1.15, cex.axis = 0.9
    )
    ## jittered points on top
    set.seed(1)
    points(
      jitter(as.integer(group), amount = 0.12), sc,
      pch = 19, cex = 0.7 #, col = grDevices::adjustcolor(clrs[as.character(group)], 0.7)
    )
    abline(h = 0, lty = 3, col = "grey60")
  }
}


pcaexplorer.plot_ftest_vs_pcdim <- function(res) {
  Y <- res$Y
  V <- res$v
  sel <- which(apply(Y, 2, function(x) max(table(x[!is.na(x)])) > 1))
  F <- matrix(NA, ncol(V), length(sel))
  P <- matrix(NA, ncol(V), length(sel))
  rownames(F) = rownames(P) = colnames(V)
  colnames(F) = colnames(P) = colnames(Y)  
  for(i in sel) {
    group <- factor(Y[,i])
    ft <- matrixTests::row_oneway_welch( t(V), group )
    F[,i] <- ft$statistic
    P[,i] <- ft$pvalue    
  }

  if(0) {
    matplot(log(1+F), type='p', pch=19, cex=1.3,
      xlab = "PC dimension",
      ylab = "F-statistic    log(F+1)"
    )
    matlines(log(1+F), lty=2, lwd=1.5)   
  } else {
    matplot(-log10(P), type='p', pch=19, cex=1.3,
      xlab = "PC dimension",
      ylab = "F-test p-value    (-logP)"
    )
    matlines(-log10(P), lty=2, lwd=1.5)  
  }
  
  title("F-test vs. PC dimension", line=1.2)
  legend("topright", legend=colnames(F), fill=1:10,
    cex=0.85, y.intersp=0.85 )
  
}


## ---------------------------------------------------------------------------
## Biplot side panels
## ---------------------------------------------------------------------------

pcaexplorer.plot_variance_proportion <- function(res) {
  aa <- paste0("PC", seq_len(ncol(res$v)))
  dd <- 100 * res$d**2 / sum(res$d**2)
  graphics::par(cex = 1.1, mar = c(3, 4, 3, 0.5))
  graphics::barplot(
    dd, names.arg = aa,
    ylab = "variance proportion (%)", width = 1,
    ylim = c(0, max(dd) * 1.1)
  )
  graphics::title("Variance vs. PC dimension")
  xx <- (seq_along(dd) - 0.5) * 1.2
  graphics::lines(xx, dd, type = "b", cex = 1.6, pch = 20)
  graphics::text(xx, dd, labels = paste0(round(dd, 2), "%"), pos = 3, cex = 1)
}

pcaexplorer.plot_variance_cumulative <- function(res) {
  aa <- paste0("PC", seq_len(ncol(res$v)))
  dd <- cumsum(100 * res$d**2 / sum(res$d**2))
  graphics::par(cex = 1.1, mar = c(3, 4, 3, 0.5))
  graphics::barplot(
    dd, names.arg = aa,
    ylab = "Cumulative variance (%)", width = 1,
    ylim = c(0, max(dd) * 1.1)
  )
  graphics::title("Variance explained")
  xx <- (seq_along(dd) - 0.5) * 1.2
  graphics::lines(xx, dd, type = "b", cex = 1.6, pch = 20)
  graphics::text(xx, dd, labels = paste0(round(dd, 2), "%"), pos = 3, cex = 1)
}

## ---------------------------------------------------------------------------
## Loadings board
## ---------------------------------------------------------------------------

## Vertical stack of top-|loading| barplots (one panel per PC).
pcaexplorer.plot_loadings_bars <- function(res, npc = 4, ntop = 40) {
  npc <- min(as.integer(npc), ncol(res$u))
  graphics::par(mfrow = c(npc, 1), mar = c(6, 4.5, 1.2, 0.8), cex = 1)
  for (k in seq_len(npc)) {
    u1 <- res$u[, k]
    sel <- utils::head(order(-abs(u1)), ntop)
    v <- sort(u1[sel])
    graphics::barplot(
      v, las = 3, cex.names = 0.95, horiz = FALSE,
      ylab = "loading", xlab = ""
    )
    graphics::title(paste0("PC", k), cex.main = 1.25, line = 0.2)
    graphics::abline(h = 0, col = "grey50", lwd = 0.8)
  }
}

## Compact loading bars for PC1–PC2 (feature board side panel).
## Optional min_fc keeps only features with max |FC| >= threshold.
pcaexplorer.plot_loadings_bars_compact <- function(res, npc = 2, ntop = 40,
                                                   normX = NULL, Y = NULL,
                                                   min_fc = 0) {
  npc <- min(as.integer(npc), ncol(res$u))
  U <- res$u
  if (is.null(normX)) normX <- res$X
  if (is.null(Y)) Y <- res$Y
  min_fc <- suppressWarnings(as.numeric(min_fc))
  if (is.finite(min_fc) && min_fc > 0 && !is.null(normX)) {
    H <- if (!is.null(res$H)) res$H else playbase::expandPhenoMatrix(Y, drop.ref = FALSE)
    genes <- intersect(rownames(U), rownames(normX))
    max_fc <- pcaexplorer.feature_max_fc(normX[genes, , drop = FALSE], H)
    keep <- pcaexplorer.filter_features_by_maxfc(genes, max_fc, min_fc)
    if (length(keep) >= 2L) {
      U <- U[keep, , drop = FALSE]
    }
  }
  graphics::par(mfrow = c(npc, 1), mar = c(6, 4.5, 1.4, 0.2), cex = 1)
  for (k in seq_len(npc)) {
    u1 <- U[, k]
    sel <- utils::head(order(-abs(u1)), ntop)
    graphics::barplot(
      sort(u1[sel]), las = 3, cex.names = 0.95, horiz = FALSE,
      xlab = "", ylab = "loading"
    )
    graphics::title(paste0("PC", k), cex.main = 1.1, line = 0.2)
  }
}

## Expression heatmap of top genes per PC (by |loading|), split into PC
## blocks via gx.splitmap; phenotype annotations on top of columns.
pcaexplorer.plot_loadings_heatmap <- function(res, npc = 4, ntop = 15) {
  X <- res$X
  Y <- res$Y
  U <- res$u
  npc <- min(as.integer(npc), ncol(U))
  pc_labs <- paste0("PC", seq_len(npc))

  blocks <- lapply(seq_len(npc), function(k) {
    u1 <- U[, k]
    g <- utils::head(names(sort(abs(u1), decreasing = TRUE)), ntop)
    g <- g[g %in% rownames(X)]
    if (!length(g)) {
      return(NULL)
    }
    list(genes = g, pc = pc_labs[k])
  })
  blocks <- Filter(Negate(is.null), blocks)
  if (!length(blocks)) {
    plot.new()
    text(0.5, 0.5, "no genes selected")
    return(invisible())
  }

  gene_vec <- unlist(lapply(blocks, `[[`, "genes"), use.names = FALSE)
  split_vec <- unlist(
    lapply(blocks, function(b) rep(b$pc, length(b$genes))),
    use.names = FALSE
  )
  split_vec <- factor(split_vec, levels = pc_labs)

  M <- X[gene_vec, , drop = FALSE]
  ## unique rownames for ComplexHeatmap; allow gene in multiple PCs
  rownames(M) <- make.unique(gene_vec, sep = " ")
  names(split_vec) <- rownames(M)

  annot <- as.data.frame(Y, stringsAsFactors = FALSE)
  nlev <- vapply(annot, function(z) length(unique(z[!is.na(z)])), integer(1))
  keep <- nlev >= 2 & nlev <= min(20L, max(2L, ncol(M) - 1L))
  if (any(keep)) {
    annot <- annot[, keep, drop = FALSE]
  } else {
    annot <- NULL
  }

  graphics::par(cex = 1.4)
  playbase::gx.splitmap(
    M,
    split = split_vec,
    splitx = NULL,
    col.annot = annot,
    scale = "row",
    nmax = nrow(M),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    order.groups = pc_labs,
    rowlab.maxlen = 30,
    rownames_width = 30 + 30,
    show_colnames = ncol(M) <= 40,
    show_rownames = 99,
    show_legend = TRUE,
    show_key = TRUE,
    annot.ht = 4,
    annot.cex = 1,
    cexRow = 1,
    cexCol = 1,
    mar = c(4, 2, 2, 8),
    main = paste0("Top genes by |loading| (PC1–PC", npc, ")")
  )
}


## ---------------------------------------------------------------------------
## Feature PCA board
## ---------------------------------------------------------------------------

## Per-feature max |fold-change| vs expanded phenotype (same metric used to
## color feature PCA points). Returns a named numeric vector.
pcaexplorer.feature_max_fc <- function(normX, H) {
  if (is.null(normX) || is.null(H) || !nrow(normX) || !ncol(H)) {
    return(numeric(0))
  }
  rx <- fast_diff_in_means(normX, t(H))
  apply(abs(rx), 1, function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) return(NA_real_)
    max(z, na.rm = TRUE)
  })
}

## Keep features with maxFC >= min_fc (drops the middle / weak-effect cloud).
pcaexplorer.filter_features_by_maxfc <- function(genes, max_fc, min_fc = 0) {
  min_fc <- suppressWarnings(as.numeric(min_fc))
  if (!is.finite(min_fc) || min_fc <= 0) {
    return(genes)
  }
  mf <- max_fc[genes]
  keep <- genes[is.finite(mf) & mf >= min_fc]
  keep
}

## Condition / phenotype direction vectors on the feature PC plane (2D or 3D).
## Same recipe as feature_pca_all: colMeans(pos * pmax(cor,0)^2), then
## unit-normalized (length is applied later via var.scale * score_span).
pcaexplorer.feature_condition_dirs <- function(pos, rx, numarrows = Inf) {
  if (is.null(pos) || is.null(rx) || !nrow(pos) || !ncol(rx)) {
    return(NULL)
  }
  genes <- intersect(rownames(pos), rownames(rx))
  if (length(genes) < 2L) {
    return(NULL)
  }
  pos <- as.matrix(pos[genes, , drop = FALSE])
  rx <- as.matrix(rx[genes, , drop = FALSE])
  nd <- ncol(pos)

  D <- matrix(NA_real_, ncol(rx), nd,
    dimnames = list(colnames(rx), colnames(pos)))
  for (i in seq_len(ncol(rx))) {
    w <- pmax(rx[, i], 0, na.rm = TRUE)**2
    if (!any(w > 0, na.rm = TRUE)) next
    ## match feature_pca_all: colMeans(pos * w), not a proper weighted mean
    dir <- colMeans(pos * w, na.rm = TRUE)
    nrm <- sqrt(sum(dir**2))
    if (!is.finite(nrm) || nrm <= 0) next
    D[i, ] <- dir / nrm
  }
  D <- D[stats::complete.cases(D), , drop = FALSE]
  if (!nrow(D)) {
    return(NULL)
  }
  nkeep <- suppressWarnings(as.integer(numarrows))
  if (is.finite(nkeep) && nkeep >= 0L && nkeep < nrow(D)) {
    if (nkeep == 0L) {
      return(D[FALSE, , drop = FALSE])
    }
    r2 <- rowSums(D^2)
    D <- D[utils::head(order(-r2), nkeep), , drop = FALSE]
  }
  D
}

## Interactive feature PCA (plotly): genes as points on PC1–PC2 with biplot
## arrows for phenotype / condition directions (same as feature_pca_all).
## Labels the top-n up/down features by loading on PC1 and on PC2.
pcaexplorer.plot_feature_pca <- function(res, normX = NULL, Y = NULL,
                                         numarrows = 20,
                                         nlabel = 5,
                                         min_fc = 0,
                                         marker.size = 6,
                                         var.scale = 0.75,
                                         main = "Feature PCA (interactive)") {
  if (is.null(normX)) normX <- res$X
  if (is.null(Y)) Y <- res$Y

  pos <- res$u[, 1:2, drop = FALSE]
  colnames(pos) <- c("PC1", "PC2")
  H <- if (!is.null(res$H)) {
    res$H
  } else {
    playbase::expandPhenoMatrix(Y, drop.ref = FALSE)
  }
  genes <- intersect(rownames(pos), rownames(normX))
  if (length(genes) < 2L) {
    return(plotly::plotly_empty(type = "scatter") %>%
      plotly::layout(title = "no features"))
  }
  pos <- pos[genes, , drop = FALSE]
  rx <- fast_diff_in_means(normX[genes, , drop = FALSE], t(H))  ## fold.change

  ## color genes by max |FC| across conditions; drop weak middle cloud
  max_fc <- apply(abs(rx), 1, function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) return(NA_real_)
    max(z, na.rm = TRUE)
  })
  keep <- pcaexplorer.filter_features_by_maxfc(rownames(pos), max_fc, min_fc)
  if (length(keep) < 2L) {
    return(plotly::plotly_empty(type = "scatter") %>%
      plotly::layout(title = paste0("no features with maxFC ≥ ", min_fc)))
  }
  pos <- pos[keep, , drop = FALSE]
  rx <- rx[keep, , drop = FALSE]
  max_fc <- max_fc[keep]

  ## top correlated condition name for hover
  top_cond <- apply(rx, 1, function(z) {
    if (!any(is.finite(z))) return("")
    colnames(rx)[which.max(replace(z, !is.finite(z), -Inf))]
  })

  df <- data.frame(
    x = pos[, 1], y = pos[, 2],
    gene = rownames(pos),
    max_fc = max_fc,
    top = top_cond,
    hover = paste0(
      "<b>", rownames(pos), "</b><br>",
      "PC1: ", sprintf("%.3f", pos[, 1]), "<br>",
      "PC2: ", sprintf("%.3f", pos[, 2]), "<br>",
      "top trait: ", top_cond, "<br>",
      "max FC: ", sprintf("%.2f", max_fc)
    ),
    stringsAsFactors = FALSE
  )

  ## top-n up / down loadings on each of PC1 and PC2
  nlab <- max(0L, as.integer(nlabel))
  label_genes <- character(0)
  if (nlab > 0L) {
    for (k in 1:2) {
      v <- pos[, k]
      v <- v[is.finite(v)]
      if (!length(v)) next
      up <- names(utils::head(sort(v, decreasing = TRUE), nlab))
      dn <- names(utils::head(sort(v, decreasing = FALSE), nlab))
      label_genes <- c(label_genes, up, dn)
    }
    label_genes <- unique(label_genes)
    label_genes <- label_genes[label_genes %in% rownames(pos)]
  }

  fig <- plotly::plot_ly(df) %>%
    plotly::add_markers(
      x = ~x, y = ~y,
      type = "scatter", mode = "markers",
      color = ~max_fc,
      colors = c("#3B4CC0", "#DDDDDD", "#B40426"),
      marker = list(
        size = marker.size,
        line = list(color = "rgba(40,40,40,0.25)", width = 0.4),
        colorbar = list(title = "max FC", len = 0.4, y = 0.8)
      ),
      text = ~hover, hoverinfo = "text",
      name = "features",
      showlegend = FALSE
    )

  ## gene labels for extreme PC1/PC2 loadings
  if (length(label_genes)) {
    dfl <- df[df$gene %in% label_genes, , drop = FALSE]
    ## textposition by quadrant so labels sit outside the cloud a bit
    tpos <- ifelse(dfl$x >= 0 & dfl$y >= 0, "top right",
            ifelse(dfl$x < 0 & dfl$y >= 0, "top left",
            ifelse(dfl$x < 0 & dfl$y < 0, "bottom left", "bottom right")))
    fig <- fig %>% plotly::add_text(
      data = dfl,
      x = ~x, y = ~y, text = ~gene,
      type = "scatter", mode = "text",
      textposition = tpos,
      textfont = list(color = "grey15", size = 11),
      hoverinfo = "skip",
      showlegend = FALSE,
      inherit = FALSE
    )
  }

  ## condition-direction arrows (biplot overlay)
  D <- pcaexplorer.feature_condition_dirs(pos, rx, numarrows = numarrows)
  arrow_ann <- list()
  if (!is.null(D) && nrow(D) && max(abs(D), na.rm = TRUE) > 0) {
    score_span <- max(abs(pos), na.rm = TRUE)
    s <- var.scale * score_span / max(abs(D), na.rm = TRUE)
    ax <- D[, 1] * s
    ay <- D[, 2] * s
    r <- sqrt(ax^2 + ay^2)
    r[r == 0] <- 1
    nudge <- 0.05 * score_span
    lx <- ax + nudge * ax / r
    ly <- ay + nudge * ay / r

    nD <- nrow(D)
    arrow_ann <- c(
      lapply(seq_len(nD), function(i) {
        list(
          x = ax[i], y = ay[i], ax = 0, ay = 0,
          xref = "x", yref = "y", axref = "x", ayref = "y",
          showarrow = TRUE, text = "",
          arrowhead = 3, arrowsize = 1.25, arrowwidth = 2.4,
          arrowcolor = "black"
        )
      }),
      lapply(seq_len(nD), function(i) {
        list(
          x = lx[i], y = ly[i],
          xref = "x", yref = "y",
          showarrow = FALSE,
          text = rownames(D)[i],
          xanchor = if (ax[i] >= 0) "left" else "right",
          yanchor = if (ay[i] >= 0) "bottom" else "top",
          font = list(color = "black", size = 15, face = 'bold'),
          bgcolor = "rgba(255,255,255,0.7)",
          borderpad = 2
        )
      })
    )
  }

  lim <- c(-1, 1) * max(abs(pos), na.rm = TRUE) * 1.15
  ## simple_white-like layout (named templates not available on plotly 4.12)
  axis_sw <- function(title_txt, ...) {
    utils::modifyList(
      list(
        title = list(text = title_txt, standoff = 8),
        showgrid = TRUE, gridcolor = "rgb(230,230,230)", gridwidth = 1,
        zeroline = TRUE, zerolinecolor = "rgb(200,200,200)", zerolinewidth = 1,
        showline = TRUE, linecolor = "black", linewidth = 1,
        mirror = TRUE, ticks = "outside", ticklen = 4, automargin = TRUE
      ),
      list(...)
    )
  }

  fig %>% plotly::layout(
    title = list(text = main, font = list(size = 16)),
    xaxis = axis_sw("feature PC1", range = lim, scaleanchor = "y", scaleratio = 1),
    yaxis = axis_sw("feature PC2", range = lim),
    annotations = arrow_ann,
    paper_bgcolor = "white",
    plot_bgcolor = "white",
    margin = list(l = 60, r = 20, t = 50, b = 50)
  )
}

## Interactive 3D feature PCA (plotly scatter3d): genes on PC1–PC3 with
## condition-direction arrows and top up/down feature labels per PC.
pcaexplorer.plot_feature_pca_3d <- function(res, normX = NULL, Y = NULL,
                                            numarrows = 20,
                                            nlabel = 5,
                                            min_fc = 0,
                                            marker.size = 3,
                                            var.scale = 0.75,
                                            main = "Feature PCA 3D") {
  if (is.null(normX)) normX <- res$X
  if (is.null(Y)) Y <- res$Y

  npc <- min(3L, ncol(res$u))
  if (npc < 3L) {
    return(plotly::plotly_empty(type = "scatter3d") %>%
      plotly::layout(title = "need at least 3 feature PCs"))
  }

  pos <- res$u[, 1:3, drop = FALSE]
  colnames(pos) <- c("PC1", "PC2", "PC3")
  H <- if (!is.null(res$H)) {
    res$H
  } else {
    playbase::expandPhenoMatrix(Y, drop.ref = FALSE)
  }
  genes <- intersect(rownames(pos), rownames(normX))
  if (length(genes) < 2L) {
    return(plotly::plotly_empty(type = "scatter3d") %>%
      plotly::layout(title = "no features"))
  }
  pos <- pos[genes, , drop = FALSE]
  rx <- fast_diff_in_means(normX[genes, , drop = FALSE], t(H))

  max_fc <- apply(abs(rx), 1, function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) return(NA_real_)
    max(z, na.rm = TRUE)
  })
  keep <- pcaexplorer.filter_features_by_maxfc(rownames(pos), max_fc, min_fc)
  if (length(keep) < 2L) {
    return(plotly::plotly_empty(type = "scatter3d") %>%
      plotly::layout(title = paste0("no features with maxFC ≥ ", min_fc)))
  }
  pos <- pos[keep, , drop = FALSE]
  rx <- rx[keep, , drop = FALSE]
  max_fc <- max_fc[keep]

  top_cond <- apply(rx, 1, function(z) {
    if (!any(is.finite(z))) return("")
    colnames(rx)[which.max(replace(z, !is.finite(z), -Inf))]
  })

  df <- data.frame(
    x = pos[, 1], y = pos[, 2], z = pos[, 3],
    gene = rownames(pos),
    max_fc = max_fc,
    top = top_cond,
    hover = paste0(
      "<b>", rownames(pos), "</b><br>",
      "PC1: ", sprintf("%.3f", pos[, 1]), "<br>",
      "PC2: ", sprintf("%.3f", pos[, 2]), "<br>",
      "PC3: ", sprintf("%.3f", pos[, 3]), "<br>",
      "top trait: ", top_cond, "<br>",
      "max FC: ", sprintf("%.2f", max_fc)
    ),
    stringsAsFactors = FALSE
  )

  ## top-n up / down loadings on each of PC1–PC3
  nlab <- max(0L, as.integer(nlabel))
  label_genes <- character(0)
  if (nlab > 0L) {
    for (k in 1:3) {
      v <- pos[, k]
      v <- v[is.finite(v)]
      if (!length(v)) next
      up <- names(utils::head(sort(v, decreasing = TRUE), nlab))
      dn <- names(utils::head(sort(v, decreasing = FALSE), nlab))
      label_genes <- c(label_genes, up, dn)
    }
    label_genes <- unique(label_genes)
    label_genes <- label_genes[label_genes %in% rownames(pos)]
  }

  fig <- plotly::plot_ly(df) %>%
    plotly::add_markers(
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d", mode = "markers",
      color = ~max_fc,
      colors = c("#3B4CC0", "#DDDDDD", "#B40426"),
      marker = list(
        size = marker.size,
        line = list(color = "rgba(40,40,40,0.25)", width = 0.3),
        colorbar = list(title = "max FC", len = 0.4, y = 0.8)
      ),
      text = ~hover, hoverinfo = "text",
      name = "features",
      showlegend = FALSE
    )

  if (length(label_genes)) {
    dfl <- df[df$gene %in% label_genes, , drop = FALSE]
    fig <- fig %>% plotly::add_text(
      data = dfl,
      x = ~x, y = ~y, z = ~z, text = ~gene,
      type = "scatter3d", mode = "text",
      textfont = list(color = "grey15", size = 11),
      hoverinfo = "skip",
      showlegend = FALSE,
      inherit = FALSE
    )
  }

  ## condition-direction arrows as scatter3d lines + tip text
  D <- pcaexplorer.feature_condition_dirs(pos, rx, numarrows = numarrows)
  if (!is.null(D) && nrow(D) && max(abs(D), na.rm = TRUE) > 0) {
    score_span <- max(abs(pos), na.rm = TRUE)
    s <- var.scale * score_span / max(abs(D), na.rm = TRUE)
    ax <- D[, 1] * s
    ay <- D[, 2] * s
    az <- D[, 3] * s
    r <- sqrt(ax^2 + ay^2 + az^2)
    r[r == 0] <- 1
    nudge <- 0.06 * score_span
    lx <- ax + nudge * ax / r
    ly <- ay + nudge * ay / r
    lz <- az + nudge * az / r

    for (i in seq_len(nrow(D))) {
      fig <- fig %>% plotly::add_trace(
        x = c(0, ax[i]), y = c(0, ay[i]), z = c(0, az[i]),
        type = "scatter3d", mode = "lines",
        line = list(color = "black", width = 6),
        showlegend = FALSE, hoverinfo = "skip", inherit = FALSE
      )
    }
    fig <- fig %>% plotly::add_text(
      x = lx, y = ly, z = lz, text = rownames(D),
      type = "scatter3d", mode = "text",
      textfont = list(color = "black", size = 14),
      showlegend = FALSE, hoverinfo = "text", inherit = FALSE
    )
  }

  fig %>% plotly::layout(
    title = list(text = main, font = list(size = 16)),
    scene = list(
      xaxis = list(title = "feature PC1"),
      yaxis = list(title = "feature PC2"),
      zaxis = list(title = "feature PC3"),
      aspectmode = "cube"
    ),
    paper_bgcolor = "white",
    margin = list(l = 0, r = 0, t = 50, b = 0)
  )
}

#' Handles NA
#' 
fast_diff_in_means <- function(X, Y) {
  # X: vars x samples (contains NAs)
  # Y: groups x samples (binary, no NAs assumed)
  
  # 1. Create a binary mask of where data actually exists
  valid_mask <- !is.na(X)
  
  # 2. Replace NAs with 0 strictly for the matrix multiplication step
  X_zeroed <- X
  X_zeroed[!valid_mask] <- 0
  
  # 3. Compute sums for group 1 and group 0
  sum_1 <- X_zeroed %*% t(Y)
  
  # Total sum of valid elements per variable
  sum_total <- matrix(rowSums(X_zeroed), nrow = nrow(X), ncol = nrow(Y))
  sum_0 <- sum_total - sum_1
  
  # 4. CRITICAL: Dynamically calculate sample sizes using the valid mask
  # This counts how many non-NA samples fall into group 1 for each variable
  n1 <- valid_mask %*% t(Y)
  
  # Count how many non-NA samples fall into group 0
  n_total <- matrix(rowSums(valid_mask), nrow = nrow(X), ncol = nrow(Y))
  n0 <- n_total - n1
  
  # 5. Compute means safely, generating NA if a group has 0 valid samples
  mean_1 <- sum_1 / n1
  mean_0 <- sum_0 / n0
  
  # Prevent 0/0 resulting in NaN instead of NA
  mean_1[n1 == 0] <- NA
  mean_0[n0 == 0] <- NA
  
  return(mean_1 - mean_0)
}

## Feature scores (res$u) colored by FC with each expanded phenotype.
pcaexplorer.plot_feature_pca_all <- function(res, normX = NULL, Y = NULL,
                                             min_fc = 0,
                                             var.scale = 0.65,
                                             show_arrows = TRUE) {
  if (is.null(normX)) normX <- res$X
  if (is.null(Y)) Y <- res$Y
  pos <- res$u[, 1:2, drop = FALSE]
  genes <- intersect(rownames(pos), rownames(normX))
  pos <- pos[genes, , drop = FALSE]
  normX <- normX[genes, , drop = FALSE]

  H <- if (!is.null(res$H)) res$H else playbase::expandPhenoMatrix(Y, drop.ref = FALSE)
  rx <- fast_diff_in_means(normX, t(H))
  max_fc <- apply(abs(rx), 1, function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) return(NA_real_)
    max(z, na.rm = TRUE)
  })
  keep <- pcaexplorer.filter_features_by_maxfc(rownames(pos), max_fc, min_fc)
  if (length(keep) < 2L) {
    plot.new()
    text(0.5, 0.5, paste0("no features with maxFC ≥ ", min_fc))
    return(invisible())
  }
  pos <- pos[keep, , drop = FALSE]
  rx <- rx[keep, , drop = FALSE]

  vs <- suppressWarnings(as.numeric(var.scale))
  if (!is.finite(vs) || vs <= 0) vs <- 0.65
  ## map slider (typical 0.1–1.5) onto arrow length in feature-PC units
  arrow_len <- 0.08 * vs

  graphics::par(
    mfrow = c(3, 3), mar = c(3.4, 4, 2.4, 1),
    mgp = c(2.4, 0.8, 0), cex = 1
  )
  if (ncol(rx) > 9) {
    graphics::par(mfrow = c(4, 4))
  }
  for (i in seq_len(min(ncol(rx), 16L))) {
    rr <- range(rx, na.rm = TRUE)  ## fix color scale range
    cc <- playbase::colorscale(c(rr, rx[, i]))
    cc <- cc[-c(1, 2)]
    graphics::plot(pos, col = cc, xlab = "PC1", ylab = "PC2", pch = 19)
    graphics::title(main = colnames(rx)[i], line = 0.8, cex.main = 1.4)

    if (isTRUE(show_arrows)) {
      dir <- colMeans(pos * pmax(rx[, i], 0, na.rm = TRUE)**2)
      nrm <- sqrt(sum(dir**2))
      if (is.finite(nrm) && nrm > 0) {
        dir <- arrow_len * dir / nrm
        graphics::arrows(0, 0, dir[1], dir[2], length = 0.15, lwd = 2.8)
      }
    }
  }
}

## Heatmap of PC × expanded-phenotype correlations (res$R).
pcaexplorer.plot_trait_correlation <- function(res, min_fc = 0) {
  R <- res$R
  if (is.null(R)) {
    plot.new()
    text(0.5, 0.5, "no trait correlations")
    return(invisible())
  }
  ## When min_fc > 0, recompute trait correlations from feature PCs of the
  ## filtered gene set is not meaningful for sample scores (res$v). Keep
  ## sample-level R; gene FC filter applies to feature scatter/loading plots.
  graphics::par(cex = 1, mar = c(4, 0.2, 3, 10), mfrow = c(1, 1))
  bluered <- gplots::colorpanel(255, low = "blue3", mid = "grey90", high = "red3")
  playbase::gx.imagemap(
    R, clust = FALSE, col = bluered, cex = 1.2,
    main = "PC-trait correlation", cex.main = 1.2
  )
}


## ---------------------------------------------------------------------------
## SD Filtering board
## ---------------------------------------------------------------------------

## Grid of PC1–PC2 scatters at successive top-SD feature cutoffs.
pcaexplorer.plot_pca_vs_topsd <- function(res, colorby_var = NULL) {
  X <- res$X
  Y <- res$Y

  nn <- c(nrow(X), 20 * 2**seq(0, 12, 2))
  nn <- utils::head(nn[nn <= nrow(X)], 9)
  nn <- sort(nn)
  sdx <- matrixStats::rowSds(X, na.rm = TRUE)
  X <- X[order(-sdx), , drop = FALSE]

  ycol <- if (!is.null(colorby_var) && colorby_var %in% colnames(Y)) {
    Y[, colorby_var]
  } else if (ncol(Y) >= 1) {
    Y[, min(2L, ncol(Y))]
  } else {
    rep(1, ncol(X))
  }
  cc <- 1 + as.integer(factor(ycol))

  graphics::par(
    mfrow = c(3, 3), mar = c(4, 4, 4, 0.5),
    cex = 0.95, mgp = c(2.3, 0.85, 0)
  )
  for (n in nn) {
    X1 <- utils::head(X, n)
    pr <- stats::prcomp(t(X1))
    graphics::plot(pr$x[, 1:2], cex = 1.7, col = cc, pch = 19)
    graphics::title(paste("N=", n))
  }
}

## Cumulative variance explained as a function of top-SD feature count,
## with the current topsd cutoff highlighted.
pcaexplorer.plot_variance_vs_topsd <- function(res, topsd = NULL) {
  X <- res$X
  nn <- c(20 * 2**seq(0, 12, 2), nrow(X))
  nn <- nn[nn <= nrow(X)]

  U <- res$u %*% diag(res$d, length(res$d))
  U <- U[order(-rowMeans(U**2)), , drop = FALSE]
  pcx <- 100 * cumsum(rowSums(U**2)) / sum(U**2)
  graphics::matplot(
    pcx, pch = ".", lwd = 3, type = "l",
    xlab = "Number of topSD features",
    ylab = "Variance explained (%)"
  )
  graphics::points(nn, pcx[nn], pch = 19, cex = 1.3)
  graphics::title("Variance vs. topSD")

  if (!is.null(topsd) && is.finite(as.numeric(topsd))) {
    topx <- as.integer(topsd)
    topx <- max(1L, min(topx, length(pcx)))
    topy <- pcx[topx]
    graphics::abline(v = topx, lty = 2, col = "red")
    graphics::abline(h = topy, lty = 2, col = "red")
  }
}

## Histogram of feature SDs with the current topsd threshold marked.
pcaexplorer.plot_sd_histogram <- function(res, topsd = NULL) {
  sdx <- matrixStats::rowSds(res$X, na.rm = TRUE)
  sdx[is.na(sdx)] <- 0

  graphics::par(mar = c(4, 4, 3, 1))
  graphics::hist(
    sdx, breaks = 100,
    main = "Histogram of SD",
    xlab = "standard deviation  (SD)"
  )
  if (!is.null(topsd) && is.finite(as.numeric(topsd))) {
    ntop <- as.integer(topsd)
    ntop <- max(1L, min(ntop, length(sdx)))
    sdthreshold <- sort(sdx, decreasing = TRUE)[ntop]
    graphics::abline(v = sdthreshold, col = "red", lty = 2, lwd = 1)
  }
}
