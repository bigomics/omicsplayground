##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Copied from board.clustering/R/clustering_plot_clusterpca.R as it stands
## after PR #1838 (biplot arrows and phenotype directions), deliberately close
## to it so the two stay diffable. What differs, and why:
##
##   - the ui/server names, because 00SourceAll.R sources every component into
##     one environment and the originals would collide;
##   - translate = FALSE, as every module in this app sets;
##   - the component selectors and the variance-explained view are dropped.
##     This app plots PC1 vs PC2 and reports variance on the axes; the
##     PC-vs-covariate heatmap beside it answers what those components are made
##     of, which is what the scree plot is usually reached for;
##   - the arrow controls are unconditional. Upstream gates them on the board's
##     layout selector, which lives in a parent namespace this app does not
##     have; arrows are dropped per panel in method_pos() instead;
##   - feature arrows are labelled "cg####### (SYMBOL)" rather than by symbol
##     alone: a CpG is the unit that was measured, and two probes in one gene
##     routinely disagree;
##   - components are computed on M-values, and the SD ranking that picks the
##     probes streams over the matrix rather than copying it.
##
## t-SNE and UMAP still come from pgx$cluster$pos, computed at upload. PCA does
## not: it is recomputed here, which is what makes loadings - and therefore
## arrows - available at all.

## Components, loadings and phenotype directions for the methylome PCA.
##
## ONE decomposition per dataset, shared by the scatter and by the
## PC-vs-covariate heatmap beside it, so the two panels report the same PC1.
##
## By default it decomposes the 1,000 most variable probes - what playbase
## (pgx.clusterSamples) and the field (minfi::mdsPlot) both do. The "Use all
## features" setting drops that cut, and the reason it exists is measured: on
## the demo cohort a top-1,000 cut leaves cell type (0.97 vs 0.96) and sex
## (0.20 vs 0.20) untouched but drops the slide/batch association from
## p 8.4e-09 to p 0.09. Batch is not a high-variance effect - it is a small
## shift over very many probes - so a variance cut hides exactly what the
## heatmap is for. Conventional and cheap by default; exhaustive when someone
## is hunting a batch effect.
##
## Cost is bounded so that holds on a real array. A 937k x 500 EPICv2 cohort is
## a 3.7 GB beta matrix and the naive route wants several copies of it, so:
##
##   - two streamed passes over probe blocks sized by cells, never the whole
##     matrix: one accumulating the n x n cross-product, one turning the
##     eigenvectors into per-probe loadings. Both are exact - addition is how a
##     cross-product decomposes over rows, and a loading is a per-probe
##     projection.
##   - above max_probes an even stride is decomposed. A stride, not a top-SD
##     cut, for the reason above; over manifest order it is stratified by
##     chromosome for free and needs no RNG.
mp_pca_components <- function(pgx, npc = 8L, reduce.sd = 1000L, max_probes = 1e5) {
  mp_require_methylomics(pgx)
  beta <- mp_beta(pgx)
  n <- ncol(beta)
  shiny::validate(shiny::need(n >= 3,
    "Too few samples for a principal component decomposition."))

  idx <- seq_len(nrow(beta))
  thinned <- FALSE
  if (isTRUE(reduce.sd > 0) && length(idx) > reduce.sd) {
    ## Rank by SD on the M scale, streamed - beta is heteroscedastic at both
    ## ends, so an SD ranking over beta selects probes sitting near 0.5 rather
    ## than probes that actually move.
    blk <- max(2000L, as.integer(2e6 / n))
    sdx <- rep(NA_real_, length(idx))
    for (from in seq(1, length(idx), by = blk)) {
      j <- from:min(from + blk - 1L, length(idx))
      m <- playbase.epigenetics::betaToM(beta[idx[j], , drop = FALSE])
      m[!is.finite(m)] <- NA
      sdx[j] <- matrixStats::rowSds(m, na.rm = TRUE)
    }
    ok <- which(is.finite(sdx) & sdx > 0)
    shiny::validate(shiny::need(length(ok) >= 3,
      "No probe varies across these samples, so there is nothing to decompose."))
    idx <- sort(idx[ok[utils::head(order(-sdx[ok]), reduce.sd)]])
  } else if (length(idx) > max_probes) {
    ## No SD cut asked for, so the only bound left is size: an even stride,
    ## stratified by chromosome for free and needing no RNG.
    thinned <- TRUE
    idx <- unique(round(seq(1, length(idx), length.out = max_probes)))
  }
  ## Small block budget: betaToM costs about 8x the block handed to it, so the
  ## block - not the array - is what that multiplier applies to.
  block <- max(2000L, as.integer(2e6 / n))
  blocks <- split(idx, ceiling(seq_along(idx) / block))

  ## Pass one: the cross-product, and the complete-probe set each block keeps.
  ## Complete probes only - an eigendecomposition tolerates no gaps, and
  ## imputing would invent the structure this is meant to detect.
  G <- matrix(0, n, n)
  keep <- vector("list", length(blocks))
  for (i in seq_along(blocks)) {
    b <- beta[blocks[[i]], , drop = FALSE]
    ok <- rowSums(is.na(b)) == 0
    keep[[i]] <- blocks[[i]][ok]
    if (!any(ok)) next
    m <- playbase.epigenetics::betaToM(b[ok, , drop = FALSE])
    G <- G + crossprod(m - rowMeans(m))
  }
  used <- sum(lengths(keep))
  shiny::validate(shiny::need(used > 100,
    "Almost every probe carries a missing value, so no complete matrix is left to decompose."))

  npc <- max(1L, min(as.integer(npc), n - 1L))
  e <- eigen(G, symmetric = TRUE)
  ev <- pmax(e$values, 0)
  pos <- e$vectors[, seq_len(npc), drop = FALSE]
  d <- sqrt(ev[seq_len(npc)])
  d[d == 0] <- 1

  ## Pass two: loadings. G = X'X, so the eigenvectors are V and the per-probe
  ## loadings are X V / d - computable a block at a time.
  loadings <- matrix(NA_real_, used, npc)
  rn <- character(used)
  at <- 0L
  for (i in seq_along(blocks)) {
    if (!length(keep[[i]])) next
    m <- playbase.epigenetics::betaToM(beta[keep[[i]], , drop = FALSE])
    m <- m - rowMeans(m)
    r <- seq.int(at + 1L, at + nrow(m))
    loadings[r, ] <- sweep(m %*% pos, 2, d, "/")
    rn[r] <- rownames(m)
    at <- at + nrow(m)
  }

  ## An SVD's sign is arbitrary. Pin the largest-magnitude loading positive so
  ## the plot does not mirror between sessions, and flip scores with it or the
  ## arrows point the wrong way.
  for (i in seq_len(npc)) {
    if (loadings[which.max(abs(loadings[, i])), i] < 0) {
      loadings[, i] <- -loadings[, i]
      pos[, i] <- -pos[, i]
    }
  }
  ## Scores, not directions. The eigenvectors are unit-norm, so returning them
  ## unscaled draws every component at the same spread whatever its variance -
  ## PC8 as wide as PC1, under axis labels stating otherwise. X V = U D, and
  ## G's eigenvalues are D^2, so d is the scale. This happens after the
  ## loadings pass because that one needs the unscaled V.
  pos <- sweep(pos, 2, d, "*")
  dimnames(pos) <- list(colnames(beta), paste0("PC", seq_len(npc)))
  dimnames(loadings) <- list(rn, paste0("PC", seq_len(npc)))

  ## Phenotype directions: each level one-hot encoded and correlated with each
  ## component. Same variables the heatmap tests, which is now the same object.
  pheno.cor <- NULL
  vars <- intersect(mp_model_vars(pgx), colnames(pgx$samples))
  if (length(vars)) {
    pheno.cor <- tryCatch({
      H <- playbase::expandPhenoMatrix(
        pgx$samples[rownames(pos), vars, drop = FALSE], drop.ref = FALSE)
      rho <- suppressWarnings(stats::cor(H, pos, use = "pairwise.complete.obs"))
      rho[stats::complete.cases(rho), , drop = FALSE]
    }, error = function(e) NULL)
  }

  list(pos = pos, loadings = loadings, pheno.cor = pheno.cor,
       varexp = round(100 * ev[seq_len(npc)] / sum(ev), 1),
       n_probes = used, thinned = thinned)
}

methylome_plot_clustpca_ui <- function(
  id,
  label = "",
  height,
  width,
  title,
  info.text,
  info.methods,
  info.references,
  info.extra_link,
  caption,
  parent
) {
  ns <- shiny::NS(id)

  ## placement = "left" on every tooltip: withTooltip() defaults to "bottom",
  ## and inside the options dropdown that opens each tip downwards, over the
  ## control below it. The dropdown hangs off the panel's top-right corner, so
  ## the open space is to its left.
  plot_opts <- shiny::tagList(
    ## Upstream reads the layout from the Clustering board's sidebar, in the
    ## parent namespace. This app has no such sidebar, so the selector lives
    ## with the panel it drives - which is also where a reader looks for it.
    withTooltip(
      shiny::selectInput(ns("clustmethod"), "Layout:",
        choices = c("PCA" = "pca", "t-SNE" = "tsne", "UMAP" = "umap"),
        selected = "pca", width = "100%"
      ),
      "PCA keeps distances and axis percentages quantitative. t-SNE and UMAP separate groups more readably, but distances between clusters carry no meaning and neither has ordered components, so the biplot arrows are dropped there.",
      placement = "left"
    ),
    withTooltip(
      shiny::selectInput(ns("hmpca.colvar"), "Color/label:", choices = NULL, width = "100%"),
      "Set colors/labels according to a given phenotype.",
      placement = "left"
    ),
    withTooltip(
      shiny::selectInput(ns("hmpca.shapevar"), "Shape:", choices = NULL, width = "100%"),
      "Set shapes according to a given phenotype.",
      placement = "left"
    ),
    withTooltip(
      shiny::selectInput(
        ns("pca_label"),
        label = "Label:",
        choices = list("group", "bottom", "sample", "<none>")
      ),
      "Place group labels as legend at the bottom or in plot as group or sample labels.",
      placement = "left"
    ),
    ## Upstream gates these on the board's layout selector, in the PARENT
    ## namespace. This app has no such input - the layout is fixed to PCA and
    ## the other embeddings are reached through "show all methods" - so the
    ## controls stand unconditional. Arrows are dropped for any non-PCA panel
    ## in method_pos(), which is where that decision actually belongs.
    shiny::conditionalPanel(
      "input.clustmethod == 'pca'",
      ns = ns,
      withTooltip(
        shiny::selectInput(ns("pca_arrows"), "Arrows:",
          choices = c("<none>", "features", "phenotypes"), width = "100%"
        ),
        "Overlay a biplot: arrows for the CpGs that define each component, or for the direction in which each sample variable increases.",
        placement = "left"
      )
    ),
    shiny::conditionalPanel(
      "input.pca_arrows != '<none>' && input.clustmethod == 'pca'",
      ns = ns,
      withTooltip(
        shiny::numericInput(ns("pca_numarrows"), "Number of arrows:",
          value = 3, min = 1, max = 50, step = 1, width = "100%"
        ),
        "How many arrows to label and draw, ranked by length. Lower numbers keep the plot readable.",
      placement = "left"
      ),
      withTooltip(
        shiny::checkboxInput(ns("pca_repel"), "repel labels", value = TRUE),
        "Move overlapping arrow labels apart. Switch off to pin each label to its own arrow tip.",
      placement = "left"
      )
    ),
    withTooltip(
      shiny::checkboxInput(ns("all_clustmethods"), "show all methods"),
      "Show an overview of all dimensionality reduction methods.",
      placement = "left"
    ),
    withTooltip(
      shiny::checkboxInput(ns("plot3d"), "plot 3D"),
      "Show 3D plot.",
      placement = "left"
    )
  )

  PlotModuleUI(
    translate = FALSE,
    ns("pltmod"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    options = plot_opts,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "clustering_prism"
  )
}

methylome_plot_clustpca_server <- function(id,
                                            pgx,
                                            selected_samples,
                                            clustmethod,
                                            pca_components,
                                            pca_dims,
                                            labeltype = shiny::reactive("feature"),
                                            watermark = FALSE,
                                            parent) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## The panel's own Layout control wins when it exists; the clustmethod
    ## argument remains the fallback, so a caller can still pin the layout and
    ## the signature stays the one upstream uses.
    layout_method <- shiny::reactive({
      m <- input$clustmethod
      if (is.null(m) || !nzchar(m)) clustmethod() else m
    })

    ## Editor: dynamic color pickers for custom palette
    output$custom_palette_ui <- shiny::renderUI({
      shiny::req(input$palette == "custom")
      colvar <- input$hmpca.colvar
      shiny::req(colvar)
      samples <- selected_samples()
      groups <- sort(unique(as.character(pgx$samples[samples, colvar])))
      custom_palette_pickers(groups, ns)
    })

    plot_data <- shiny::reactive({
      samples <- selected_samples()
      ## drop pca2d.varexp and friends: only the position matrices are tabular
      cluster.pos <- Filter(is.matrix, pgx$cluster$pos)
      ## the board recomputes PCA on the fly, so export those instead of the
      ## stored PC1/PC2 only
      cluster.pos <- cluster.pos[grep("^pca", names(cluster.pos), invert = TRUE)]
      cluster.pos$pca <- pca_components()$pos
      for (m in names(cluster.pos)) {
        colnames(cluster.pos[[m]]) <- paste0(m, ".", colnames(cluster.pos[[m]]))
      }
      all.pos <- do.call(cbind, cluster.pos)
      all.pos <- all.pos[samples, ]
      pd <- list(pos = all.pos)
      return(pd)
    })


    ## arrows scaled to fill the score cloud (loadings and correlations live on
    ## a much smaller scale than the scores), with the label just past the tip
    scale_arrows <- function(arrows, pos) {
      s <- 0.65 * max(abs(pos)) / max(abs(arrows))
      ax <- arrows[, 1] * s
      ay <- arrows[, 2] * s
      r <- sqrt(ax^2 + ay^2)
      r[r == 0] <- 1
      nudge <- 0.05 * max(abs(pos))
      list(
        x = ax, y = ay,
        lx = ax + nudge * ax / r, ly = ay + nudge * ay / r,
        text = rownames(arrows)
      )
    }

    create_plot <- function(pgx, pos, xlab, ylab, method, colvar, shapevar, label, cex, palette = "muted_light", custom_colors = NULL, arrows = NULL) {
      do3d <- (ncol(pos) == 3)
      sel <- rownames(pos)
      df <- cbind(pos, pgx$samples[sel, , drop = FALSE])

      textvar <- NULL
      if (colvar %in% colnames(df)) {
        colvar <- factor(df[, colvar])
      }

      if (shapevar %in% colnames(df)) {
        shapevar <- factor(df[, shapevar])
      }

      ann.text <- rep(" ", nrow(df))

      label.samples <- (label == "sample")

      if (!do3d && label.samples) ann.text <- rownames(df)
      if (!is.null(colvar)) textvar <- factor(colvar)

      symbols <- c(
        "circle", "square", "star", "triangle-up", "triangle-down", "pentagon",
        "bowtie", "hexagon", "asterisk", "hash", "cross", "triangle-left",
        "triangle-right", "+", c(15:0)
      )

      Y <- cbind("sample" = rownames(pos), pgx$Y[sel, ])
      tt.info <- apply(Y, 1, function(y) paste0(colnames(Y), ": <b>", y, "</b></br>", collapse = ""))
      tt.info <- I(as.character(tt.info))
      cex1 <- c(1.0, 0.8, 0.6)[1 + 1 * (nrow(pos) > 30) + 1 * (nrow(pos) > 200)]
      clrs.length <- length(unique(colvar))
      if (!is.null(custom_colors) && length(custom_colors) >= clrs.length) {
        clrs <- custom_colors[1:clrs.length]
      } else {
        clrs <- rep(omics_pal_d(palette = palette)(8), ceiling(clrs.length / 8))[1:clrs.length]
      }

      if (do3d) {
        plt <- plotly::plot_ly(df, mode = "markers") %>%
          plotly::add_markers(
            x = df[, 1],
            y = df[, 2],
            z = df[, 3],
            type = "scatter3d",
            color = colvar,
            colors = clrs,
            marker = list(
              size = 6 * cex1 * cex,
              line = list(color = "grey10", width = 0.1)
            ),
            symbol = shapevar,
            symbols = symbols,
            text = tt.info
          ) %>%
          plotly::add_annotations(
            x = pos[, 1],
            y = pos[, 2],
            z = pos[, 3],
            text = ann.text,
            showarrow = FALSE
          )
        ## add cluster annotation labels
        if (0 && length(unique(colvar)) > 1) {
          ## add cluster annotation labels
          grp.pos <- apply(pos, 2, function(x) tapply(x, colvar, median))
          cex2 <- ifelse(length(grp.pos) > 20, 0.8, 1)
          plt <- plt %>% plotly::add_annotations(
            x = grp.pos[, 1], y = grp.pos[, 2], z = grp.pos[, 3],
            text = rownames(grp.pos),
            font = list(size = 24 * cex2 * cex, color = "#555"),
            showarrow = FALSE
          )
        }
        if (label == "<none>") {
          plt <- plt %>%
            plotly::layout(showlegend = FALSE)
        }
      } else {
        plt <- plotly::plot_ly(
          df,
          mode = "markers",
          hovertemplate = "</br>%{text}<extra></extra>"
        ) %>%
          plotly::add_markers(
            x = df[, 1],
            y = df[, 2],
            type = "scattergl",
            color = colvar,
            colors = clrs,
            marker = list(
              size = 16 * cex1 * cex,
              line = list(color = "grey20", width = 0.6)
            ),
            symbol = shapevar,
            symbols = symbols,
            text = tt.info
          )

        ## Sample labels go in a text TRACE, not annotations: plotly.repel
        ## replaces the whole annotation array in the browser, which would
        ## silently erase them whenever arrow labels are being repelled.
        if (label.samples) {
          plt <- plt %>% plotly::add_text(
            x = pos[, 1], y = pos[, 2], text = rownames(df),
            textposition = "top center",
            textfont = list(size = 12 * cex, color = unname(omics_colors("super_dark_grey"))),
            showlegend = FALSE, hoverinfo = "skip", inherit = FALSE
          )
        }

        plt <- plt %>% plotly::layout(
          xaxis = list(title = xlab),
          yaxis = list(title = ylab)
        )

        if (!is.null(arrows) && nrow(arrows) && max(abs(arrows)) > 0) {
          aa <- scale_arrows(arrows, pos)
          ## saturated accent: the points run on muted_light (desaturated and
          ## lightened), so the biplot overlay reads as an overlay and not as
          ## another group, and its labels stay apart from the grey sample ones
          arrow.clr <- unname(omics_colors("terra_cotta"))
          ## derived, not written out: a literal rgba() would go stale the day
          ## terra_cotta changes
          arrow.seg <- sprintf(
            "rgba(%s,0.6)", paste(grDevices::col2rgb(arrow.clr), collapse = ",")
          )
          ## Arrows are drawn as traces, NOT as annotations: plotly.repel takes
          ## ownership of layout.annotations at render time (plotly-repel.js
          ## relayouts the whole array), which wipes anything annotation-drawn.
          ## Shafts are one NA-separated line trace, heads are markers rotated
          ## along each arrow (angle is the bearing, clockwise from up).
          plt <- plt %>%
            plotly::add_trace(
              x = as.vector(rbind(0, aa$x, NA)),
              y = as.vector(rbind(0, aa$y, NA)),
              type = "scatter", mode = "lines",
              line = list(color = arrow.clr, width = 1.8),
              showlegend = FALSE, hoverinfo = "skip", inherit = FALSE
            ) %>%
            plotly::add_trace(
              x = aa$x, y = aa$y,
              type = "scatter", mode = "markers",
              marker = list(
                symbol = "arrow", size = 12, color = arrow.clr,
                angle = atan2(aa$x, aa$y) * 180 / pi
              ),
              text = aa$text, hoverinfo = "text",
              showlegend = FALSE, inherit = FALSE
            )
          ## Labels are ALWAYS drawn as annotations, nudged just past their own
          ## tip. plotly.repel is then layered on top when asked for: its JS
          ## relayouts the whole annotation array in the browser, so on screen
          ## the repelled labels REPLACE these rather than doubling up.
          ## Wherever that JS does not run -- kaleido image export, or
          ## plotly::subplot() stripping the render hook in "show all methods"
          ## -- the annotations survive, so labels are never silently lost.
          plt <- plt %>% plotly::add_annotations(
            x = aa$lx, y = aa$ly, text = aa$text,
            showarrow = FALSE,
            font = list(color = arrow.clr, size = 13 * cex),
            bgcolor = "rgba(255,255,255,0.65)", borderpad = 2
          )
          if (isTRUE(input$pca_repel) && requireNamespace("plotly.repel", quietly = TRUE)) {
            plt <- plotly.repel::add_text_repel(
              plt,
              data = data.frame(x = aa$x, y = aa$y, text = aa$text),
              x = ~x, y = ~y, text = ~text,
              box_padding = 0.4, point_padding = 0.3, max_overlaps = 20,
              font = list(color = arrow.clr, size = 13 * cex),
              ## leader lines are connective tissue, not data: same hue, faded
              segment = list(color = arrow.seg, width = 1)
            )
          }
        }

        ## add group/cluster annotation labels
        if (label == "inside") {
          plt <- plt %>%
            plotly::layout(legend = list(x = 0.05, y = 0.95))
        } else if (label == "bottom") {
          plt <- plt %>%
            plotly::layout(legend = list(orientation = "h"))
        } else if (label == "group") {
          if (!is.null(textvar) && length(unique(textvar)) > 1) {
            grp.pos <- apply(pos, 2, function(x) tapply(x, as.character(textvar), median))
            cex2 <- 1
            if (length(grp.pos) > 20) cex2 <- 0.8
            if (length(grp.pos) > 50) cex2 <- 0.6
            ## text trace, not annotations: see the sample-label note above
            plt <- plt %>% plotly::add_text(
              x = grp.pos[, 1],
              y = grp.pos[, 2],
              text = paste0("<b>", rownames(grp.pos), "</b>"),
              textfont = list(size = 24 * cex2 * cex, color = "#555"),
              showlegend = FALSE, hoverinfo = "skip", inherit = FALSE
            )
          }
          plt <- plt %>%
            plotly::layout(showlegend = FALSE)
        } else if (label == "sample") {
          plt <- plt %>%
            plotly::layout(showlegend = FALSE)
        } else if (label == "<none>") {
          plt <- plt %>%
            plotly::layout(showlegend = FALSE)
        }
      }
      return(plt)
    }

    create_ggplot <- function(pgx, pos, xlab, ylab, method, colvar, label, cex, palette = "muted_light", custom_colors = NULL, arrows = NULL) {
      sel <- rownames(pos)
      df <- cbind(pos, pgx$samples[sel, , drop = FALSE])

      if (colvar %in% colnames(df)) {
        colvar_fct <- factor(df[, colvar])
      } else {
        colvar_fct <- factor(rep("all", nrow(df)))
      }

      clrs.length <- length(unique(colvar_fct))
      if (!is.null(custom_colors) && length(custom_colors) >= clrs.length) {
        clrs <- custom_colors[1:clrs.length]
      } else {
        clrs <- rep(omics_pal_d(palette = palette)(8), ceiling(clrs.length / 8))[1:clrs.length]
      }

      showlabels <- (label %in% c("group", "inside"))

      p <- playbase::pgx.scatterPlotXY.GGPLOT(
        pos,
        var = colvar_fct,
        col = clrs,
        cex = cex,
        xlab = xlab,
        ylab = ylab,
        title = toupper(method),
        cex.title = 1.2,
        cex.clust = 1.1,
        label.clusters = showlabels,
        legend = !(label %in% c("<none>", "group", "sample"))
      )

      if (!is.null(arrows) && nrow(arrows) && max(abs(arrows)) > 0) {
        aa <- scale_arrows(arrows, pos)
        df.a <- data.frame(x = aa$x, y = aa$y, lx = aa$lx, ly = aa$ly, text = aa$text)
        ## same accent as the plotly path, see the note there
        arrow.clr <- unname(omics_colors("terra_cotta"))
        p <- p +
          ggplot2::geom_segment(
            data = df.a,
            ggplot2::aes(x = 0, y = 0, xend = x, yend = y),
            arrow = ggplot2::arrow(length = ggplot2::unit(0.18, "cm"), type = "closed"),
            colour = arrow.clr, linewidth = 0.5, inherit.aes = FALSE
          ) +
          if (isTRUE(input$pca_repel)) {
            ggrepel::geom_text_repel(
              data = df.a,
              ggplot2::aes(x = x, y = y, label = text),
              colour = arrow.clr, size = 3 * cex, inherit.aes = FALSE,
              box.padding = 0.4, point.padding = 0.3,
              segment.colour = arrow.clr, segment.alpha = 0.6, max.overlaps = 20
            )
          } else {
            ggplot2::geom_text(
              data = df.a,
              ggplot2::aes(x = lx, y = ly, label = text),
              colour = arrow.clr, size = 3 * cex, inherit.aes = FALSE
            )
          }
      }
      p
    }

    ## Scree data: percentage of variance per component plus its running total.
    ## Both series are percentages, so they share one axis. The percentages come
    ## from pca_components() and are relative to the total variance, so they do
    ## not sum to 100 over the five components we compute -- that is correct.
    ## Biplot arrow endpoints on the two selected components, in loading or
    ## correlation units. Ranked by length and cut to the top N; the caller
    ## scales them to the score cloud.
    arrow_ends <- function(pc, k) {
      mode <- input$pca_arrows
      if (is.null(mode) || mode == "<none>") {
        return(NULL)
      }
      A <- if (mode == "features") pc$loadings else pc$pheno.cor
      if (is.null(A) || max(k) > ncol(A)) {
        return(NULL)
      }
      A <- A[stats::complete.cases(A[, k, drop = FALSE]), k, drop = FALSE]
      if (!nrow(A)) {
        return(NULL)
      }
      ## free-text numeric: empty or nonsense while typing must not blank the plot
      n <- suppressWarnings(as.integer(input$pca_numarrows))
      if (length(n) != 1 || is.na(n)) n <- 3
      n <- min(max(1, n), 50) ## the UI max is advisory: typed values bypass it
      A <- A[head(order(-(A[, 1]^2 + A[, 2]^2)), n), , drop = FALSE]
      if (mode == "features") {
        ## probe2symbol returns NULL when the requested label column is absent
        ## (reachable while labeltype still holds the previous dataset's
        ## choice). Blanking the rownames would blank the labels downstream.
        sym <- playbase::probe2symbol(rownames(A), pgx$genes, labeltype(), fill_na = TRUE)
        ## Keep the probe id and append the gene only when it adds something:
        ## the CpG is the unit that was measured, and two probes in the same
        ## gene routinely move in opposite directions.
        if (length(sym) == nrow(A)) {
          sym <- sub(";.*", "", as.character(sym))
          add <- !is.na(sym) & nzchar(sym) & sym != rownames(A)
          rownames(A)[add] <- sprintf("%s (%s)", rownames(A)[add], sym[add])
        }
      }
      A
    }

    ## Positions and axis labels for one layout method. PCA is recomputed on
    ## the fly (pgx only stores PC1-PC3) so any pair of PC1-PC5 can be plotted.
    ## t-SNE and UMAP off the SAME components as everything else, computed
    ## here rather than read from pgx$cluster$pos.
    ##
    ## This is what the platform does too - pgx.clusterBigMatrix embeds its
    ## 50-component reduction, not the raw matrix - so the only thing that
    ## changes is which decomposition feeds it. That matters: the stored
    ## coordinates are frozen at upload against the top 1,000 probes, so with
    ## "Use all features" on, switching layout used to drop the user silently
    ## back to a 1,000-probe view of the same cohort.
    ##
    ## Seeded, because both are stochastic and a plot that reshuffles on every
    ## redraw cannot be compared with the one in someone's slide deck. The
    ## stored coordinates additionally carry a 1e-3 jitter that nothing seeds.
    embeddings <- shiny::reactive({
      X <- pca_components()$pos
      n <- nrow(X)
      ## Rtsne requires 3 * perplexity < n - 1; playbase lands near (n-1)/4.
      perp <- max(2, min(30, floor((n - 1) / 4)))
      emb <- function(f) tryCatch(f(), error = function(e) NULL)
      list(
        tsne2d = emb(function() { set.seed(42)
          Rtsne::Rtsne(X, dims = 2, perplexity = perp, pca = FALSE,
                       check_duplicates = FALSE)$Y }),
        tsne3d = emb(function() { set.seed(42)
          Rtsne::Rtsne(X, dims = 3, perplexity = perp, pca = FALSE,
                       check_duplicates = FALSE)$Y }),
        umap2d = emb(function() { set.seed(42)
          uwot::umap(X, n_components = 2, n_neighbors = max(2, min(15, n - 1))) }),
        umap3d = emb(function() { set.seed(42)
          uwot::umap(X, n_components = 3, n_neighbors = max(2, min(15, n - 1))) })
      )
    })

    method_pos <- function(m, samples, do3d = FALSE) {
      arrows <- NULL
      pc <- pca_components()
      if (m == "pca") {
        k <- if (do3d) seq_len(min(3, ncol(pc$pos))) else pca_dims()
        pos <- pc$pos[samples, k, drop = FALSE]
        labs <- paste0("PC", k, " (", pc$varexp[k], "%)")
        if (!do3d) arrows <- arrow_ends(pc, k)
      } else {
        m1 <- paste0(m, ifelse(do3d, "3d", "2d"))
        pos <- embeddings()[[m1]]
        shiny::validate(shiny::need(!is.null(pos),
          sprintf("%s could not be computed for this cohort - too few samples for its neighbourhood size.",
                  toupper(m))))
        rownames(pos) <- rownames(pc$pos)
        pos <- pos[samples, , drop = FALSE]
        labs <- paste0(toupper(m), seq_len(ncol(pos)))
      }
      list(pos = pos, xlab = labs[1], ylab = labs[2], arrows = arrows)
    }

    create_plotlist <- function() {
      samples <- selected_samples()
      options <- input$hmpca_options
      colvar <- input$hmpca.colvar
      shapevar <- input$hmpca.shapevar
      clustmethod <- layout_method()
      label <- input$pca_label

      shiny::validate(shiny::need(
        length(samples) > 1,
        "Filtering too restrictive. Please change 'Filter samples' settings."
      ))
      shiny::req(samples, colvar, shapevar, clustmethod, legend)

      methods <- layout_method()
      ## All three are computed in this panel now, so the list is fixed rather
      ## than read off whatever the pgx happens to store.
      if (input$all_clustmethods) methods <- c("pca", "tsne", "umap")
      do3d <- (input$plot3d)
      multiplot <- length(methods) > 1

      ## Editor: palette and custom colors
      groups <- sort(unique(as.character(pgx$samples[samples, colvar])))
      n_groups <- length(groups)
      custom_colors <- resolve_palette_colors(input, n_groups, fallback_colors = omics_pal_d("muted_light")(n_groups))
      palette <- if (!is.null(input$palette)) input$palette else "muted_light"
      if (palette %in% c("original", "default", "custom", "")) palette <- "muted_light"

      plist <- list()
      for (i in 1:length(methods)) {
        m <- methods[i]
        mp <- method_pos(m, samples, do3d = do3d)
        plist[[i]] <- create_plot(
          pgx = pgx,
          pos = mp$pos,
          xlab = mp$xlab,
          ylab = mp$ylab,
          method = m,
          colvar = colvar,
          shapevar = shapevar,
          label = label,
          cex = ifelse(length(methods) > 1, 0.6, 1),
          palette = palette,
          custom_colors = custom_colors,
          arrows = mp$arrows
        )
      }
      plist
    }

    create_gg_plotlist <- function() {
      samples <- selected_samples()
      colvar <- input$hmpca.colvar
      clustmethod <- layout_method()
      label <- input$pca_label

      shiny::validate(shiny::need(
        length(samples) > 1,
        "Filtering too restrictive. Please change 'Filter samples' settings."
      ))
      shiny::req(samples, colvar, clustmethod)

      methods <- layout_method()
      ## All three are computed in this panel now, so the list is fixed rather
      ## than read off whatever the pgx happens to store.
      if (input$all_clustmethods) methods <- c("pca", "tsne", "umap")

      groups <- sort(unique(as.character(pgx$samples[samples, colvar])))
      n_groups <- length(groups)
      custom_colors <- resolve_palette_colors(input, n_groups, fallback_colors = omics_pal_d("muted_light")(n_groups))
      palette <- if (!is.null(input$palette)) input$palette else "muted_light"
      if (palette %in% c("original", "default", "custom", "")) palette <- "muted_light"

      gp <- extract_ggprism_params(input)

      plts <- list()
      for (i in seq_along(methods)) {
        m <- methods[i]
        mp <- method_pos(m, samples)

        p <- create_ggplot(
          pgx = pgx,
          pos = mp$pos,
          xlab = mp$xlab,
          ylab = mp$ylab,
          method = m,
          colvar = colvar,
          label = label,
          cex = ifelse(length(methods) > 1, 0.6, 1),
          palette = palette,
          custom_colors = custom_colors,
          arrows = mp$arrows
        )
        p <- apply_ggprism_theme(p, gp)
        p <- apply_editor_theme(p, input)
        plts[[length(plts) + 1]] <- p
      }
      plts
    }

    plot.RENDER <- reactive({
      gp <- extract_ggprism_params(input)
      do3d <- isTRUE(input$plot3d)

      if (gp$use_ggprism && !do3d) {
        plts <- create_gg_plotlist()
        shiny::req(length(plts) > 0)
        if (length(plts) == 1) {
          ggplot_as_plotly_image(plts[[1]], width = 6, height = 6)
        } else {
          nc <- ceiling(sqrt(length(plts)))
          combined <- patchwork::wrap_plots(plts, ncol = nc)
          nr <- ceiling(length(plts) / nc)
          ggplot_as_plotly_image(combined, width = nc * 4, height = nr * 4)
        }
      } else {
        fig <- if (length(create_plotlist()) == 1) {
          create_plotlist()[[1]]
        } else {
          plist <- create_plotlist()
          nc <- ceiling(sqrt(length(plist)))
          plotly::subplot(plist, nrows = nc, margin = 0.04)
        }
        apply_plotly_editor_theme(fig, input)
      }
    })

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 170), ## resolution of plots
      pdf.width = 8,
      pdf.height = 8,
      add.watermark = watermark,
      parent_session = session
    )
  })
}
