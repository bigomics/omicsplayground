##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## EWAS. The model is fitted here rather than read from the pgx, so the user
## can adjust for cell composition, batch, age and anything else on the sample
## sheet - see methylome_model.R for why that is not optional.
##
## The settings panel lives in methylome_ewas_inputs.R.


## ------------------------------------------------------------------ shared --

## Per-sample group labels for a contrast. pgx$contrasts holds the actual
## labels; exp.matrix only holds signed weights, so prefer the former.
mp_contrast_groups <- function(pgx, cmp) {
  ct <- pgx$contrasts
  if (!is.null(ct) && cmp %in% colnames(ct)) {
    g <- as.character(ct[, cmp]); names(g) <- rownames(ct)
    g[!is.na(g) & g != ""]
  } else {
    em <- pgx$model.parameters$exp.matrix
    if (is.null(em) || !cmp %in% colnames(em)) return(NULL)
    w <- em[, cmp]
    stats::setNames(ifelse(w > 0, "group 1", "group 0"), names(w))[w != 0]
  }
}

## Attach the array annotation to a fitted table. pgx$genes$chr is a cytoband
## ("16q12.2"), not a bare chromosome, so the arm and band are stripped as
## board.epigenomics does.
mp_ewas_annotate <- function(pgx, tab) {
  g <- pgx$genes[tab$probe, , drop = FALSE]
  cc <- grep("^chr$|chrom", tolower(colnames(g)))[1]
  pc <- grep("^pos$|position", tolower(colnames(g)))[1]
  isl <- mp_context_col(pgx, "island")
  loc <- mp_context_col(pgx, "gene")
  idx <- match(tab$probe, rownames(pgx$genes))
  tab$gene <- if ("symbol" %in% colnames(g)) as.character(g$symbol) else NA
  tab$chr  <- if (!is.na(cc)) sub("(p|q|cen).*", "", sub("^chr", "", as.character(g[[cc]]))) else NA
  tab$pos  <- if (!is.na(pc)) suppressWarnings(as.numeric(sub(";.*", "", g[[pc]]))) else NA
  tab$island   <- if (!is.null(isl)) sub(";.*", "", isl[idx]) else NA
  tab$location <- if (!is.null(loc)) sub(";.*", "", loc[idx]) else NA
  tab
}

mp_ewas_sig <- function(d, thresh) {
  stat <- if (identical(thresh$type, "p")) d$p else d$q
  ok <- !is.na(stat) & stat <= thresh$value
  if (!is.null(thresh$min_dbeta) && thresh$min_dbeta > 0) {
    ok <- ok & !is.na(d$dbeta) & abs(d$dbeta) >= thresh$min_dbeta
  }
  ok
}

## A q-value cut-off has no fixed position on a -log10(p) axis, so use the
## largest p that still passes; with no hits there is no line.
mp_ewas_hline <- function(d, thresh) {
  if (identical(thresh$type, "p")) return(-log10(thresh$value))
  sig <- !is.na(d$q) & d$q <= thresh$value
  if (!any(sig)) return(NA_real_)
  -log10(max(d$p[sig], na.rm = TRUE))
}

mp_lambda <- function(p) {
  stats::median(stats::qchisq(1 - p, 1), na.rm = TRUE) / stats::qchisq(0.5, 1)
}

## Empirical-null bias and inflation. Applied as a diagnostic only: bacon on
## an unadjusted model would launder confounding into well-calibrated-looking
## p-values, which is worse than an honest lambda.
mp_bacon <- function(tstat) {
  if (!requireNamespace("bacon", quietly = TRUE)) return(NULL)
  tryCatch({
    bc <- bacon::bacon(teststatistics = as.numeric(tstat[is.finite(tstat)]))
    list(bias = as.numeric(bacon::bias(bc)),
         inflation = as.numeric(bacon::inflation(bc)))
  }, error = function(e) NULL)
}

## ---------------------------------------------------------------- Manhattan --

methylome_plot_manhattan_ui <- function(id, title, caption, info.text, info.methods,
                                        height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "a")
}

methylome_plot_manhattan_server <- function(id, r.ewas, r.thresh, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      d <- r.ewas()$data
      th <- r.thresh()
      d$sig <- mp_ewas_sig(d, th)
      d <- d[!is.na(d$chr) & !is.na(d$pos) & is.finite(d$pos) & is.finite(d$p), ]
      shiny::validate(shiny::need(nrow(d) > 0, "No probes carry chromosome and position."))
      chrs <- intersect(c(as.character(1:22), "X", "Y"), unique(d$chr))
      d <- d[d$chr %in% chrs, ]
      d$chr <- factor(d$chr, levels = chrs)
      d <- d[order(d$chr, d$pos), ]
      mx <- vapply(split(d$pos, d$chr), function(z) if (length(z)) max(z) else 0, numeric(1))
      off <- c(0, utils::head(cumsum(as.numeric(mx) * 1.02), -1)); names(off) <- chrs
      d$gx <- d$pos + off[as.character(d$chr)]
      op <- graphics::par(mar = c(4, 4.5, 1, 1), las = 1); on.exit(graphics::par(op))
      cols <- ifelse(as.integer(d$chr) %% 2 == 0, "#a9b6c4", "#7f8f9f")
      plot(d$gx, -log10(d$p), pch = 19, cex = .35, col = cols, xaxt = "n",
           xlab = "chromosome", ylab = expression(-log[10](p)))
      graphics::points(d$gx[d$sig], -log10(d$p[d$sig]), pch = 19, cex = .6,
        col = ifelse(d$dbeta[d$sig] > 0 | is.na(d$dbeta[d$sig]), MP_PAL$hyper, MP_PAL$hypo))
      graphics::axis(1, at = tapply(d$gx, d$chr, stats::median, na.rm = TRUE),
                     labels = levels(d$chr), tick = FALSE, cex.axis = .7)
      h <- mp_ewas_hline(d, th)
      if (is.finite(h)) {
        graphics::abline(h = h, lty = 2, col = MP_PAL$grey)
        graphics::mtext(sprintf("%s <= %g  (%s hits)", th$type, th$value,
                                format(sum(d$sig), big.mark = ",")),
                        side = 3, adj = 1, cex = .75, col = MP_PAL$grey)
      } else {
        graphics::mtext(sprintf("no probe passes %s <= %g", th$type, th$value),
                        side = 3, adj = 1, cex = .75, col = MP_PAL$grey)
      }
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive(r.ewas()$data), renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot, res = c(90, 130),
      pdf.width = 10, pdf.height = 5, add.watermark = watermark)
  })
}

## ----------------------------------------------------------------------- QQ --

methylome_plot_qq_ui <- function(id, title, caption, info.text, info.methods,
                                 height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "b")
}

methylome_plot_qq_server <- function(id, r.ewas, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      res <- r.ewas(); d <- res$data
      pv <- sort(d$p[!is.na(d$p)])
      shiny::req(length(pv) > 0)
      op <- graphics::par(mar = c(4, 4.5, 1, 1), las = 1); on.exit(graphics::par(op))
      plot(-log10(stats::ppoints(length(pv))), -log10(pv), pch = 19, cex = .35,
           col = MP_PAL$grey, xlab = "expected", ylab = "observed")
      graphics::abline(0, 1, lty = 2, col = MP_PAL$hypo)
      lab <- sprintf("lambda = %.3f", mp_lambda(d$p))
      if (!is.null(res$bacon)) {
        lab <- c(lab,
                 sprintf("bacon inflation = %.3f", res$bacon$inflation),
                 sprintf("bacon bias = %+.3f", res$bacon$bias))
      }
      graphics::legend("topleft", bty = "n", legend = lab)
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive(r.ewas()$data), renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot, res = c(90, 130),
      pdf.width = 6, pdf.height = 6, add.watermark = watermark)
  })
}

## --------------------------------------------------------- hit context enrich --

methylome_plot_enrichment_ui <- function(id, title, caption, info.text, info.methods,
                                         height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "c")
}

methylome_plot_enrichment_server <- function(id, r.ewas, r.thresh, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      d <- r.ewas()$data
      d$sig <- mp_ewas_sig(d, r.thresh())
      d <- d[!is.na(d$island), ]
      sig <- d$sig
      shiny::validate(shiny::need(sum(sig) >= 5,
        "Too few hits at this threshold to test context enrichment."))
      lv <- intersect(c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"),
                      unique(d$island))
      ## Background is the probes actually tested, not the whole array.
      or <- vapply(lv, function(k) {
        a <- sum(sig & d$island == k); b <- sum(sig) - a
        c1 <- sum(!sig & d$island == k); d1 <- sum(!sig) - c1
        ((a + .5) / (b + .5)) / ((c1 + .5) / (d1 + .5))
      }, numeric(1))
      data.frame(context = lv, odds_ratio = as.numeric(or),
                 log2_or = log2(as.numeric(or)), stringsAsFactors = FALSE)
    })
    plot.RENDER <- function() {
      e <- plot_data(); shiny::req(nrow(e) > 0)
      op <- graphics::par(mar = c(4, 6.5, 0.5, 1.5), las = 1, mgp = c(2.2, 0.6, 0))
      on.exit(graphics::par(op))
      cols <- ifelse(e$log2_or >= 0, MP_PAL$hyper, MP_PAL$hypo)
      graphics::barplot(rev(e$log2_or), horiz = TRUE, col = rev(cols), border = NA,
        names.arg = rev(e$context), cex.names = 0.85, cex.axis = 0.85,
        xlab = "log2 odds ratio vs tested background")
      graphics::abline(v = 0, col = "#444444")
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = plot_data, renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot, res = c(90, 130),
      pdf.width = 6, pdf.height = 5, add.watermark = watermark)
  })
}

## ------------------------------------------------------------- hits table ---

methylome_table_hits_ui <- function(id, title, info.text, caption, height, width) {
  ns <- shiny::NS(id)
  TableModuleUI(ns("tblmod"), title = title, info.text = info.text,
                caption = caption, height = height, width = width, label = "d")
}

methylome_table_hits_server <- function(id, r.ewas, r.thresh, scrollY = "22vh") {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      d <- r.ewas()$data
      sig <- mp_ewas_sig(d, r.thresh())
      shiny::validate(shiny::need(any(sig),
        "No CpG passes the current threshold. Relax the cut-off in the settings panel."))
      d <- d[sig, , drop = FALSE]
      d <- d[order(d$p), , drop = FALSE]
      out <- data.frame(
        CpG = d$probe,
        Gene = ifelse(is.na(d$gene) | d$gene == "", "-", d$gene),
        Chr = d$chr, Position = d$pos,
        Context = ifelse(is.na(d$island), "-", d$island),
        Region = ifelse(is.na(d$location) | d$location == "", "-", d$location),
        `Delta beta` = d$dbeta,
        Direction = ifelse(is.na(d$dbeta), "-", ifelse(d$dbeta > 0, "hyper", "hypo")),
        `P value` = signif(d$p, 3), FDR = signif(d$q, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      rownames(out) <- NULL
      out
    })
    render <- function(sy) {
      dt <- table_data(); shiny::req(dt)
      DT::datatable(dt, class = "compact hover", rownames = FALSE,
        extensions = c("Buttons", "Scroller"), plugins = "scrollResize",
        selection = list(mode = "single", target = "row"),
        options = list(dom = "lfrtip", scroller = TRUE, scrollX = TRUE,
                       scrollY = sy, scrollResize = TRUE, deferRender = TRUE,
                       order = list(list(8, "asc")))) |>
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") |>
        DT::formatStyle("Direction",
          color = DT::styleEqual(c("hyper", "hypo"), c(MP_PAL$hyper, MP_PAL$hypo)),
          fontWeight = "bold") |>
        DT::formatStyle("Delta beta",
          background = DT::styleColorBar(c(-1, 1) * max(abs(dt$`Delta beta`), na.rm = TRUE),
                                         "#dce7f2"),
          backgroundSize = "98% 60%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center")
    }
    TableModuleServer("tblmod", func = function() render(scrollY),
                      func2 = function() render("60vh"),
                      csvFunc = table_data, selector = "none")
  })
}

## ------------------------------------------------------- per-CpG stripcharts --
##
## The panel the field expects next to a Manhattan (plotCpg in the Bioconductor
## workflow, meffil.ewas.cpg.plot in meffil): the only view showing whether a
## hit is outlier-driven and at what absolute methylation level it sits.

methylome_plot_stripcharts_ui <- function(id, title, caption, info.text, info.methods,
                                          height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "e")
}

methylome_plot_stripcharts_server <- function(id, pgx, r.ewas, r.thresh,
                                              r.topn = shiny::reactive(6),
                                              watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      p <- pgx(); res <- r.ewas()
      d <- res$data
      d <- d[mp_ewas_sig(d, r.thresh()), , drop = FALSE]
      shiny::validate(shiny::need(nrow(d) > 0,
        "No CpG passes the current threshold. Relax the cut-off in the settings panel."))
      n <- r.topn(); if (is.null(n) || is.na(n) || n < 1) n <- 6
      d <- d[order(d$p), , drop = FALSE][seq_len(min(n, nrow(d))), , drop = FALSE]
      grp <- mp_contrast_groups(p, res$contrast)
      shiny::validate(shiny::need(!is.null(grp), "No group labels for this contrast."))
      X <- mp_beta(p)
      ss <- intersect(names(grp), colnames(X))
      list(d = d, X = X[d$probe, ss, drop = FALSE], grp = factor(grp[ss]))
    })
    plot.RENDER <- function() {
      res <- plot_data(); d <- res$d; X <- res$X; grp <- res$grp
      n <- nrow(d); nc <- min(3, n); nr <- ceiling(n / nc)
      op <- graphics::par(mfrow = c(nr, nc), mar = c(2.6, 3.2, 2.2, 0.6),
                          las = 1, mgp = c(2, 0.6, 0), cex.axis = 0.85)
      on.exit(graphics::par(op))
      lv <- levels(grp); cols <- c("#4575b4", "#d73027")[seq_along(lv)]
      for (i in seq_len(n)) {
        y <- as.numeric(X[i, ])
        graphics::plot(NA, xlim = c(0.5, length(lv) + 0.5), ylim = c(0, 1),
                       xaxt = "n", xlab = "", ylab = "beta")
        graphics::axis(1, at = seq_along(lv), labels = lv, tick = FALSE)
        graphics::abline(h = c(0.2, 0.8), lty = 3, col = "#c9ccc7")
        for (j in seq_along(lv)) {
          yy <- y[grp == lv[j]]
          graphics::points(jitter(rep(j, length(yy)), amount = 0.13), yy,
                           pch = 19, cex = 0.7,
                           col = grDevices::adjustcolor(cols[j], 0.65))
          graphics::segments(j - 0.28, mean(yy, na.rm = TRUE),
                             j + 0.28, mean(yy, na.rm = TRUE),
                             lwd = 2.5, col = cols[j])
        }
        lab <- if (!is.na(d$gene[i]) && d$gene[i] != "") {
          paste0(d$probe[i], "  (", d$gene[i], ")")
        } else d$probe[i]
        graphics::title(main = lab, cex.main = 0.95, font.main = 1)
        graphics::mtext(sprintf("d-beta %+.3f   q %.2g", d$dbeta[i], d$q[i]),
                        side = 3, line = 0.1, cex = 0.62, col = "#697586")
      }
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive({
        r <- plot_data(); cbind(probe = rownames(r$X), as.data.frame(r$X))
      }),
      renderFunc = shiny::renderPlot, renderFunc2 = shiny::renderPlot,
      res = c(90, 130), pdf.width = 9, pdf.height = 6, add.watermark = watermark)
  })
}
