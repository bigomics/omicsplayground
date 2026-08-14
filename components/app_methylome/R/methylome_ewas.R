##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## EWAS. The only screen here that reads pgx$gx.meta - the differential
## methylation already computed at upload, which no methylation board has
## displayed before.

## ---------------------------------------------------------------- settings --

## Threshold drives the Manhattan line, which probes count as hits, the
## context enrichment and the hit table, so it lives in the tab settings
## rather than as a per-plot option.
methylome_ewas_inputs <- function(id) {
  ns <- shiny::NS(id)
  bigdash::tabSettings(
    withTooltip(
      shiny::selectInput(
        ns("thresh_type"), "Threshold on:",
        choices = c("FDR (q-value)" = "q", "Nominal p-value" = "p"),
        selected = "q"
      ),
      "Call hits on the multiple-testing corrected q-value, or on the raw p-value.",
      placement = "top"
    ),
    withTooltip(
      shiny::numericInput(
        ns("thresh_value"), "Cut-off:",
        value = 0.05, min = 0, max = 1, step = 0.01
      ),
      "Probes at or below this value are called hits and listed in the table.",
      placement = "top"
    ),
    withTooltip(
      shiny::numericInput(
        ns("min_dbeta"), "Minimum |delta-beta|:",
        value = 0, min = 0, max = 1, step = 0.01
      ),
      paste(
        "Effect-size filter on the beta difference between the two groups.",
        "A large M-value change can still be a negligible change in methylation,",
        "so significance alone is a weak filter. 0 disables it."
      ),
      placement = "top"
    )
  )
}

## ------------------------------------------------------------------ shared --

## Note pgx$genes$chr is a cytoband ("16q12.2"), not a bare chromosome, so the
## arm and band are stripped exactly as board.epigenomics does.
mp_ewas_table <- function(pgx, contrast = NULL) {
  mp_require_methylomics(pgx)
  shiny::validate(shiny::need(!is.null(pgx$gx.meta),
    "This dataset has no differential methylation result."))
  cmp <- contrast
  if (is.null(cmp) || !cmp %in% names(pgx$gx.meta$meta)) cmp <- names(pgx$gx.meta$meta)[1]
  m <- pgx$gx.meta$meta[[cmp]]
  g <- pgx$genes[rownames(m), , drop = FALSE]
  cc <- grep("^chr$|chrom", tolower(colnames(g)))[1]
  pc <- grep("^pos$|position", tolower(colnames(g)))[1]
  isl <- mp_context_col(pgx, "island")
  loc <- mp_context_col(pgx, "gene")
  idx <- match(rownames(m), rownames(pgx$genes))

  data.frame(
    probe = rownames(m),
    gene = if ("symbol" %in% colnames(g)) as.character(g$symbol) else NA,
    p = as.numeric(m$meta.p), q = as.numeric(m$meta.q), fx = as.numeric(m$meta.fx),
    dbeta = mp_delta_beta(pgx, cmp)[rownames(m)],
    chr = if (!is.na(cc)) sub("(p|q|cen).*", "", sub("^chr", "", as.character(g[[cc]]))) else NA,
    pos = if (!is.na(pc)) suppressWarnings(as.numeric(sub(";.*", "", g[[pc]]))) else NA,
    island = if (!is.null(isl)) sub(";.*", "", isl[idx]) else NA,
    location = if (!is.null(loc)) sub(";.*", "", loc[idx]) else NA,
    stringsAsFactors = FALSE
  )
}

## The stored effect size is a logFC on M-values, which is not a quantity
## anyone reads. Recompute the difference in mean beta between the two groups
## of the contrast, which is what the field reports.
mp_delta_beta <- function(pgx, contrast) {
  X <- mp_beta(pgx)
  em <- pgx$model.parameters$exp.matrix
  if (is.null(em) || !contrast %in% colnames(em)) {
    return(stats::setNames(rep(NA_real_, nrow(X)), rownames(X)))
  }
  w <- em[, contrast]
  g1 <- names(w)[w > 0]
  g0 <- names(w)[w < 0]
  g1 <- intersect(g1, colnames(X))
  g0 <- intersect(g0, colnames(X))
  if (!length(g1) || !length(g0)) {
    return(stats::setNames(rep(NA_real_, nrow(X)), rownames(X)))
  }
  round(rowMeans(X[, g1, drop = FALSE], na.rm = TRUE) -
          rowMeans(X[, g0, drop = FALSE], na.rm = TRUE), 4)
}

## Which probes are hits under the current settings.
mp_ewas_sig <- function(d, thresh) {
  stat <- if (identical(thresh$type, "p")) d$p else d$q
  ok <- !is.na(stat) & stat <= thresh$value
  if (!is.null(thresh$min_dbeta) && thresh$min_dbeta > 0) {
    ok <- ok & !is.na(d$dbeta) & abs(d$dbeta) >= thresh$min_dbeta
  }
  ok
}

## Where to draw the line on a -log10(p) axis. A q-value cut-off has no fixed
## p, so use the largest p that still passes; with no hits there is no line.
mp_ewas_hline <- function(d, thresh) {
  if (identical(thresh$type, "p")) return(-log10(thresh$value))
  sig <- !is.na(d$q) & d$q <= thresh$value
  if (!any(sig)) return(NA_real_)
  -log10(max(d$p[sig], na.rm = TRUE))
}

mp_lambda <- function(p) {
  stats::median(stats::qchisq(1 - p, 1), na.rm = TRUE) / stats::qchisq(0.5, 1)
}

## ---------------------------------------------------------------- Manhattan --

methylome_plot_manhattan_ui <- function(id, title, caption, info.text, info.methods,
                                        height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "a"
  )
}

methylome_plot_manhattan_server <- function(id, pgx, r.contrast = shiny::reactive(NULL),
                                            r.thresh, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx())
      mp_ewas_table(pgx(), r.contrast())
    })
    plot.RENDER <- function() {
      d <- plot_data()
      th <- r.thresh()
      ## Carry significance as a column so it survives every subset and
      ## reorder below without a parallel vector to keep in step.
      d$sig <- mp_ewas_sig(d, th)
      d <- d[!is.na(d$chr) & !is.na(d$pos) & is.finite(d$pos) & is.finite(d$p), ]
      shiny::validate(shiny::need(nrow(d) > 0, "No probes carry chromosome and position."))
      chrs <- intersect(c(as.character(1:22), "X", "Y"), unique(d$chr))
      d <- d[d$chr %in% chrs, ]
      d$chr <- factor(d$chr, levels = chrs)
      d <- d[order(d$chr, d$pos), ]
      mx <- vapply(split(d$pos, d$chr), function(z) if (length(z)) max(z) else 0, numeric(1))
      off <- c(0, utils::head(cumsum(as.numeric(mx) * 1.02), -1))
      names(off) <- chrs
      d$gx <- d$pos + off[as.character(d$chr)]
      op <- graphics::par(mar = c(4, 4.5, 1, 1), las = 1); on.exit(graphics::par(op))
      cols <- ifelse(as.integer(d$chr) %% 2 == 0, "#a9b6c4", "#7f8f9f")
      plot(d$gx, -log10(d$p), pch = 19, cex = .35, col = cols, xaxt = "n",
           xlab = "chromosome", ylab = expression(-log[10](p)))
      graphics::points(d$gx[d$sig], -log10(d$p[d$sig]), pch = 19, cex = .6,
        col = ifelse(d$dbeta[d$sig] > 0 | is.na(d$dbeta[d$sig]),
                     MP_PAL$hyper, MP_PAL$hypo))
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
      csvFunc = plot_data, renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot, res = c(90, 130),
      pdf.width = 10, pdf.height = 5, add.watermark = watermark)
  })
}

## ----------------------------------------------------------------------- QQ --

methylome_plot_qq_ui <- function(id, title, caption, info.text, info.methods,
                                 height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "b"
  )
}

methylome_plot_qq_server <- function(id, pgx, r.contrast = shiny::reactive(NULL),
                                     watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx()); mp_ewas_table(pgx(), r.contrast())
    })
    plot.RENDER <- function() {
      d <- plot_data()
      pv <- sort(d$p[!is.na(d$p)])
      shiny::req(length(pv) > 0)
      op <- graphics::par(mar = c(4, 4.5, 1, 1), las = 1); on.exit(graphics::par(op))
      plot(-log10(stats::ppoints(length(pv))), -log10(pv), pch = 19, cex = .35,
           col = MP_PAL$grey, xlab = "expected", ylab = "observed")
      graphics::abline(0, 1, lty = 2, col = MP_PAL$hypo)
      graphics::legend("topleft", bty = "n",
                       legend = sprintf("lambda = %.3f", mp_lambda(d$p)))
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = plot_data, renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot, res = c(90, 130),
      pdf.width = 6, pdf.height = 6, add.watermark = watermark)
  })
}

## --------------------------------------------------------- hit context enrich --

methylome_plot_enrichment_ui <- function(id, title, caption, info.text, info.methods,
                                         height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "c"
  )
}

methylome_plot_enrichment_server <- function(id, pgx, r.contrast = shiny::reactive(NULL),
                                             r.thresh, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx())
      d <- mp_ewas_table(pgx(), r.contrast())
      d$sig <- mp_ewas_sig(d, r.thresh())
      d <- d[!is.na(d$island), ]
      sig <- d$sig
      shiny::validate(shiny::need(sum(sig) >= 5,
        "Too few hits at this threshold to test context enrichment."))
      lv <- intersect(c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"),
                      unique(d$island))
      ## Background is the probes actually tested in this contrast, not the
      ## whole array, or every category comes out enriched.
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
      ## Six horizontal categories in a short panel: trim the margin and the
      ## label size so none of them get clipped.
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

methylome_table_hits_server <- function(id, pgx, r.contrast = shiny::reactive(NULL),
                                        r.thresh, scrollY = "22vh") {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      shiny::req(pgx())
      d <- mp_ewas_table(pgx(), r.contrast())
      sig <- mp_ewas_sig(d, r.thresh())
      shiny::validate(shiny::need(any(sig),
        "No CpG passes the current threshold. Relax the cut-off in the settings panel."))
      d <- d[sig, , drop = FALSE]
      d <- d[order(d$p), , drop = FALSE]
      out <- data.frame(
        CpG = d$probe,
        Gene = ifelse(is.na(d$gene) | d$gene == "", "-", d$gene),
        Chr = d$chr,
        Position = d$pos,
        Context = ifelse(is.na(d$island), "-", d$island),
        Region = ifelse(is.na(d$location) | d$location == "", "-", d$location),
        `Delta beta` = d$dbeta,
        Direction = ifelse(is.na(d$dbeta), "-",
                           ifelse(d$dbeta > 0, "hyper", "hypo")),
        `P value` = signif(d$p, 3),
        FDR = signif(d$q, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      rownames(out) <- NULL
      out
    })

    render <- function(scrollY) {
      dt <- table_data()
      shiny::req(dt)
      DT::datatable(dt,
        class = "compact hover", rownames = FALSE,
        extensions = c("Buttons", "Scroller"), plugins = "scrollResize",
        selection = list(mode = "single", target = "row"),
        options = list(
          dom = "lfrtip", scroller = TRUE, scrollX = TRUE,
          scrollY = scrollY, scrollResize = TRUE, deferRender = TRUE,
          order = list(list(8, "asc"))
        )
      ) |>
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

    TableModuleServer(
      "tblmod",
      func = function() render(scrollY),
      func2 = function() render("60vh"),
      csvFunc = table_data,
      selector = "none"
    )
  })
}
