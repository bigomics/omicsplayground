##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## EWAS. The only screen here that reads pgx$gx.meta - the differential
## methylation already computed at upload, which no methylation board has
## displayed before.

## Shared by the three plots. Note pgx$genes$chr is a cytoband ("16q12.2"),
## not a bare chromosome, so the arm and band are stripped exactly as
## board.epigenomics does.
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
  data.frame(
    probe = rownames(m),
    p = as.numeric(m$meta.p), q = as.numeric(m$meta.q), fx = as.numeric(m$meta.fx),
    chr = if (!is.na(cc)) sub("(p|q|cen).*", "", sub("^chr", "", as.character(g[[cc]]))) else NA,
    pos = if (!is.na(pc)) suppressWarnings(as.numeric(sub(";.*", "", g[[pc]]))) else NA,
    island = if (!is.null(isl)) sub(";.*", "", isl[match(rownames(m), rownames(pgx$genes))]) else NA,
    stringsAsFactors = FALSE
  )
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
                                            watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx())
      mp_ewas_table(pgx(), r.contrast())
    })
    plot.RENDER <- function() {
      d <- plot_data()
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
      sig <- d$q < 0.05 & !is.na(d$q)
      graphics::points(d$gx[sig], -log10(d$p[sig]), pch = 19, cex = .6,
        col = ifelse(d$fx[sig] > 0, MP_PAL$hyper, MP_PAL$hypo))
      graphics::axis(1, at = tapply(d$gx, d$chr, stats::median, na.rm = TRUE),
                     labels = levels(d$chr), tick = FALSE, cex.axis = .7)
      graphics::abline(h = -log10(1e-7), lty = 2, col = MP_PAL$grey)
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
                                             watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx())
      d <- mp_ewas_table(pgx(), r.contrast())
      d <- d[!is.na(d$island), ]
      sig <- d$q < 0.05 & !is.na(d$q)
      shiny::validate(shiny::need(sum(sig) >= 5,
        "Too few significant CpGs to test context enrichment."))
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
      op <- graphics::par(mar = c(4, 8, 1, 2), las = 1); on.exit(graphics::par(op))
      cols <- ifelse(e$log2_or >= 0, MP_PAL$hyper, MP_PAL$hypo)
      graphics::barplot(rev(e$log2_or), horiz = TRUE, col = rev(cols), border = NA,
        names.arg = rev(e$context), xlab = "log2 odds ratio vs tested background")
      graphics::abline(v = 0, col = "#444444")
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = plot_data, renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot, res = c(90, 130),
      pdf.width = 6, pdf.height = 5, add.watermark = watermark)
  })
}
