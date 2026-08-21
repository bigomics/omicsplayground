##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## "Did batch eat my experiment?" - asked before the EWAS rather than after it.
##
## The scatter beside this panel shows coordinates; coordinates alone are
## decoration, because the eye reads separation and cannot read what caused it.
## This tests every principal component against every sample column the model
## would accept as a covariate, so a slide, a plate or a sex effect that
## dominates a component is visible as a cell, not inferred from a cloud.
##
## It tests the SAME components the scatter draws - mp_pca_components() runs
## once per dataset and both panels read it, so the two cannot report a
## different PC1. That shared decomposition runs over every complete probe
## rather than the top 1,000 by SD precisely so this panel keeps working: a
## variance cut leaves cell type and sex untouched and drops the batch
## association from p 8.4e-09 to p 0.09 on the demo cohort.
mp_pc_covar <- function(pgx, pc) {
  mp_require_methylomics(pgx)
  pcs <- pc$pos
  varexp <- pc$varexp
  k <- ncol(pcs)
  n <- nrow(pcs)
  shiny::validate(shiny::need(n >= 6,
    "Too few samples to test: a principal component needs a cohort, not a handful."))

  ## The same columns the EWAS offers as covariates, and for the same reason:
  ## this panel answers "would adjusting for this change my answer". Ids,
  ## constants and free text are already gone from that list.
  vars <- mp_model_vars(pgx)
  smp <- rownames(pcs)
  types <- setNames(rep(NA_character_, length(vars)), vars)
  p <- r2 <- matrix(NA_real_, k, length(vars),
                    dimnames = list(colnames(pcs), vars))
  for (v in vars) {
    x <- pgx$samples[smp, v]
    num <- suppressWarnings(as.numeric(as.character(x)))
    ## A numeric column with few distinct values is a label wearing numbers -
    ## a slide id, a plate, a batch. Correlating a PC against it would report a
    ## trend through arbitrary codes; ten distinct values is where we stop
    ## believing the numbers mean magnitude.
    is_cont <- sum(!is.na(num)) >= 6 && length(unique(num[!is.na(num)])) > 10
    if (is_cont) {
      types[v] <- "continuous"
      ok <- !is.na(num)
      if (sum(ok) < 6 || stats::sd(num[ok]) == 0) next
      for (i in seq_len(k)) {
        ct <- tryCatch(stats::cor.test(pcs[ok, i], num[ok]), error = function(e) NULL)
        if (is.null(ct)) next
        p[i, v] <- ct$p.value
        r2[i, v] <- unname(ct$estimate)^2
      }
    } else {
      types[v] <- "categorical"
      g <- as.character(x)
      ok <- !is.na(g) & g != ""
      g <- droplevels(factor(g[ok]))
      if (nlevels(g) < 2 || nlevels(g) > 12) next
      for (i in seq_len(k)) {
        y <- pcs[ok, i]
        kt <- tryCatch(stats::kruskal.test(y ~ g), error = function(e) NULL)
        if (is.null(kt)) next
        p[i, v] <- kt$p.value
        ## Eta-squared: the share of that PC's spread the grouping accounts
        ## for. Kruskal says "not by chance"; this says "and by how much",
        ## which is the number that decides whether to adjust.
        r2[i, v] <- summary(stats::lm(y ~ g))$r.squared
      }
    }
  }

  ## One BH across the whole grid: every cell in this panel is looked at, so
  ## correcting per column would under-count the tests actually performed.
  q <- p
  q[] <- stats::p.adjust(p, method = "BH")
  list(pc = pcs, varexp = varexp, p = p, q = q, r2 = r2,
       types = types, n = n, n_probes = pc$n_probes, thinned = pc$thinned)
}

methylome_plot_pccovar_ui <- function(id, title, caption, info.text, info.methods,
                                      height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    translate = FALSE,
    ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "b"
  )
}

methylome_plot_pccovar_server <- function(id, pgx, pca_components, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    ## One decomposition per dataset, reused by the redraw and the CSV.
    r.pccov <- shiny::reactive({
      p <- pgx()
      shiny::req(p)
      mp_pc_covar(p, pca_components())
    })

    plot.RENDER <- function() {
      r <- r.pccov()
      shiny::validate(shiny::need(ncol(r$p) > 0,
        "No sample column can be tested against the components: every column is an identifier, a constant, or free text with too many values."))
      k <- nrow(r$p); nv <- ncol(r$p)
      z <- pmin(-log10(r$p), 16)
      z[is.na(z)] <- 0
      ## PC1 on top: image() draws y upward, so the rows go in reversed.
      ord <- rev(seq_len(k))
      ramp <- grDevices::colorRampPalette(c("#f7f7f7", "#f4c9c4", MP_PAL$hyper))(64)
      op <- graphics::par(mar = c(9, 6.4, 1.2, 1.2), las = 1, xpd = NA)
      on.exit(graphics::par(op))
      graphics::image(x = seq_len(nv), y = seq_len(k),
                      z = t(z[ord, , drop = FALSE]), col = ramp, zlim = c(0, 16),
                      xlab = "", ylab = "", axes = FALSE)
      graphics::axis(1, at = seq_len(nv), labels = colnames(r$p),
                     las = 2, cex.axis = 0.85, tick = FALSE)
      graphics::axis(2, at = seq_len(k),
                     labels = sprintf("PC%d (%.1f%%)", ord, r$varexp[ord]),
                     cex.axis = 0.85, tick = FALSE)
      graphics::box(col = "#c8d2dc")
      for (i in seq_len(k)) {
        for (j in seq_len(nv)) {
          if (is.na(r$r2[i, j])) next
          graphics::text(j, which(ord == i), sprintf("%.2f", r$r2[i, j]),
                         cex = 0.78, font = if (r$q[i, j] < 0.05) 2 else 1,
                         col = if (z[i, j] > 9) "white" else "#33414f")
        }
      }
      ## Two short lines, not one long one: this panel sits in five of twelve
      ## columns and a single-line key ran off the right edge.
      graphics::mtext("cell = variance of the component explained",
                      side = 1, line = 6.9, cex = 0.78, col = MP_PAL$grey, adj = 0)
      graphics::mtext("shading = -log10 p   |   bold = FDR < 5%",
                      side = 1, line = 7.9, cex = 0.78, col = MP_PAL$grey, adj = 0)
    }

    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive({
        r <- r.pccov()
        data.frame(PC = rep(rownames(r$p), ncol(r$p)),
                   variable = rep(colnames(r$p), each = nrow(r$p)),
                   type = rep(r$types[colnames(r$p)], each = nrow(r$p)),
                   varexp = round(rep(r$varexp, ncol(r$p)), 2),
                   r2 = round(as.vector(r$r2), 4),
                   p = signif(as.vector(r$p), 3),
                   q = signif(as.vector(r$q), 3),
                   stringsAsFactors = FALSE)
      }),
      renderFunc = shiny::renderPlot, renderFunc2 = shiny::renderPlot,
      res = c(90, 130), pdf.width = 9, pdf.height = 6, add.watermark = watermark)
  })
}
