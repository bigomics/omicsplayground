##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Regions and pathways: a real DMR caller, and gene-set testing corrected for
## the probes-per-gene bias.
##
## dmrff rather than DMRcate: dmrff runs on EWAS summary statistics, which is
## exactly what the model step already produces, and its entire dependency
## list is parallel + ggplot2. DMRcate would drag in bsseq, Gviz, biomaRt,
## missMethyl and an ExperimentHub package that downloads at runtime.
##
## The dependency argument is the whole argument. An earlier version of this
## comment claimed "a 2024 BMC Bioinformatics comparison put dmrff at the
## highest power with a controlled false-positive rate" - there is no such
## paper, and the published DMR benchmarks (Mallik 2019, Brief Bioinform
## 20:2224; Zheng 2022, Epigenetics 17:2241) do not include dmrff at all.
## dmrff itself has never been journal-published either: it is still a preprint
## (Suderman 2018, bioRxiv 508556). Worth knowing before quoting it in a paper.

## ------------------------------------------------------------------- DMRs --

mp_can_dmrff <- function() requireNamespace("dmrff", quietly = TRUE)

#' Call DMRs from the fitted model.
#'
#' dmrff needs the effect estimate, its standard error and the p-value per
#' probe, plus chromosome and position. limma does not return the standard
#' error directly, but se = logFC / t and t = logFC / se, so it comes back
#' out of the moderated t statistic.
mp_call_dmrs <- function(pgx, res, maxgap = 500) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  shiny::validate(shiny::need(mp_can_dmrff(),
    "The dmrff package is not installed, so region calling is unavailable."))
  d <- res$data
  fit <- res$meta
  ## dmrff combines a signed estimate and its standard error across neighbouring
  ## probes. An F-test across three or more groups produces neither, and any
  ## stand-in for them would be an invented direction, so refuse.
  shiny::validate(shiny::need(!isTRUE(fit$anova),
    "Region calling needs a signed effect and its standard error per probe. This model is an F-test across three or more groups, which has neither. Fit a two-group contrast or a continuous outcome to call regions."))
  ## dmrff joins neighbouring probes by genomic position. A collapsed fit has
  ## already merged them, and a feature's position is the median of the probes
  ## that went into it - joining those would be calling regions of regions.
  shiny::validate(shiny::need(identical(fit$collapse %||% "probe", "probe"),
    "Region calling works on per-probe statistics, and this model was fitted on collapsed features. Set 'Test at' back to probe level to call regions."))
  shiny::validate(shiny::need(!is.null(fit$t),
    "This model did not return moderated t statistics."))

  keep <- !is.na(d$chr) & !is.na(d$pos) & is.finite(d$pos) &
    is.finite(d$p) & is.finite(fit$t) & fit$t != 0
  d <- d[keep, , drop = FALSE]
  tt <- fit$t[keep]
  shiny::validate(shiny::need(nrow(d) > 100, "Too few positioned probes to call regions."))

  est <- d$logFC_M
  se <- abs(est / tt)
  se[!is.finite(se) | se <= 0] <- NA

  beta <- mp_beta(pgx)
  ## The samples the model was actually fitted on, not every sample carrying a
  ## contrast label: incomplete covariates drop rows, and dmrff must see the
  ## same set the summary statistics came from.
  ss <- intersect(res$meta$samples, colnames(beta))
  M <- playbase.epigenetics::betaToM(beta[d$probe, ss, drop = FALSE])

  ok <- !is.na(se)
  out <- tryCatch(
    dmrff::dmrff(estimate = est[ok], se = se[ok], p.value = d$p[ok],
                 methylation = M[ok, , drop = FALSE],
                 chr = paste0("chr", d$chr[ok]), pos = d$pos[ok],
                 maxgap = maxgap, verbose = FALSE),
    error = function(e) {
      warning("[methylome] dmrff failed: ", conditionMessage(e)); NULL
    })
  shiny::validate(shiny::need(!is.null(out), "Region calling failed on this dataset."))
  ## Single-CpG "regions" are just DMPs and are already in the hits table.
  out <- out[!is.na(out$n) & out$n > 1, , drop = FALSE]
  out
}

## Genes overlapping each region, from the annotation already loaded.
mp_dmr_genes <- function(pgx, dmrs, d) {
  if (!nrow(dmrs)) return(character(0))
  vapply(seq_len(nrow(dmrs)), function(i) {
    hit <- d$chr == sub("^chr", "", dmrs$chr[i]) &
      d$pos >= dmrs$start[i] & d$pos <= dmrs$end[i]
    g <- unique(d$gene[hit])
    g <- g[!is.na(g) & g != ""]
    if (!length(g)) "-" else paste(utils::head(g, 4), collapse = ", ")
  }, character(1))
}

methylome_table_dmr_ui <- function(id, title, info.text, caption, height, width) {
  ns <- shiny::NS(id)
  TableModuleUI(
    translate = FALSE,ns("tblmod"), title = title, info.text = info.text,
                caption = caption, height = height, width = width, label = "a")
}

methylome_table_dmr_server <- function(id, pgx, r.regions, scrollY = "26vh") {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      dd <- r.regions()
      shiny::validate(shiny::need(!is.null(dd) && nrow(dd$dmrs) > 0,
        "No regions called yet. Fit a model on the Manhattan & hits sub-tab first, then press Call regions in the settings panel."))
      x <- dd$dmrs
      data.frame(
        Region = sprintf("%s:%s-%s", x$chr, format(x$start, big.mark = ","),
                         format(x$end, big.mark = ",")),
        Genes = dd$genes,
        `CpGs` = x$n,
        `Width (bp)` = x$end - x$start,
        `Effect (M)` = signif(x$estimate, 3),
        `P value` = signif(x$p.value, 3),
        `Adjusted p` = signif(x$p.adjust, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })
    render <- function(sy) {
      dt <- table_data(); shiny::req(dt)
      DT::datatable(dt, class = "compact hover", rownames = FALSE,
        selection = list(mode = "single", target = "row"),
        extensions = c("Buttons", "Scroller"), plugins = "scrollResize",
        options = list(dom = "lfrtip", scroller = TRUE, scrollX = TRUE,
                       scrollY = sy, scrollResize = TRUE, deferRender = TRUE,
                       order = list(list(5, "asc")))) |>
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }
    ## Returned so the caller can drive the region plot from the clicked row -
    ## rows_selected indexes the data as passed to DT, which is the same order
    ## as dmrs, so it needs no lookup however the table is sorted on screen.
    TableModuleServer("tblmod", func = function() render(scrollY),
                      func2 = function() render("60vh"),
                      csvFunc = table_data, selector = "single")
  })
}

## ---------------------------------------------------------- region detail --
##
## A DMR as a table row is not interpretable: a run of ten CpGs moving together
## and one strong CpG dragging its neighbours produce similar statistics. This
## is the plot the field draws instead (DMRcate::DMR.plot, coMET) - methylation
## against genomic position across the region, so coherence is visible.

## Probes inside a region, plus flanking context so the region has an edge.
mp_region_probes <- function(res, dmr, pad = 0.3) {
  d <- res$data
  w <- max(500, (dmr$end - dmr$start) * pad)
  k <- !is.na(d$chr) & d$chr == sub("^chr", "", dmr$chr) & !is.na(d$pos) &
    d$pos >= dmr$start - w & d$pos <= dmr$end + w
  d <- d[k, , drop = FALSE]
  d[order(d$pos), , drop = FALSE]
}

## Island context shading for the track under the axis.
MP_ISLAND_COL <- c(Island = "#2c6b3f", N_Shore = "#7bab84", S_Shore = "#7bab84",
                   N_Shelf = "#c3ddc7", S_Shelf = "#c3ddc7", OpenSea = "#eef1ee")

methylome_plot_dmrregion_ui <- function(id, title, caption, info.text, info.methods,
                                        height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    translate = FALSE,ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "c")
}

methylome_plot_dmrregion_server <- function(id, pgx, r.ewas, r.regions, r.pick,
                                            watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      rg <- r.regions()
      shiny::validate(shiny::need(!is.null(rg) && nrow(rg$dmrs) > 0,
        "No regions to draw. Fit a model on the Manhattan & hits sub-tab, then press Call regions in the settings panel."))
      res <- r.ewas()
      ## No row clicked yet: show the strongest region rather than whichever
      ## one dmrff happened to return first.
      i <- suppressWarnings(as.integer(r.pick()))
      if (!length(i) || is.na(i)) i <- which.min(rg$dmrs$p.value)
      i <- max(1L, min(i, nrow(rg$dmrs)))
      dm <- rg$dmrs[i, ]
      d <- mp_region_probes(res, dm)
      shiny::validate(shiny::need(nrow(d) >= 2,
        "Too few positioned probes around this region to draw it."))
      X <- mp_beta(pgx())
      ss <- intersect(res$meta$samples, colnames(X))
      list(dm = dm, gene = rg$genes[i], d = d,
           B = X[d$probe, ss, drop = FALSE],
           grp = droplevels(res$meta$strata[ss]),
           continuous = isTRUE(res$meta$continuous))
    })

    plot.RENDER <- function() {
      z <- plot_data(); d <- z$d; B <- z$B; grp <- z$grp; dm <- z$dm
      lv <- levels(grp)
      cols <- grDevices::colorRampPalette(c(MP_PAL$hypo, MP_PAL$hyper))(max(2, length(lv)))
      op <- graphics::par(mar = c(4.2, 4.4, 2.4, 1), las = 1, mgp = c(2.6, 0.6, 0))
      on.exit(graphics::par(op))

      xr <- range(d$pos)
      ## Room under zero for the island track.
      graphics::plot(NA, xlim = xr, ylim = c(-0.1, 1), yaxt = "n",
                     xlab = sprintf("chr%s position (bp)", sub("^chr", "", dm$chr)),
                     ylab = "beta")
      graphics::axis(2, at = seq(0, 1, 0.25))
      ## The called region, so flanking context is visibly outside it.
      graphics::rect(dm$start, 0, dm$end, 1, col = "#f2f5f8", border = NA)
      graphics::abline(h = c(0.2, 0.8), lty = 3, col = "#c9ccc7")

      ## Every sample, faint, then the stratum means over the top.
      for (j in seq_along(lv)) {
        ss <- names(grp)[grp == lv[j]]
        for (s in ss) {
          graphics::lines(d$pos, B[, s], col = grDevices::adjustcolor(cols[j], 0.12))
        }
      }
      for (j in seq_along(lv)) {
        mu <- rowMeans(B[, names(grp)[grp == lv[j]], drop = FALSE], na.rm = TRUE)
        graphics::lines(d$pos, mu, col = cols[j], lwd = 2.6)
        graphics::points(d$pos, mu, col = cols[j], pch = 19, cex = 0.8)
      }

      ## Island context track.
      isl <- ifelse(is.na(d$island), "OpenSea", d$island)
      tc <- ifelse(isl %in% names(MP_ISLAND_COL), MP_ISLAND_COL[isl], "#eef1ee")
      graphics::rect(d$pos - diff(xr) * 0.004, -0.075,
                     d$pos + diff(xr) * 0.004, -0.03, col = tc, border = NA)
      graphics::mtext("island context", side = 1, line = -1.1, adj = 0,
                      cex = 0.65, col = MP_PAL$grey)

      graphics::title(main = sprintf("chr%s:%s-%s   %s   %d CpGs",
                                     sub("^chr", "", dm$chr),
                                     format(dm$start, big.mark = ","),
                                     format(dm$end, big.mark = ","),
                                     z$gene, dm$n),
                      cex.main = 1, font.main = 1)
      graphics::legend("topright", bty = "n", horiz = TRUE, cex = 0.8,
                       legend = lv, col = cols[seq_along(lv)], lwd = 2.6,
                       title = if (z$continuous) "outcome tertile" else NULL)
    }

    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive({
        z <- plot_data()
        cbind(probe = z$d$probe, pos = z$d$pos, island = z$d$island,
              as.data.frame(z$B))
      }),
      renderFunc = shiny::renderPlot, renderFunc2 = shiny::renderPlot,
      res = c(90, 130), pdf.width = 9, pdf.height = 5.5, add.watermark = watermark)
  })
}

## ------------------------------------------------------------- enrichment --

mp_can_gometh <- function() requireNamespace("missMethyl", quietly = TRUE)

#' Gene-set testing on the hit list, corrected for probes per gene.
#'
#' Genes carry very different numbers of probes, so a naive hypergeometric
#' test returns the same probe-dense developmental gene sets whatever the
#' biology (Geeleher 2013). gometh reweights with Wallenius' noncentral
#' hypergeometric and also handles CpGs mapping to several genes.
mp_run_gometh <- function(sig_cpg, all_cpg, collection = "GO", array_type = "450K") {
  shiny::validate(shiny::need(mp_can_gometh(),
    "The missMethyl package is not installed, so bias-corrected enrichment is unavailable."))
  shiny::validate(shiny::need(length(sig_cpg) >= 10,
    "Too few hits at this threshold for a meaningful gene-set test."))
  out <- tryCatch(
    missMethyl::gometh(sig.cpg = sig_cpg, all.cpg = all_cpg,
                       collection = collection, array.type = array_type,
                       prior.prob = TRUE),
    error = function(e) {
      warning("[methylome] gometh failed: ", conditionMessage(e)); NULL
    })
  shiny::validate(shiny::need(!is.null(out), "Gene-set testing failed on this dataset."))
  out$TERM_ID <- rownames(out)
  out[order(out$P.DE), , drop = FALSE]
}

methylome_table_enrich_ui <- function(id, title, info.text, caption, height, width) {
  ns <- shiny::NS(id)
  TableModuleUI(
    translate = FALSE,ns("tblmod"), title = title, info.text = info.text,
                caption = caption, height = height, width = width, label = "b")
}

methylome_table_enrich_server <- function(id, r.enrich, scrollY = "26vh") {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      e <- r.enrich()
      shiny::validate(shiny::need(!is.null(e) && nrow(e) > 0,
        "No gene-set result yet. Fit a model on the Manhattan & hits sub-tab, then press Test gene sets in the settings panel."))
      keep <- utils::head(e, 300)
      data.frame(
        Term = if ("TERM" %in% colnames(keep)) keep$TERM else keep$TERM_ID,
        Ontology = if ("ONTOLOGY" %in% colnames(keep)) keep$ONTOLOGY else "-",
        `Genes in set` = keep$N,
        `Hits` = keep$DE,
        `P value` = signif(keep$P.DE, 3),
        FDR = signif(keep$FDR, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })
    render <- function(sy) {
      dt <- table_data(); shiny::req(dt)
      DT::datatable(dt, class = "compact hover", rownames = FALSE, selection = "none",
        extensions = c("Buttons", "Scroller"), plugins = "scrollResize",
        options = list(dom = "lfrtip", scroller = TRUE, scrollX = TRUE,
                       scrollY = sy, scrollResize = TRUE, deferRender = TRUE)) |>
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }
    TableModuleServer("tblmod", func = function() render(scrollY),
                      func2 = function() render("60vh"),
                      csvFunc = table_data, selector = "none")
  })
}
