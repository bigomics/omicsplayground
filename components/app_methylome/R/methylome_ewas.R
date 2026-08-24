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
## `annot` overrides pgx$genes: a collapsed fit tests features that are not
## probes, so its annotation was built by mp_collapse_beta() and pgx$genes has
## no row for them.
mp_ewas_annotate <- function(pgx, tab, annot = NULL) {
  src <- if (is.null(annot)) pgx$genes else annot
  g <- src[tab$probe, , drop = FALSE]
  ## Every field comes off the same frame the features were looked up in.
  ## Island and context used to be read from pgx$genes by row index, which is
  ## NA for every feature of a collapsed fit - those rows are not probes.
  col <- function(pat) {
    hit <- grep(pat, colnames(src), ignore.case = TRUE)
    if (!length(hit)) NULL else as.character(g[[colnames(src)[hit[1]]]])
  }
  cc <- col("^chr$|chrom")
  pc <- col("^pos$|position")
  tab$gene <- if ("symbol" %in% colnames(g)) as.character(g$symbol) else NA
  tab$chr  <- if (!is.null(cc)) sub("(p|q|cen).*", "", sub("^chr", "", cc)) else NA
  tab$pos  <- if (!is.null(pc)) suppressWarnings(as.numeric(sub(";.*", "", pc))) else NA
  isl <- col("Relation_to_Island")
  loc <- col("genomic_location|RefGene_Group")
  tab$island   <- if (!is.null(isl)) sub(";.*", "", isl) else NA
  tab$location <- if (!is.null(loc)) sub(";.*", "", loc) else NA
  ## Only a collapsed fit has this. It has to reach the table: probes per gene
  ## region run from 3 to over a thousand, and a feature built from 400 probes
  ## has a much smaller residual variance than one built from 3 - so the reader
  ## needs to see which kind of feature a hit is.
  if ("n_probes" %in% colnames(g)) tab$n_probes <- as.integer(g$n_probes)
  ## Collapsed fits carry the extra fields the table needs to be honest about
  ## what an aggregated row is: the span rather than a midpoint, the partition
  ## the feature actually belongs to, and how many probes back the context.
  for (k in c("region", "pos_start", "pos_end", "island_n")) {
    if (k %in% colnames(g)) tab[[k]] <- g[[k]]
  }
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
  ## Drawn from the same rule that colours the points and counts the hits -
  ## using the q cut-off alone put the line somewhere no plotted hit sat
  ## whenever the effect-size filter was on.
  sig <- mp_ewas_sig(d, thresh)
  if (!any(sig)) {
    return(if (identical(thresh$type, "p")) -log10(thresh$value) else NA_real_)
  }
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
    out <- list(bias = as.numeric(bacon::bias(bc)),
                inflation = as.numeric(bacon::inflation(bc)))
    ## bacon estimates an empirical null by fitting a three-component mixture,
    ## which assumes most tests are null. A contrast where they are not - a sex
    ## EWAS on a dataset carrying the sex chromosomes puts 36% of probes past
    ## q 0.05 and takes t to -142 - leaves the EM without a null to find, and it
    ## returns NaN rather than failing. Report nothing instead of printing NaN
    ## at the reader: lambda is beside it and is still meaningful.
    if (!all(vapply(out, is.finite, logical(1)))) return(NULL)
    out
  }, error = function(e) NULL)
}

## ---------------------------------------------------------------- Manhattan --

methylome_plot_manhattan_ui <- function(id, title, caption, info.text, info.methods,
                                        height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    translate = FALSE,ns("pltmod"), plotlib = "base", title = title, caption = caption,
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
      ## Count before the plot drops what it cannot place. A feature with no
      ## chromosome or position is still a hit; it just cannot be drawn. The
      ## corner annotation reads as THE hit count, and counting after the drop
      ## put it 31% below the hits table on the same screen - 6,534 against
      ## 9,422 - with no hint that the two meant different things.
      n_hits <- sum(d$sig)
      n_plotted <- sum(d$sig & !is.na(d$chr) & !is.na(d$pos))
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
      ## An F-test across three or more groups has no direction, so hits get one
      ## accent colour rather than a hyper/hypo split that would be invented.
      hitcol <- if (isTRUE(r.ewas()$meta$anova)) MP_PAL$bad else
        ifelse(d$dbeta[d$sig] > 0 | is.na(d$dbeta[d$sig]), MP_PAL$hyper, MP_PAL$hypo)
      graphics::points(d$gx[d$sig], -log10(d$p[d$sig]), pch = 19, cex = .6, col = hitcol)
      graphics::axis(1, at = tapply(d$gx, d$chr, stats::median, na.rm = TRUE),
                     labels = levels(d$chr), tick = FALSE, cex.axis = .7)
      h <- mp_ewas_hline(d, th)
      if (is.finite(h)) {
        graphics::abline(h = h, lty = 2, col = MP_PAL$grey)
        graphics::mtext(sprintf("%s <= %g  (%s hits%s)", th$type, th$value,
                                format(n_hits, big.mark = ","),
                                if (n_plotted < n_hits)
                                  sprintf("; %s drawn, the rest carry no position",
                                          format(n_plotted, big.mark = ",")) else ""),
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
  PlotModuleUI(
    translate = FALSE,ns("pltmod"), plotlib = "base", title = title, caption = caption,
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
  PlotModuleUI(
    translate = FALSE,ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "c")
}

methylome_plot_enrichment_server <- function(id, r.ewas, r.thresh, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      res <- r.ewas()
      ## Refuses on a collapsed fit. The island relation of an aggregated
      ## feature is a majority vote over its probes, and a gene body routinely
      ## spans island, shore, shelf and open sea - on this array 12% of features
      ## carry a label under 60% of their own probes agree with. An odds ratio
      ## computed over those labels measures the vote, not where the hits sit.
      shiny::validate(shiny::need(identical(res$meta$collapse, "probe"),
        "Genomic context is a per-probe property, and this model was fitted on gene regions whose context is a majority vote over their probes. Set 'Test at' to probe level for this panel."))
      d <- res$data
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

## ---------------------------------------------------------- EWAS Catalog ---

## Deep links back into ewascatalog.org. Both query forms were checked against
## the live site: ?cpg= returns that CpG's association table, ?trait= returns
## every CpG reported for the trait. Following a link is the browser's job, so
## this adds no egress requirement to the container - which is why bundling the
## data and linking out are not in tension.
MP_EWAS_SITE <- "https://www.ewascatalog.org/"

## The arrow icon is the link, the text beside it stays plain - same as the
## signature board does with geneset accessions. Making the label itself the
## anchor turns every trait name into something that navigates away when you
## meant to select the row.
##
## Trait strings come out of a downloaded file and land in a cell rendered with
## escape = FALSE, so the value is URL-encoded here and any label the caller
## prints beside it must be HTML-escaped there. Neither is optional: without
## them a trait containing a quote or a tag would be injected into the page.
mp_ewas_link <- function(param, value) {
  url <- paste0(MP_EWAS_SITE, "?", param, "=",
                vapply(value, utils::URLencode, character(1), reserved = TRUE, USE.NAMES = FALSE))
  paste0("<a href='", url, "' target='_blank' rel='noopener' title='Open in the EWAS Catalog'>",
         "<i class='fa-solid fa-arrow-up-right-from-square weblink'></i></a>")
}

## The trait list as the table shows it, each trait a link to its catalog page
## and the trailing "+N more" a link to the CpG's own page - which is where the
## rest of them are.
mp_catalog_reported_html <- function(reported, probe) {
  vapply(seq_along(reported), function(i) {
    s <- reported[i]
    if (is.na(s) || !nzchar(s)) {
      ## Inline rather than a class: styles.min.css is a compiled artefact and
      ## one grey span is not worth a sass rebuild in the loop.
      return("<span style='color:#9aa5b1'>not reported</span>")
    }
    parts <- strsplit(s, "; ", fixed = TRUE)[[1]]
    more <- grepl("^\\+[0-9]+ more$", parts)
    ## "Smoking (46)" -> plain text, then the icon that opens its trait page.
    linked <- vapply(parts[!more], function(p) {
      tr <- sub(" \\([0-9]+\\)$", "", p)
      paste0(htmltools::htmlEscape(p), mp_ewas_link("trait", tr))
    }, character(1), USE.NAMES = FALSE)
    if (any(more)) {
      ## The rest of the traits are on the CpG's own page, so that is where the
      ## remainder points.
      linked <- c(linked, paste0(htmltools::htmlEscape(parts[more][1]),
                                 mp_ewas_link("cpg", probe[i])))
    }
    paste(linked, collapse = "; ")
  }, character(1), USE.NAMES = FALSE)
}

## Per-probe summary of the catalog slice: the formatted trait list the table
## shows, plus the two counts worth sorting on. Kept here rather than in
## playbase.epigenetics because the shape is what this table renders, not a
## property of the catalog.
##
## `studies` counts distinct study analyses across all of a CpG's traits, so it
## is "how hard has anyone looked at this CpG", which is the number you sort
## ascending to put your novel hits on top.
mp_catalog_summary <- function(traits, probes, max_traits = 3) {
  empty <- data.frame(probe = probes, reported = "", studies = 0L, n_traits = 0L,
                      stringsAsFactors = FALSE)
  if (is.null(traits) || !nrow(traits)) return(empty)
  d <- traits[order(as.integer(traits$cpg), -traits$n), , drop = FALSE]
  key <- as.character(d$cpg)
  r <- rle(key)
  top <- sequence(r$lengths) <= max_traits
  txt <- vapply(split(sprintf("%s (%d)", as.character(d$trait[top]), d$n[top]),
                      factor(key[top], levels = r$values)),
                paste, character(1), collapse = "; ")
  more <- r$lengths - max_traits
  txt <- ifelse(more > 0, paste0(txt, "; +", more, " more"), txt)
  i <- match(probes, r$values)
  data.frame(
    probe = probes,
    reported = ifelse(is.na(i), "", txt[i]),
    studies = ifelse(is.na(i), 0L, as.integer(vapply(split(d$n, factor(key, levels = r$values)),
                                                     sum, numeric(1))[i])),
    n_traits = ifelse(is.na(i), 0L, r$lengths[i]),
    stringsAsFactors = FALSE
  )
}

## ------------------------------------------------- promoter vs body ---

## The one question only the collapsed fit can answer. Promoter methylation
## silences while gene-body methylation correlates positively with expression
## (Yang 2014, Cancer Cell 26:577; Jones 2012, Nat Rev Genet 13:484), so a gene
## whose two halves move in opposite directions is a specific finding - and it
## is exactly the pattern that makes averaging a whole gene meaningless. On a
## full 450K array about two thirds of genes carry both features, so this is a
## view of most of the hit list, not a corner of it.
mp_divergence <- function(res, thresh) {
  m <- res$meta
  shiny::validate(shiny::need(identical(m$collapse, "gene_region"),
    "This panel compares each gene's promoter against its body, so it needs a fit collapsed to gene regions. Set 'Test at' to gene region."))
  ## An F-test reports an unsigned beta range, so "opposite directions" has no
  ## meaning on that branch - refuse rather than compare magnitudes and call it
  ## a direction.
  shiny::validate(shiny::need(!isTRUE(m$anova),
    "This model is an F-test across three or more groups, which has no direction, so promoter and body cannot be compared for opposing change. Fit a two-group contrast or a continuous outcome."))
  d <- res$data
  d$sym <- sub(":(promoter|body)$", "", d$probe)
  d$reg <- sub("^.*:", "", d$probe)
  pr <- d[d$reg == "promoter", , drop = FALSE]
  bo <- d[d$reg == "body", , drop = FALSE]
  g <- intersect(pr$sym, bo$sym)
  shiny::validate(shiny::need(length(g) > 0,
    "No gene in this dataset has both a promoter and a gene-body feature that survived the probe floor."))
  i <- match(g, pr$sym); j <- match(g, bo$sym)
  out <- data.frame(
    gene = g,
    dbeta_promoter = pr$dbeta[i], dbeta_body = bo$dbeta[j],
    q_promoter = pr$q[i], q_body = bo$q[j],
    p_promoter = pr$p[i], p_body = bo$p[j],
    n_promoter = if (is.null(pr$n_probes)) NA_integer_ else pr$n_probes[i],
    n_body = if (is.null(bo$n_probes)) NA_integer_ else bo$n_probes[j],
    stringsAsFactors = FALSE
  )
  sig_p <- mp_ewas_sig(data.frame(p = out$p_promoter, q = out$q_promoter,
                                  dbeta = out$dbeta_promoter), thresh)
  sig_b <- mp_ewas_sig(data.frame(p = out$p_body, q = out$q_body,
                                  dbeta = out$dbeta_body), thresh)
  out$both_sig <- sig_p & sig_b
  ## Sign, not magnitude: the finding is that the two halves disagree.
  out$opposite <- !is.na(out$dbeta_promoter) & !is.na(out$dbeta_body) &
    sign(out$dbeta_promoter) != sign(out$dbeta_body)
  ## Body minus promoter, which is the published direction: Li 2019 defines
  ## MeGDP as "the numerical DNA methylation difference between gene body and
  ## promoter regions" (Genome Res 29:270). Computing it the other way round
  ## would have been the same number wearing a different sign and a made-up
  ## name, so it follows theirs.
  out$megdp <- out$dbeta_body - out$dbeta_promoter
  out[order(-out$both_sig, -out$opposite, -abs(out$megdp)), , drop = FALSE]
}

methylome_plot_divergence_ui <- function(id, title, caption, info.text, info.methods,
                                         height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    translate = FALSE,ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "a")
}

methylome_plot_divergence_server <- function(id, r.ewas, r.thresh, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive(mp_divergence(r.ewas(), r.thresh()))
    plot.RENDER <- function() {
      e <- plot_data(); shiny::req(nrow(e) > 0)
      op <- graphics::par(mar = c(4.2, 4.4, 2.2, 1), las = 1, mgp = c(2.4, 0.6, 0))
      on.exit(graphics::par(op))
      lim <- max(abs(c(e$dbeta_promoter, e$dbeta_body)), na.rm = TRUE) * 1.05
      ## Grey is everything, blue is significant in both, red is significant in
      ## both AND disagreeing - the only class that is a finding.
      col <- rep("#d5dbe1", nrow(e))
      col[e$both_sig] <- MP_PAL$hypo
      col[e$both_sig & e$opposite] <- MP_PAL$hyper
      cex <- ifelse(e$both_sig, 0.9, 0.45)
      plot(e$dbeta_promoter, e$dbeta_body, pch = 19, cex = cex, col = col,
           xlim = c(-lim, lim), ylim = c(-lim, lim),
           xlab = "delta beta, promoter", ylab = "delta beta, gene body")
      graphics::abline(h = 0, v = 0, col = "#9aa5b1")
      graphics::abline(0, 1, lty = 3, col = "#c2cad2")
      ## The off-diagonal quadrants are the point, so they are named on the plot.
      graphics::mtext("promoter up / body down", side = 3, adj = 1, cex = .7,
                      col = MP_PAL$grey)
      graphics::mtext("promoter down / body up", side = 1, adj = 0, cex = .7,
                      col = MP_PAL$grey, line = 2.6)
      lab <- utils::head(e[e$both_sig & e$opposite, , drop = FALSE], 8)
      if (nrow(lab)) {
        graphics::text(lab$dbeta_promoter, lab$dbeta_body, lab$gene,
                       pos = 3, cex = 0.7, col = "#333333", xpd = NA)
      }
      graphics::mtext(sprintf("%s genes with both features  -  %s significant in both, %s of those opposing",
                              format(nrow(e), big.mark = ","),
                              format(sum(e$both_sig), big.mark = ","),
                              format(sum(e$both_sig & e$opposite), big.mark = ",")),
                      side = 3, adj = 0, cex = .72, col = MP_PAL$grey)
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = plot_data, renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot, res = c(90, 130),
      pdf.width = 7, pdf.height = 6, add.watermark = watermark)
  })
}

methylome_table_divergence_ui <- function(id, title, info.text, caption, height, width) {
  ns <- shiny::NS(id)
  TableModuleUI(
    translate = FALSE,ns("tblmod"), title = title, info.text = info.text,
                caption = caption, height = height, width = width, label = "b")
}

methylome_table_divergence_server <- function(id, r.ewas, r.thresh, scrollY = "22vh") {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      e <- mp_divergence(r.ewas(), r.thresh())
      data.frame(
        Gene = e$gene,
        `Promoter` = round(e$dbeta_promoter, 4),
        `Body` = round(e$dbeta_body, 4),
        `MeGDP` = round(e$megdp, 4),
        `Opposing` = ifelse(e$opposite, "yes", "-"),
        `FDR promoter` = e$q_promoter,
        `FDR body` = e$q_body,
        `Probes promoter` = e$n_promoter,
        `Probes body` = e$n_body,
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })
    render <- function(sy) {
      dt <- table_data(); shiny::req(dt)
      DT::datatable(dt, class = "compact hover", rownames = FALSE,
        extensions = c("Scroller"), plugins = "scrollResize", selection = "none",
        options = list(dom = "lfrtip", scroller = TRUE, scrollX = TRUE,
                       scrollY = sy, scrollResize = TRUE, deferRender = TRUE)) |>
        DT::formatSignif(c("FDR promoter", "FDR body"), digits = 3) |>
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") |>
        DT::formatStyle("Opposing",
          color = DT::styleEqual("yes", MP_PAL$hyper), fontWeight = "bold") |>
        DT::formatStyle(c("Promoter", "Body"),
          color = DT::styleInterval(0, c(MP_PAL$hypo, MP_PAL$hyper)))
    }
    TableModuleServer("tblmod", func = function() render(scrollY),
                      func2 = function() render("60vh"),
                      csvFunc = table_data, selector = "none")
  })
}

## ----------------------------------------------------- trait frequency (C) ---

## What the hit list, taken together, is already known for. The per-CpG column
## cannot answer this: it keeps the top traits per CpG and collapses the rest,
## so a trait sitting fourth on a thousand CpGs is invisible there and top of
## the chart here.
methylome_plot_traitfreq_ui <- function(id, title, caption, info.text, info.methods,
                                        height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    translate = FALSE,ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    options = shiny::tagList(
      shiny::sliderInput(ns("topn"), "Traits shown:", min = 5, max = 30, value = 15, step = 5)
    ),
    height = height, width = width, label = "e")
}

methylome_plot_traitfreq_server <- function(id, r.catalog, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      cc <- r.catalog()
      ## Refuse with the catalog's own reason - "collapsed fit" and "no hits"
      ## are different problems and the user needs to know which one this is.
      shiny::validate(shiny::need(isTRUE(cc$ok), cc$reason))
      shiny::validate(shiny::need(nrow(cc$traits) > 0,
        "None of the CpGs passing this threshold is in the EWAS Catalog."))
      ## Count CpGs per trait, not studies: "how much of my hit list is this
      ## trait" is the question. A single heavily studied CpG would otherwise
      ## dominate the chart on its own.
      tab <- sort(table(droplevels(cc$traits$trait)), decreasing = TRUE)
      data.frame(trait = names(tab), n_cpgs = as.integer(tab),
                 pct = 100 * as.integer(tab) / length(cc$probes),
                 stringsAsFactors = FALSE)
    })
    plot.RENDER <- function() {
      e <- plot_data(); shiny::req(nrow(e) > 0)
      cc <- r.catalog(); shiny::req(isTRUE(cc$ok))
      k <- min(nrow(e), if (is.null(input$topn)) 15L else as.integer(input$topn))
      e <- utils::head(e, k)
      op <- graphics::par(mar = c(4, 12, 2.2, 1.5), las = 1, mgp = c(2.2, 0.6, 0))
      on.exit(graphics::par(op))
      graphics::barplot(rev(e$n_cpgs), horiz = TRUE, col = MP_PAL$hyper, border = NA,
        names.arg = rev(substr(e$trait, 1, 34)), cex.names = 0.7, cex.axis = 0.8,
        xlab = "CpGs in your hit list reported for this trait")
      ## The denominator, always. Without it the bars read as "these traits
      ## explain my result" when most of the hit list may simply be absent from
      ## the catalog - which is not the same as novel.
      n_known <- length(unique(as.character(cc$traits$cpg)))
      graphics::mtext(sprintf("%s of %s hits appear in the catalog (%.0f%%)",
                              format(n_known, big.mark = ","),
                              format(length(cc$probes), big.mark = ","),
                              100 * n_known / length(cc$probes)),
                      side = 3, adj = 1, cex = .75, col = MP_PAL$grey)
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = plot_data, renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot, res = c(90, 130),
      pdf.width = 6, pdf.height = 5, add.watermark = watermark)
  })
}

## ------------------------------------------------------- trait breakdown ---

## The tabular twin of the barplot, and the answer to "which traits exactly?".
## Self-contained: everything it needs is the hit list's own catalog slice, so
## it needs no selection made on another sub-tab.
methylome_table_traits_ui <- function(id, title, info.text, caption, height, width) {
  ns <- shiny::NS(id)
  TableModuleUI(
    translate = FALSE,ns("tblmod"), title = title, info.text = info.text,
                caption = caption, height = height, width = width, label = "f")
}

methylome_table_traits_server <- function(id, r.catalog, scrollY = "22vh") {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      cc <- r.catalog()
      shiny::validate(shiny::need(isTRUE(cc$ok), cc$reason))
      shiny::validate(shiny::need(nrow(cc$traits) > 0,
        "None of the CpGs passing this threshold is in the EWAS Catalog."))
      d <- cc$traits
      n_hits <- length(cc$probes)
      ## Two different counts, and the difference matters: how many of YOUR
      ## hits carry the trait, against how many study analyses reported it at
      ## all. A trait can be heavily studied on one CpG, or lightly studied
      ## across hundreds of them.
      agg <- stats::aggregate(cbind(cpgs = rep(1L, nrow(d)), studies = d$n),
                              by = list(Trait = as.character(d$trait)), FUN = sum)
      agg <- agg[order(-agg$cpgs, -agg$studies), , drop = FALSE]
      data.frame(
        Trait = agg$Trait,
        CpGs = as.integer(agg$cpgs),
        `% of hits` = round(100 * agg$cpgs / n_hits, 1),
        Studies = as.integer(agg$studies),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })
    render <- function(sy) {
      dt <- table_data(); shiny::req(dt)
      show <- dt
      show$Trait <- paste0(htmltools::htmlEscape(show$Trait),
                           mp_ewas_link("trait", show$Trait))
      DT::datatable(show, class = "compact hover", rownames = FALSE, escape = FALSE,
        extensions = c("Scroller"), plugins = "scrollResize", selection = "none",
        options = list(dom = "lfrtip", scroller = TRUE, scrollX = TRUE,
                       scrollY = sy, scrollResize = TRUE, deferRender = TRUE,
                       order = list(list(1, "desc")))) |>
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") |>
        DT::formatStyle("CpGs",
          background = DT::styleColorBar(c(0, max(dt$CpGs, 1L)), "#dce7f2"),
          backgroundSize = "98% 60%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center") |>
        DT::formatStyle("Studies",
          background = DT::styleColorBar(c(0, max(dt$Studies, 1L)), "#e6dcea"),
          backgroundSize = "98% 60%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center")
    }
    TableModuleServer("tblmod", func = function() render(scrollY),
                      func2 = function() render("60vh"),
                      csvFunc = table_data, selector = "none")
  })
}

## ------------------------------------------------------------- hits table ---

methylome_table_hits_ui <- function(id, title, info.text, caption, height, width) {
  ns <- shiny::NS(id)
  TableModuleUI(
    translate = FALSE,ns("tblmod"), title = title, info.text = info.text,
                caption = caption, height = height, width = width, label = "d")
}

methylome_table_hits_server <- function(id, r.ewas, r.thresh, r.catalog,
                                        scrollY = "22vh") {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      res <- r.ewas(); d <- res$data
      anova_fit <- isTRUE(res$meta$anova)
      sig <- mp_ewas_sig(d, r.thresh())
      shiny::validate(shiny::need(any(sig),
        "No CpG passes the current threshold. Relax the cut-off in the settings panel."))
      d <- d[sig, , drop = FALSE]
      d <- d[order(d$p), , drop = FALSE]
      out <- data.frame(
        CpG = d$probe,
        Gene = ifelse(is.na(d$gene) | d$gene == "", "-", d$gene),
        Chr = d$chr,
        Position = if (is.null(d$pos_start)) d$pos else
          ifelse(is.na(d$pos_start), NA_character_,
                 ifelse(d$pos_start == d$pos_end, format(d$pos_start, big.mark = ",", trim = TRUE),
                        paste0(format(d$pos_start, big.mark = ",", trim = TRUE), "-",
                               format(d$pos_end, big.mark = ",", trim = TRUE)))),
        ## On a collapsed fit these three describe an aggregate, so each says
        ## how well it does: Context carries the share of probes that actually
        ## carry it, Region is the partition the feature belongs to rather than
        ## the commonest RefGene sub-group, and Position is the span.
        Context = if (is.null(d$island_n)) {
          ifelse(is.na(d$island), "-", d$island)
        } else {
          ifelse(is.na(d$island), "-",
                 sprintf("%s %d/%d", d$island, d$island_n, d$n_probes))
        },
        Region = if (!is.null(d$region)) d$region else
          ifelse(is.na(d$location) | d$location == "", "-", d$location),
        ## Unsigned across k groups, so it is named for what it is.
        effect = d$dbeta,
        Direction = if (anova_fit) "-" else
          ifelse(is.na(d$dbeta), "-", ifelse(d$dbeta > 0, "hyper", "hypo")),
        `P value` = d$p, FDR = d$q,
        check.names = FALSE, stringsAsFactors = FALSE
      )
      colnames(out)[colnames(out) == "effect"] <-
        if (anova_fit) "Beta range" else "Delta beta"
      if (anova_fit && !is.null(d$Fstat)) out$F <- signif(d$Fstat, 3)
      ## Which traits each hit is already known for. Absent from the catalog is
      ## absent from the catalog, not novel biology, so it says so in words.
      ## Catalog columns only when the catalog can actually answer. A collapsed
      ## fit still gets its full hit list; it just has no CpG ids to look up.
      cc <- r.catalog()
      if (isTRUE(cc$ok)) {
        cat <- mp_catalog_summary(cc$traits, d$probe, max_traits = 2)
        out$Reported <- cat$reported
        out$Studies <- cat$studies
      }
      ## Collapsed fits only: how many probes went into each feature.
      if (!is.null(d$n_probes)) out$Probes <- as.integer(d$n_probes)
      ## What the eye needs first is what the model found: effect, direction and
      ## significance, straight after the identity. Position and context sit in
      ## the middle; the catalog evidence goes last, where it can be as wide as
      ## it likes without pushing the statistics off screen.
      eff_name <- if (anova_fit) "Beta range" else "Delta beta"
      ## "CpG" is the wrong word for a collapsed feature; the column is named
      ## for whatever was actually tested.
      id_name <- if (identical(res$meta$collapse, "probe")) "CpG" else "Gene region"
      colnames(out)[colnames(out) == "CpG"] <- id_name
      front <- c(id_name, "Gene", eff_name, "Direction",
                 if (anova_fit) "F", "P value", "FDR",
                 if ("Probes" %in% colnames(out)) "Probes")
      back <- intersect(c("Studies", "Reported"), colnames(out))
      mid <- setdiff(colnames(out), c(front, back))
      out <- out[, c(front, mid, back), drop = FALSE]
      rownames(out) <- NULL
      out
    })
    render <- function(sy) {
      dt <- table_data(); shiny::req(dt)
      ## Links are built here, not in table_data, so the CSV download stays
      ## plain text rather than a column of anchor tags.
      has_cat <- "Reported" %in% colnames(dt)
      if (has_cat) dt$Reported <- mp_catalog_reported_html(dt$Reported, dt$CpG)
      ## Only a real CpG has a catalog page to link to.
      if (has_cat) dt$CpG <- paste0(htmltools::htmlEscape(dt$CpG),
                                    mp_ewas_link("cpg", dt$CpG))
      ## Probes-per-feature, bar-scaled so a 400-probe body reads differently
      ## from a 3-probe promoter at a glance.
      has_np <- "Probes" %in% colnames(dt)
      eff <- if ("Beta range" %in% colnames(dt)) "Beta range" else "Delta beta"
      DT::datatable(dt, class = "compact hover", rownames = FALSE, escape = FALSE,
        extensions = c("Buttons", "Scroller"), plugins = "scrollResize",
        selection = list(mode = "single", target = "row"),
        options = list(dom = "lfrtip", scroller = TRUE, scrollX = TRUE,
                       scrollY = sy, scrollResize = TRUE, deferRender = TRUE,
                       ## Reported is the widest column by far. It is not
                       ## character-truncated: the cell is now anchor markup and
                       ## trunc_display_row() substr's the raw string, which
                       ## would cut a tag in half. Two traits plus a "+N more"
                       ## link bounds it instead, and the row click has the rest.
                       columnDefs = if (has_cat) list(list(
                         targets = which(colnames(dt) == "Reported") - 1L,
                         width = "260px")) else list(),
                       order = list(list(which(colnames(dt) == "P value") - 1L, "asc")))) |>
        ## signif() in the data would leave a double that DT re-serialises at
        ## full precision - 6.049999999999999e-64 - so the rounding happens at
        ## render and the column stays numeric, which is what DT sorts on.
        DT::formatSignif(intersect(c("P value", "FDR", "F"), colnames(dt)), digits = 3) |>
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") |>
        DT::formatStyle("Direction",
          color = DT::styleEqual(c("hyper", "hypo"), c(MP_PAL$hyper, MP_PAL$hypo)),
          fontWeight = "bold") |>
        DT::formatStyle(eff,
          background = DT::styleColorBar(c(-1, 1) * max(abs(dt[[eff]]), na.rm = TRUE),
                                         "#dce7f2"),
          backgroundSize = "98% 60%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center") |>
        ## How hard anyone has looked at this CpG. Sort ascending for the hits
        ## nobody has reported yet.
        DT::formatStyle(if (has_np) "Probes" else eff,
          background = DT::styleColorBar(
            if (has_np) c(0, max(dt$Probes, 1L)) else
              c(-1, 1) * max(abs(dt[[eff]]), na.rm = TRUE), "#e8e2d6"),
          backgroundSize = "98% 60%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center") |>
        DT::formatStyle(if (has_cat) "Studies" else eff,
          background = DT::styleColorBar(
            if (has_cat) c(0, max(dt$Studies, 1L)) else
              c(-1, 1) * max(abs(dt[[eff]]), na.rm = TRUE),
            if (has_cat) "#e6dcea" else "#dce7f2"),
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
  PlotModuleUI(
    translate = FALSE,ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "e")
}

methylome_plot_stripcharts_server <- function(id, pgx, r.ewas, r.thresh,
                                              r.topn = shiny::reactive(6),
                                              watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      p <- pgx(); res <- r.ewas(); m <- res$meta
      d <- res$data
      d <- d[mp_ewas_sig(d, r.thresh()), , drop = FALSE]
      shiny::validate(shiny::need(nrow(d) > 0,
        "No CpG passes the current threshold. Relax the cut-off in the settings panel."))
      n <- r.topn(); if (is.null(n) || is.na(n) || n < 1) n <- 6
      d <- d[order(d$p), , drop = FALSE][seq_len(min(n, nrow(d))), , drop = FALSE]
      X <- mp_feature_beta(p, m, d$probe)
      d <- d[d$probe %in% rownames(X), , drop = FALSE]
      shiny::validate(shiny::need(nrow(d) > 0, "No feature left to draw."))
      ss <- intersect(m$samples, colnames(X))
      shiny::validate(shiny::need(length(ss) > 2, "No samples left for this model."))
      list(d = d, X = X[d$probe, ss, drop = FALSE], grp = droplevels(m$strata[ss]),
           x = if (isTRUE(m$continuous)) m$x[ss] else NULL,
           anova = isTRUE(m$anova), outcome = res$contrast)
    })
    plot.RENDER <- function() {
      res <- plot_data(); d <- res$d; X <- res$X; grp <- res$grp
      n <- nrow(d); nc <- min(3, n); nr <- ceiling(n / nc)
      op <- graphics::par(mfrow = c(nr, nc), mar = c(3.4, 3.2, 2.2, 0.6),
                          las = 1, mgp = c(2, 0.6, 0), cex.axis = 0.85)
      on.exit(graphics::par(op))
      lv <- levels(grp)
      cols <- grDevices::colorRampPalette(c(MP_PAL$hypo, MP_PAL$hyper))(max(2, length(lv)))
      lab_i <- function(i) {
        if (!is.na(d$gene[i]) && d$gene[i] != "") {
          paste0(d$probe[i], "  (", d$gene[i], ")")
        } else d$probe[i]
      }
      for (i in seq_len(n)) {
        y <- as.numeric(X[i, ])
        if (!is.null(res$x)) {
          ## Continuous outcome: the regression itself is the honest picture,
          ## and shows whether the slope rides on a couple of extreme samples.
          graphics::plot(res$x, y, pch = 19, cex = 0.75, ylim = c(0, 1),
                         col = grDevices::adjustcolor(MP_PAL$hypo, 0.6),
                         xlab = res$outcome, ylab = "beta")
          graphics::abline(h = c(0.2, 0.8), lty = 3, col = "#c9ccc7")
          good <- is.finite(res$x) & is.finite(y)
          if (sum(good) > 2) {
            graphics::abline(stats::lm(y[good] ~ res$x[good]),
                             col = MP_PAL$hyper, lwd = 1.6)
          }
        } else {
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
        }
        graphics::title(main = lab_i(i), cex.main = 0.95, font.main = 1)
        ## The k-group effect is a range, so a signed "+0.12" would read as a
        ## direction the F-test never established.
        graphics::mtext(
          if (res$anova) sprintf("beta range %.3f   q %.2g", d$dbeta[i], d$q[i])
          else sprintf("d-beta %+.3f   q %.2g", d$dbeta[i], d$q[i]),
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

## ------------------------------------------------------------------ volcano --
##
## Rated near-universal in both the population-EWAS and the cancer literature,
## and the one view that shows the effect-size/significance trade-off directly.
## Both axes come from the fitted model, so this cannot be borrowed from the
## platform's own volcano, which draws the unadjusted result stored at upload.

methylome_plot_volcano_ui <- function(id, title, caption, info.text, info.methods,
                                      height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    translate = FALSE,ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "b")
}

methylome_plot_volcano_server <- function(id, r.ewas, r.thresh, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      res <- r.ewas(); th <- r.thresh()
      ## Refuse rather than approximate: the x axis of a volcano is a signed
      ## effect, and a difference among three or more groups does not have one.
      ## Plotting the unsigned beta range here would look like a volcano and
      ## mean something else.
      shiny::validate(shiny::need(!isTRUE(res$meta$anova),
        "This model is an F-test across three or more groups, which has no single direction - so there is no volcano to draw. The stripcharts on the QQ & context sub-tab show what each group actually does at the top CpGs."))
      d <- res$data
      d <- d[is.finite(d$p) & is.finite(d$dbeta), , drop = FALSE]
      shiny::validate(shiny::need(nrow(d) > 0, "No probe has both an effect size and a p-value."))
      sig <- mp_ewas_sig(d, th)
      op <- graphics::par(mar = c(4.2, 4.5, 1.4, 1), las = 1, mgp = c(2.4, 0.6, 0))
      on.exit(graphics::par(op))
      ## For a continuous outcome delta-beta is a slope, so say so on the axis.
      xlab <- if (isTRUE(res$meta$continuous)) {
        paste0("delta-beta per unit ", res$contrast)
      } else {
        "delta-beta"
      }
      graphics::plot(d$dbeta, -log10(d$p), pch = 19, cex = 0.35, col = "#c2ccd6",
                     xlab = xlab, ylab = expression(-log[10](p)))
      graphics::abline(v = 0, lty = 3, col = MP_PAL$grey)
      if (any(sig)) {
        graphics::points(d$dbeta[sig], -log10(d$p[sig]), pch = 19, cex = 0.6,
          col = ifelse(d$dbeta[sig] > 0, MP_PAL$hyper, MP_PAL$hypo))
      }
      h <- mp_ewas_hline(d, th)
      if (is.finite(h)) graphics::abline(h = h, lty = 2, col = MP_PAL$grey)
      if (!is.null(th$min_dbeta) && th$min_dbeta > 0) {
        graphics::abline(v = c(-1, 1) * th$min_dbeta, lty = 2, col = MP_PAL$grey)
      }
      if (any(sig)) {
        s <- d[sig, , drop = FALSE]
        s <- s[order(s$p), , drop = FALSE]
        s <- utils::head(s, 8)
        lab <- ifelse(is.na(s$gene) | s$gene == "", s$probe, sub(";.*", "", s$gene))
        graphics::text(s$dbeta, -log10(s$p), labels = lab, pos = 3, cex = 0.7,
                       col = "#33404d", xpd = NA)
      }
      graphics::mtext(sprintf("%s hits", format(sum(sig), big.mark = ",")),
                      side = 3, adj = 1, cex = 0.75, col = MP_PAL$grey)
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive(r.ewas()$data), renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot, res = c(90, 130),
      pdf.width = 6, pdf.height = 6, add.watermark = watermark)
  })
}

## ------------------------------------------------------------- hits heatmap --
##
## Samples x top CpGs, samples ordered by the model's own strata. The platform's
## clustering heatmap selects features by variance across the whole matrix; here
## the selection IS the analysis - these are the CpGs the adjusted model called -
## so it cannot be borrowed either.

methylome_plot_hitmap_ui <- function(id, title, caption, info.text, info.methods,
                                     height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    translate = FALSE,ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "c")
}

methylome_plot_hitmap_server <- function(id, pgx, r.ewas, r.thresh,
                                         n.max = 50, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      p <- pgx(); res <- r.ewas(); m <- res$meta
      d <- res$data
      sig <- mp_ewas_sig(d, r.thresh())
      shiny::validate(shiny::need(sum(sig) >= 2,
        "Fewer than two CpGs pass the current threshold - nothing to cluster."))
      d <- d[sig, , drop = FALSE]
      d <- utils::head(d[order(d$p), , drop = FALSE], n.max)
      X <- mp_feature_beta(p, m, d$probe)
      d <- d[d$probe %in% rownames(X), , drop = FALSE]
      shiny::validate(shiny::need(nrow(d) >= 2,
        "Fewer than two features could be traced back to the methylation matrix."))
      ss <- intersect(m$samples, colnames(X))
      shiny::validate(shiny::need(length(ss) > 1, "No samples left for this model."))
      grp <- droplevels(m$strata[ss])
      ord <- order(grp)
      list(B = X[d$probe, ss[ord], drop = FALSE], grp = grp[ord], d = d)
    })
    plot.RENDER <- function() {
      z <- plot_data(); B <- z$B; grp <- z$grp
      lv <- levels(grp)
      gcol <- grDevices::colorRampPalette(c(MP_PAL$hypo, MP_PAL$hyper))(max(2, length(lv)))
      op <- graphics::par(mar = c(1.5, 8, 3.4, 1), las = 1)
      on.exit(graphics::par(op))
      nr <- nrow(B); nc <- ncol(B)
      ## Beta on a fixed [0,1] scale rather than a row z-score: the absolute
      ## methylation level is the interpretable quantity, and a z-score makes a
      ## 2% difference look identical to a 60% one.
      graphics::image(x = seq_len(nc), y = seq_len(nr),
                      z = t(B[rev(seq_len(nr)), , drop = FALSE]),
                      zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "",
                      col = grDevices::colorRampPalette(
                        c("#2166ac", "#f7f7f7", "#b2182b"))(64))
      lab <- ifelse(is.na(z$d$gene) | z$d$gene == "", z$d$probe,
                    paste0(sub(";.*", "", z$d$gene), "  ", z$d$probe))
      graphics::axis(2, at = seq_len(nr), labels = rev(lab), tick = FALSE,
                     cex.axis = 0.5, line = -0.6)
      ## Group bar above the matrix.
      bh <- max(1, nr * 0.035)
      graphics::rect(seq_len(nc) - 0.5, nr + 0.6, seq_len(nc) + 0.5, nr + 0.6 + bh,
                     col = gcol[as.integer(grp)], border = NA, xpd = NA)
      graphics::legend("top", horiz = TRUE, bty = "n", inset = c(0, -0.13),
                       legend = lv, fill = gcol[seq_along(lv)], border = NA,
                       cex = 0.75, xpd = NA,
                       title = if (isTRUE(r.ewas()$meta$continuous)) "outcome tertile" else NULL)
      graphics::mtext(sprintf("%d samples", nc), side = 1, adj = 1, cex = 0.7,
                      col = MP_PAL$grey, line = 0.2)
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive({
        z <- plot_data(); cbind(probe = rownames(z$B), as.data.frame(z$B))
      }),
      renderFunc = shiny::renderPlot, renderFunc2 = shiny::renderPlot,
      res = c(90, 130), pdf.width = 8, pdf.height = 7, add.watermark = watermark)
  })
}
