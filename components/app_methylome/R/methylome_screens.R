##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Methylome Profiler — the four screens. Each screen is a plain Shiny module so
## the same code serves the bigdash shell inside Playground and the standalone
## harness in ../test/app.R.
##
## Everything here reads the loaded pgx and nothing else: pgx$X (beta),
## pgx$genes (450K/EPIC annotation) and pgx$gx.meta (the EWAS already computed
## at upload). No screen needs a contrast except EWAS.

MP_PAL <- list(
  hypo = "#4575b4", hyper = "#d73027", grey = "#8697a6",
  ok = "#18753a", bad = "#a8202c", warn = "#8a5a06"
)

## ---------------------------------------------------------------- helpers ---

mp_beta <- function(pgx) {
  X <- pgx$X
  if (max(X, na.rm = TRUE) > 1 || min(X, na.rm = TRUE) < 0) X <- playbase::mToBeta(X)
  X
}

## Fraction of probes outside the intermediate band. A healthy methylome is
## bimodal at ~0.2/0.8, so this sits near 0.85-0.95; a flat sample drops.
mp_bimodality <- function(X) {
  round(1 - colMeans(X > 0.3 & X < 0.7, na.rm = TRUE), 3)
}

## Imprinted DMRs sit near 0.5. Deviation means normalisation went wrong.
mp_imprint_drift <- function(X) {
  imp <- tryCatch({
    e <- new.env(); utils::data("iDMR", package = "wateRmelon", envir = e)
    get("iDMR", envir = e)
  }, error = function(e) character(0))
  ii <- intersect(rownames(X), imp)
  if (length(ii) < 10) return(rep(NA_real_, ncol(X)))
  round(colMeans(abs(X[ii, , drop = FALSE] - 0.5), na.rm = TRUE), 3)
}

mp_clocks <- function(X) {
  ## agep() calls data("age_coefficients") with no package= argument, so the
  ## coefficients only resolve when wateRmelon is attached. Calling it
  ## namespace-qualified silently yields NA ages.
  if (!"package:wateRmelon" %in% search()) {
    suppressPackageStartupMessages(library(wateRmelon))
  }
  A <- tryCatch(wateRmelon::agep(X, method = "all"), error = function(e) NULL)
  if (is.null(A)) return(NULL)
  A <- as.data.frame(A)
  keep <- grep("\\.age$", colnames(A), value = TRUE)
  out <- A[, keep, drop = FALSE]
  colnames(out) <- sub("\\..*", "", keep)
  round(out, 1)
}

mp_sex <- function(X) {
  chr <- NULL
  s <- tryCatch(wateRmelon::estimateSex(X, do_plot = FALSE), error = function(e) NULL)
  if (is.null(s) || !"predicted_sex" %in% names(s)) return(NULL)
  s
}

## Does the dataset even carry sex chromosomes? GEO matrices are often shipped
## with filterXY already applied, in which case we refuse to predict.
mp_has_xy <- function(pgx) {
  cc <- grep("^chr$|chrom", tolower(colnames(pgx$genes)))
  if (!length(cc)) return(FALSE)
  v <- as.character(pgx$genes[[cc[1]]])
  any(grepl("^(chr)?[XY]$", v, ignore.case = TRUE))
}

mp_context_col <- function(pgx, what = c("island", "gene")) {
  what <- match.arg(what)
  cn <- colnames(pgx$genes)
  hit <- if (what == "island") grep("Relation_to_Island", cn, ignore.case = TRUE)
         else grep("genomic_location|RefGene_Group", cn, ignore.case = TRUE)
  if (!length(hit)) return(NULL)
  as.character(pgx$genes[[cn[hit[1]]]])
}

mp_card <- function(title, ..., right = NULL) {
  shiny::div(
    class = "mp-card",
    shiny::div(class = "mp-card-head", title, shiny::span(class = "mp-sp"), right),
    shiny::div(class = "mp-card-body", ...)
  )
}

mp_tile <- function(k, v, d) {
  shiny::div(class = "mp-tile",
    shiny::div(class = "mp-k", k), shiny::div(class = "mp-v", v), shiny::div(class = "mp-d", d))
}

## ------------------------------------------------------- 1. sample ledger ---

methylome_ledger_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("tiles")),
    shiny::div(class = "mp-row2",
      mp_card("Per-sample checks", DT::DTOutput(ns("tbl"))),
      mp_card("Beta distribution — all samples", shiny::plotOutput(ns("dens"), height = "300px"))
    )
  )
}

methylome_ledger_server <- function(id, pgx) {
  shiny::moduleServer(id, function(input, output, session) {
    dat <- shiny::reactive({
      shiny::req(pgx())
      p <- pgx(); X <- mp_beta(p)
      bim <- mp_bimodality(X)
      drift <- mp_imprint_drift(X)
      miss <- round(100 * colMeans(is.na(X)), 2)
      clk <- mp_clocks(X)
      hasxy <- mp_has_xy(p)
      sex <- if (hasxy) mp_sex(X) else NULL
      data.frame(
        Sample = colnames(X),
        Bimodality = bim,
        `Missing %` = miss,
        `DNAm age` = if (!is.null(clk)) clk[[1]] else NA,
        `Sex predicted` = if (!is.null(sex)) as.character(sex$predicted_sex) else "not available",
        `Imprint drift` = drift,
        ## Absolute cut-offs do not transfer between platforms or tissues, so
        ## flag relative to the cohort: 3 MADs above the median is an outlier
        ## here, not a universal threshold.
        Verdict = {
          out_drift <- !is.na(drift) &
            drift > median(drift, na.rm = TRUE) + 3 * mad(drift, na.rm = TRUE)
          out_bim <- bim < median(bim, na.rm = TRUE) - 3 * mad(bim, na.rm = TRUE)
          ifelse(out_drift | out_bim, "CHECK", "PASS")
        },
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })

    output$tiles <- shiny::renderUI({
      d <- dat()
      shiny::div(class = "mp-row4",
        mp_tile("Samples", nrow(d), "in this dataset"),
        mp_tile("Passing", sum(d$Verdict == "PASS"), "all checks clear"),
        mp_tile("Flagged", sum(d$Verdict == "CHECK"), "bimodality or drift"),
        mp_tile("Sex prediction",
                if (all(d$`Sex predicted` == "not available")) "n/a" else "ok",
                if (all(d$`Sex predicted` == "not available")) "no X/Y probes in data" else "from X/Y probes")
      )
    })

    output$tbl <- DT::renderDT({
      d <- dat()
      DT::datatable(d, rownames = FALSE, selection = "none",
        options = list(pageLength = 8, dom = "tp", scrollX = TRUE)) |>
        DT::formatStyle("Verdict", color = DT::styleEqual(c("PASS", "CHECK"),
          c(MP_PAL$ok, MP_PAL$warn)), fontWeight = "bold")
    })

    output$dens <- shiny::renderPlot({
      p <- pgx(); X <- mp_beta(p)
      set.seed(1)
      ii <- sample(nrow(X), min(6000, nrow(X)))
      op <- par(mar = c(4, 4, 1, 1), las = 1); on.exit(par(op))
      plot(NA, xlim = c(0, 1), ylim = c(0, 4.2), xlab = "Beta value", ylab = "Density")
      for (j in seq_len(min(24, ncol(X)))) {
        lines(density(X[ii, j], na.rm = TRUE, from = 0, to = 1),
              col = adjustcolor(MP_PAL$grey, 0.55), lwd = 1)
      }
      abline(v = c(0.2, 0.8), lty = 2, col = c(MP_PAL$hypo, MP_PAL$hyper))
      legend("topright", bty = "n", lty = 2, col = c(MP_PAL$hypo, MP_PAL$hyper),
             legend = c("hypo 0.2", "hyper 0.8"))
    })
  })
}

## ------------------------------------------------------ 2. epigenetic age ---

methylome_age_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("tiles")),
    shiny::div(class = "mp-row2",
      mp_card("Clock agreement", shiny::plotOutput(ns("pairs"), height = "330px")),
      mp_card("DNAm age by group", shiny::plotOutput(ns("box"), height = "330px"))
    ),
    mp_card("Per-clock coverage", DT::DTOutput(ns("cov")))
  )
}

methylome_age_server <- function(id, pgx) {
  shiny::moduleServer(id, function(input, output, session) {
    clocks <- shiny::reactive({
      shiny::req(pgx())
      X <- mp_beta(pgx())
      if (!"package:wateRmelon" %in% search()) suppressPackageStartupMessages(library(wateRmelon))
      raw <- as.data.frame(wateRmelon::agep(X, method = "all"))
      A <- mp_clocks(X)
      ## A clock fitted on a fraction of its probes returns a confident, wrong
      ## number (we have seen -13.5 years). Anything missing probes is withheld
      ## rather than shown with a caveat.
      miss <- grep("n_missing$", colnames(raw), value = TRUE)
      nmiss <- if (length(miss)) vapply(raw[miss], max, numeric(1)) else rep(0, ncol(A))
      names(nmiss) <- sub("\\..*", "", miss)
      ok <- names(nmiss)[nmiss == 0]
      list(age = A, raw = raw, nmiss = nmiss, ok = intersect(colnames(A), ok))
    })

    output$tiles <- shiny::renderUI({
      cl <- clocks(); A <- cl$age; shiny::req(A)
      ok <- cl$ok
      shiny::validate(shiny::need(length(ok) > 0,
        "No clock has complete probe coverage in this dataset."))
      h <- A[[ok[1]]]
      conc <- if (length(ok) >= 2) {
        cm <- cor(A[, ok, drop = FALSE]); round(median(cm[upper.tri(cm)]), 2)
      } else "n/a"
      shiny::div(class = "mp-row4",
        mp_tile("Median DNAm age", round(median(h), 1), paste0(ok[1], ", years")),
        mp_tile("Range", paste(round(range(h), 1), collapse = " – "), "years"),
        mp_tile("Clocks usable", paste0(length(ok), " of ", ncol(A)), "full probe coverage"),
        mp_tile("Concordance", conc, "median pairwise r")
      )
    })

    output$pairs <- shiny::renderPlot({
      cl <- clocks(); A <- cl$age; ok <- cl$ok; shiny::req(A)
      shiny::validate(shiny::need(length(ok) >= 2,
        "Need two clocks with full coverage to compare."))
      x <- A[[ok[1]]]; y <- A[[ok[2]]]
      op <- par(mar = c(4, 4, 1, 1), las = 1); on.exit(par(op))
      plot(x, y, pch = 19, col = adjustcolor(MP_PAL$hypo, .7), cex = 1.1,
           xlab = paste(ok[1], "age"), ylab = paste(ok[2], "age"))
      abline(0, 1, lty = 2, col = MP_PAL$grey)
      legend("topleft", bty = "n", legend = sprintf("r = %.2f", cor(x, y)))
    })

    output$box <- shiny::renderPlot({
      cl <- clocks(); A <- cl$age; p <- pgx(); shiny::req(A, length(cl$ok) > 0)
      A <- A[, cl$ok, drop = FALSE]
      grp <- NULL
      for (cn in colnames(p$samples)) {
        v <- as.character(p$samples[[cn]])
        if (length(unique(v[!is.na(v)])) == 2) { grp <- v; gname <- cn; break }
      }
      op <- par(mar = c(4, 4, 1, 1), las = 1); on.exit(par(op))
      if (is.null(grp)) { plot.new(); text(.5, .5, "No two-level phenotype"); return() }
      boxplot(A[[1]] ~ factor(grp), col = c("#dce7f2", "#f2dce2"),
              xlab = gname, ylab = "DNAm age (Horvath)", border = "#444444")
      points(jitter(as.numeric(factor(grp)), 0.6), A[[1]],
             pch = 19, col = adjustcolor(MP_PAL$grey, .6), cex = .8)
    })

    output$cov <- DT::renderDT({
      R <- clocks()$raw; shiny::req(R)
      ages <- grep("\\.age$", colnames(R), value = TRUE)
      miss <- grep("n_missing$", colnames(R), value = TRUE)
      nm <- if (length(miss)) as.integer(sapply(R[miss], max)) else rep(0L, length(ages))
      d <- data.frame(
        Clock = sub("\\..*", "", ages),
        `Median age` = ifelse(nm == 0, sprintf("%.1f", sapply(R[ages], median)), "withheld"),
        `Probes missing` = nm,
        Status = ifelse(nm == 0, "full coverage",
                        paste0(nm, " probes absent — value not shown")),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(d, rownames = FALSE, selection = "none",
                    options = list(dom = "t", paging = FALSE)) |>
        DT::formatStyle("Status", color = DT::styleEqual("full coverage", MP_PAL$ok,
                                                         default = MP_PAL$warn))
    })
  })
}

## -------------------------------------------------- 3. methylome character ---

methylome_character_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("tiles")),
    shiny::div(class = "mp-row2",
      mp_card("Mean beta by relation to CpG island", shiny::plotOutput(ns("island"), height = "330px")),
      mp_card("Mean beta by position in gene", shiny::plotOutput(ns("gene"), height = "330px"))
    )
  )
}

methylome_character_server <- function(id, pgx) {
  shiny::moduleServer(id, function(input, output, session) {
    prof <- function(ctx, X) {
      shiny::req(ctx)
      ctx1 <- sub(";.*", "", ctx)
      ctx1[is.na(ctx1) | ctx1 == ""] <- "unannotated"
      mb <- rowMeans(X, na.rm = TRUE)
      tapply(mb, ctx1, mean, na.rm = TRUE)
    }
    output$tiles <- shiny::renderUI({
      p <- pgx(); X <- mp_beta(p)
      isl <- prof(mp_context_col(p, "island"), X)
      mb <- rowMeans(X, na.rm = TRUE)
      shiny::div(class = "mp-row4",
        mp_tile("Island methylation", round(isl[["Island"]], 3), "mean beta, CpG islands"),
        mp_tile("Open sea", round(isl[["OpenSea"]], 3), "mean beta"),
        mp_tile("Global mean", round(mean(mb), 3), "all probes"),
        mp_tile("Hypermethylated", paste0(round(100 * mean(mb > 0.8), 1), "%"), "probes beta > 0.8")
      )
    })
    ## Colour encodes the beta value itself (blue hypo -> red hyper), not the
    ## row order, so the bars read the same way as the ideogram scale.
    barp <- function(v) {
      ramp <- colorRamp(c(MP_PAL$hypo, "#dcd6c8", MP_PAL$hyper))
      cols <- apply(ramp(pmin(pmax(v, 0), 1)), 1, function(z) rgb(z[1], z[2], z[3], maxColorValue = 255))
      op <- par(mar = c(4, 8, 1, 2), las = 1); on.exit(par(op))
      bp <- barplot(rev(v), horiz = TRUE, xlim = c(0, 1), col = rev(cols),
                    border = NA, xlab = "mean beta")
      text(rev(v) - 0.04, bp, sprintf("%.2f", rev(v)), col = "white", font = 2, cex = .9)
      abline(v = c(0.2, 0.8), lty = 3, col = MP_PAL$grey)
    }
    output$island <- shiny::renderPlot({
      p <- pgx(); v <- prof(mp_context_col(p, "island"), mp_beta(p))
      ord <- c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea")
      barp(v[intersect(ord, names(v))])
    })
    output$gene <- shiny::renderPlot({
      p <- pgx(); v <- prof(mp_context_col(p, "gene"), mp_beta(p))
      ord <- c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR", "unannotated")
      barp(v[intersect(ord, names(v))])
    })
  })
}

## ---------------------------------------------------------------- 4. EWAS ---

methylome_ewas_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("tiles")),
    mp_card("Manhattan — genome-wide significance",
            shiny::plotOutput(ns("manhattan"), height = "300px"),
            right = shiny::uiOutput(ns("pick"), inline = TRUE)),
    shiny::div(class = "mp-row2",
      mp_card("QQ and inflation", shiny::plotOutput(ns("qq"), height = "300px")),
      mp_card("Where the hits sit — enrichment vs tested background",
              shiny::plotOutput(ns("enrich"), height = "300px"))
    )
  )
}

methylome_ewas_server <- function(id, pgx) {
  shiny::moduleServer(id, function(input, output, session) {
    output$pick <- shiny::renderUI({
      p <- pgx(); shiny::req(p$gx.meta)
      shiny::selectInput(session$ns("contrast"), NULL,
                         choices = names(p$gx.meta$meta), width = "260px")
    })
    res <- shiny::reactive({
      p <- pgx(); shiny::req(p$gx.meta)
      cmp <- input$contrast; if (is.null(cmp)) cmp <- names(p$gx.meta$meta)[1]
      m <- p$gx.meta$meta[[cmp]]
      g <- p$genes[rownames(m), , drop = FALSE]
      cc <- grep("^chr$|chrom", tolower(colnames(g)))[1]
      pc <- grep("^pos$|position", tolower(colnames(g)))[1]
      data.frame(
        probe = rownames(m),
        p = as.numeric(m$meta.p), q = as.numeric(m$meta.q), fx = as.numeric(m$meta.fx),
        ## pgx$genes$chr is a cytoband ("16q12.2"), not a bare chromosome, so
        ## strip the arm and band exactly as board.epigenomics does.
        chr = if (!is.na(cc)) sub("(p|q|cen).*", "", sub("^chr", "", as.character(g[[cc]]))) else NA,
        pos = if (!is.na(pc)) suppressWarnings(as.numeric(sub(";.*", "", g[[pc]]))) else NA,
        island = sub(";.*", "", mp_context_col(p, "island")[match(rownames(m), rownames(p$genes))]),
        stringsAsFactors = FALSE
      )
    })
    output$tiles <- shiny::renderUI({
      d <- res()
      lam <- round(median(qchisq(1 - d$p, 1), na.rm = TRUE) / qchisq(0.5, 1), 3)
      shiny::div(class = "mp-row4",
        mp_tile("Probes tested", format(nrow(d), big.mark = ","), "after filtering"),
        mp_tile("FDR < 0.05", format(sum(d$q < 0.05, na.rm = TRUE), big.mark = ","), "significant CpGs"),
        mp_tile("Hypermethylated", sum(d$q < 0.05 & d$fx > 0, na.rm = TRUE), "of the hits"),
        mp_tile("Inflation", paste0("λ = ", lam), "inflation or true signal")
      )
    })
    output$manhattan <- shiny::renderPlot({
      d <- res()
      d <- d[!is.na(d$chr) & !is.na(d$pos) & is.finite(d$pos) & is.finite(d$p), ]
      shiny::validate(shiny::need(nrow(d) > 0, "No probes carry chromosome and position."))
      ## Only chromosomes actually present, so the cumulative offset can never
      ## pick up an NA from an empty factor level.
      chrs <- intersect(c(as.character(1:22), "X", "Y"), unique(d$chr))
      d <- d[d$chr %in% chrs, ]
      d$chr <- factor(d$chr, levels = chrs)
      d <- d[order(d$chr, d$pos), ]
      mx <- vapply(split(d$pos, d$chr), function(z) if (length(z)) max(z) else 0, numeric(1))
      off <- c(0, head(cumsum(as.numeric(mx) * 1.02), -1))
      names(off) <- chrs
      d$gx <- d$pos + off[as.character(d$chr)]
      op <- par(mar = c(4, 4.5, 1, 1), las = 1); on.exit(par(op))
      cols <- ifelse(as.integer(d$chr) %% 2 == 0, "#a9b6c4", "#7f8f9f")
      plot(d$gx, -log10(d$p), pch = 19, cex = .35, col = cols, xaxt = "n",
           xlab = "chromosome", ylab = expression(-log[10](p)))
      sig <- d$q < 0.05
      points(d$gx[sig], -log10(d$p[sig]), pch = 19, cex = .6,
             col = ifelse(d$fx[sig] > 0, MP_PAL$hyper, MP_PAL$hypo))
      mid <- tapply(d$gx, d$chr, median, na.rm = TRUE)
      axis(1, at = mid, labels = levels(d$chr), tick = FALSE, cex.axis = .7)
      abline(h = -log10(1e-7), lty = 2, col = MP_PAL$grey)
    })
    output$qq <- shiny::renderPlot({
      d <- res(); pv <- sort(d$p[!is.na(d$p)])
      ex <- -log10(ppoints(length(pv))); ob <- -log10(pv)
      lam <- median(qchisq(1 - d$p, 1), na.rm = TRUE) / qchisq(0.5, 1)
      op <- par(mar = c(4, 4.5, 1, 1), las = 1); on.exit(par(op))
      plot(ex, ob, pch = 19, cex = .35, col = MP_PAL$grey,
           xlab = "expected", ylab = "observed")
      abline(0, 1, lty = 2, col = MP_PAL$hypo)
      legend("topleft", bty = "n", legend = sprintf("λ = %.3f", lam))
    })
    output$enrich <- shiny::renderPlot({
      d <- res(); d <- d[!is.na(d$island), ]
      sig <- d$q < 0.05 & !is.na(d$q)
      if (sum(sig) < 5) { plot.new(); text(.5, .5, "Too few hits to test enrichment"); return() }
      ord <- c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea")
      lv <- intersect(ord, unique(d$island))
      or <- sapply(lv, function(k) {
        a <- sum(sig & d$island == k); b <- sum(sig) - a
        c1 <- sum(!sig & d$island == k); dd <- sum(!sig) - c1
        ((a + .5) / (b + .5)) / ((c1 + .5) / (dd + .5))
      })
      op <- par(mar = c(4, 8, 1, 2), las = 1); on.exit(par(op))
      cols <- ifelse(or >= 1, MP_PAL$hyper, MP_PAL$hypo)
      barplot(rev(log2(or)), horiz = TRUE, col = rev(cols), border = NA,
              names.arg = rev(lv),
              xlab = "log2 odds ratio vs tested background")
      abline(v = 0, col = "#444444")
    })
  })
}
