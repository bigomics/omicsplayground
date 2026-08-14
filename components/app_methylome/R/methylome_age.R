##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Epigenetic age. Clocks come from methylclock when it is installed, which
## covers ten of them, and fall back to the five in wateRmelon otherwise, so
## the board still works on an image without the heavier dependency.
##
## A clock fitted on a fraction of its probes returns a confident wrong number
## (we have seen -13.5 years), so anything below the coverage floor is
## withheld rather than shown with a caveat.

MP_HAS_METHYLCLOCK <- function() requireNamespace("methylclock", quietly = TRUE)

## Families for the settings panel. Names are methylclock's own clock ids.
MP_CLOCK_FAMILIES <- list(
  "Chronological age" = c("Horvath", "skinHorvath", "Hannum", "BNN", "BLUP", "EN"),
  "Paediatric"        = c("PedBE", "Wu"),
  "Biological age"    = c("Levine", "TL")
)
MP_CLOCK_ALL <- unlist(MP_CLOCK_FAMILIES, use.names = FALSE)

## What each clock is, for the settings tooltip and the coverage table.
MP_CLOCK_NOTE <- c(
  Horvath     = "Multi-tissue, 353 CpGs (Horvath 2013)",
  skinHorvath = "Skin & blood, 391 CpGs (Horvath 2018)",
  Hannum      = "Whole blood, 71 CpGs (Hannum 2013)",
  BNN         = "Bayesian neural net on the Horvath/Hannum CpGs",
  BLUP        = "Best linear unbiased prediction, 319k CpGs",
  EN          = "Elastic net, 514 CpGs (Zhang 2019)",
  PedBE       = "Paediatric buccal epigenetic, 94 CpGs (McEwen 2019)",
  Wu          = "Paediatric blood, 111 CpGs (Wu 2019)",
  Levine      = "PhenoAge, mortality-trained, 513 CpGs (Levine 2018)",
  TL          = "Telomere length estimate, 140 CpGs (Lu 2019)"
)

## ---------------------------------------------------------------- settings --

methylome_age_inputs <- function(id) {
  ns <- shiny::NS(id)
  bigdash::tabSettings(
    withTooltip(
      shiny::checkboxGroupInput(
        ns("clock_families"), "Clock families:",
        choices = names(MP_CLOCK_FAMILIES),
        selected = c("Chronological age", "Biological age")
      ),
      "Toggle a whole family on or off. Individual clocks below stay editable afterwards.",
      placement = "top"
    ),
    withTooltip(
      shiny::checkboxGroupInput(
        ns("clocks"), "Clocks:",
        choices = MP_CLOCK_ALL,
        selected = c(MP_CLOCK_FAMILIES[["Chronological age"]],
                     MP_CLOCK_FAMILIES[["Biological age"]])
      ),
      "Individual clocks. This list is what is actually computed; the families above just set it.",
      placement = "top"
    ),
    withTooltip(
      shiny::sliderInput(
        ns("min_cov"), "Minimum probe coverage:",
        min = 0.5, max = 1, value = 0.8, step = 0.05
      ),
      "A clock with fewer than this fraction of its CpGs present is withheld rather than estimated from partial data.",
      placement = "top"
    ),
    ## Fitting ten clocks takes ~20s, so the settings above are staged and
    ## only applied on demand rather than recomputing on every tick.
    shiny::div(
      style = "margin-top: 6px;",
      shiny::actionButton(ns("recompute_clocks"), "Recompute clocks",
                          class = "btn btn-primary btn-sm", width = "100%"),
      shiny::uiOutput(ns("clock_stale"))
    ),

    ## ---- acceleration panel ----
    ## Which clock and which phenotype used to be whatever came first in the
    ## sample sheet, so on a cohort whose first categorical column is the slide
    ## or plate the headline read "age acceleration by phenotype" over a t-test
    ## against batch. Both are now chosen.
    shiny::tags$hr(style = "margin:12px 0;"),
    withTooltip(
      shiny::selectInput(ns("acc_clock"), "Acceleration clock:", choices = NULL),
      "Which clock the acceleration panel uses. Only clocks that passed the coverage floor are listed.",
      placement = "top"
    ),
    withTooltip(
      shiny::selectInput(ns("acc_pheno"), "Acceleration phenotype:", choices = NULL),
      "Which sample variable the acceleration is split by.",
      placement = "top"
    ),
    withTooltip(
      shiny::checkboxInput(ns("acc_intrinsic"), "Adjust for cell composition", FALSE),
      "Intrinsic acceleration: regress out the estimated cell proportions as well as chronological age, so a shift in blood composition is not reported as an accelerated methylome. Estimate the proportions first on the Cell composition screen.",
      placement = "top"
    )
  )
}

## ------------------------------------------------------------------ shared --

## Returns: age = data.frame of ages per usable clock, cov = per-clock
## coverage, acc = age acceleration when a chronological age column exists.
mp_clock_set <- function(pgx, clocks = MP_CLOCK_ALL, min_cov = 0.8) {
  X <- mp_beta(pgx)
  chron <- mp_chronological_age(pgx)

  if (MP_HAS_METHYLCLOCK()) {
    df <- data.frame(CpGName = rownames(X), X, check.names = FALSE)
    want <- intersect(clocks, MP_CLOCK_ALL)
    if (!length(want)) want <- MP_CLOCK_ALL
    args <- list(x = df, clocks = want, cell.count = FALSE, min.perc = min_cov)
    ## DNAmAge errors on age = NULL, so only pass it when we actually have one.
    if (!is.null(chron)) args$age <- chron
    res <- tryCatch(do.call(methylclock::DNAmAge, args),
                    error = function(e) {
                      warning("[methylome] DNAmAge failed: ", conditionMessage(e))
                      NULL
                    })
    if (!is.null(res)) {
      res <- as.data.frame(res)
      ages <- res[, intersect(want, colnames(res)), drop = FALSE]
      acc  <- res[, grep("ageAcc", colnames(res), value = TRUE), drop = FALSE]
      rownames(ages) <- colnames(X)
      if (ncol(acc)) rownames(acc) <- colnames(X)
      ok <- vapply(ages, function(v) any(!is.na(v)), logical(1))
      return(list(age = round(ages[, ok, drop = FALSE], 1),
                  all = ages, acc = acc, chron = chron,
                  source = "methylclock",
                  withheld = names(ok)[!ok]))
    }
  }

  ## Fallback: wateRmelon's five.
  A <- mp_clocks(X)
  if (is.null(A)) return(NULL)
  raw <- as.data.frame(wateRmelon::agep(X, method = "all"))
  miss <- grep("n_missing$", colnames(raw), value = TRUE)
  nm <- if (length(miss)) vapply(raw[miss], max, numeric(1)) else rep(0, ncol(A))
  names(nm) <- sub("\\..*", "", miss)
  ok <- intersect(colnames(A), names(nm)[nm == 0])
  list(age = A[, ok, drop = FALSE], all = A, acc = NULL, chron = chron,
       source = "wateRmelon", withheld = setdiff(colnames(A), ok))
}

## First numeric sample column that looks like a chronological age.
mp_chronological_age <- function(pgx) {
  s <- pgx$samples
  if (is.null(s)) return(NULL)
  cand <- grep("^age$|chron|age_", tolower(colnames(s)), value = FALSE)
  cand <- cand[!grepl("dnam|accel|group", tolower(colnames(s))[cand])]
  for (i in cand) {
    v <- suppressWarnings(as.numeric(as.character(s[[i]])))
    if (sum(!is.na(v)) >= 3 && stats::sd(v, na.rm = TRUE) > 0) return(v)
  }
  NULL
}

## ------------------------------------------------- DNAm vs chronological age --

methylome_plot_agecor_ui <- function(id, title, caption, info.text, info.methods,
                                     height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "a"
  )
}

methylome_plot_agecor_server <- function(id, r.clockset, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- r.clockset
    plot.RENDER <- function() {
      cl <- plot_data()
      shiny::validate(shiny::need(!is.null(cl) && ncol(cl$age) > 0,
        "No clock has enough probe coverage at this setting."))
      shiny::validate(shiny::need(!is.null(cl$chron),
        "This dataset has no chronological age column, so acceleration cannot be shown. DNAm age itself is on the other panels."))
      a <- cl$age; ch <- cl$chron
      k <- min(4, ncol(a))
      op <- graphics::par(mfrow = c(ifelse(k > 2, 2, 1), ifelse(k > 1, 2, 1)),
                          mar = c(3.8, 3.8, 2, 0.8), las = 1, mgp = c(2.3, 0.6, 0))
      on.exit(graphics::par(op))
      for (i in seq_len(k)) {
        y <- a[[i]]
        graphics::plot(ch, y, pch = 19, cex = 1, col = grDevices::adjustcolor(MP_PAL$hypo, .65),
             xlab = "chronological age", ylab = paste(colnames(a)[i], "age"))
        graphics::abline(0, 1, lty = 2, col = MP_PAL$grey)
        good <- !is.na(ch) & !is.na(y)
        if (sum(good) > 2) {
          graphics::abline(stats::lm(y[good] ~ ch[good]), col = MP_PAL$hyper, lwd = 1.6)
          graphics::title(main = sprintf("%s   r = %.2f", colnames(a)[i],
                                         stats::cor(ch[good], y[good])),
                          cex.main = 1, font.main = 1)
        }
      }
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive({
        cl <- plot_data(); cbind(chronological_age = cl$chron, cl$age)
      }),
      renderFunc = shiny::renderPlot, renderFunc2 = shiny::renderPlot,
      res = c(90, 130), pdf.width = 8, pdf.height = 7, add.watermark = watermark)
  })
}

## ------------------------------------------------------------ clock agreement --

methylome_plot_clocks_ui <- function(id, title, caption, info.text, info.methods,
                                     height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "b"
  )
}

methylome_plot_clocks_server <- function(id, r.clockset, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- r.clockset
    plot.RENDER <- function() {
      cl <- plot_data()
      shiny::validate(shiny::need(!is.null(cl) && ncol(cl$age) >= 2,
        "Select at least two clocks with enough coverage to compare them."))
      a <- cl$age
      op <- graphics::par(mar = c(4, 4, 1, 1), las = 1); on.exit(graphics::par(op))
      graphics::pairs(a, pch = 19, cex = 0.55,
                      col = grDevices::adjustcolor(MP_PAL$hypo, 0.6),
                      gap = 0.25, oma = c(2.5, 2.5, 1.5, 1.5),
                      lower.panel = function(x, y, ...) {
                        graphics::points(x, y, pch = 19, cex = 0.55,
                                         col = grDevices::adjustcolor(MP_PAL$hypo, 0.6))
                        graphics::abline(0, 1, lty = 3, col = MP_PAL$grey)
                      },
                      upper.panel = function(x, y, ...) {
                        graphics::text(mean(range(x, na.rm = TRUE)),
                                       mean(range(y, na.rm = TRUE)),
                                       sprintf("%.2f", stats::cor(x, y, use = "complete.obs")),
                                       cex = 1.1, col = MP_PAL$hyper)
                      })
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive(plot_data()$age),
      renderFunc = shiny::renderPlot, renderFunc2 = shiny::renderPlot,
      res = c(90, 130), pdf.width = 7, pdf.height = 7, add.watermark = watermark)
  })
}

## --------------------------------------------------- acceleration by phenotype --

methylome_plot_agegroup_ui <- function(id, title, caption, info.text, info.methods,
                                       height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "c"
  )
}

methylome_plot_agegroup_server <- function(id, pgx, r.clockset,
                                           r.pheno = shiny::reactive(NULL),
                                           r.clock = shiny::reactive(NULL),
                                           r.cells = shiny::reactive(NULL),
                                           r.intrinsic = shiny::reactive(FALSE),
                                           watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx())
      p <- pgx(); cl <- r.clockset()
      shiny::validate(shiny::need(!is.null(cl) && ncol(cl$age) > 0,
        "No clock has enough probe coverage at this setting."))
      gname <- r.pheno()
      shiny::validate(shiny::need(!is.null(gname) && nzchar(gname) &&
                                    gname %in% colnames(p$samples),
        "Pick a phenotype for the acceleration panel in the settings."))
      k <- r.clock()
      if (is.null(k) || !k %in% colnames(cl$age)) k <- colnames(cl$age)[1]
      ## Intrinsic acceleration needs the proportions; asking for it without
      ## having estimated them is a prompt, not a silent fall back to raw.
      cf <- NULL
      if (isTRUE(r.intrinsic())) {
        cf <- r.cells()
        shiny::validate(shiny::need(!is.null(cf) && nrow(cf) > 0,
          "Estimate cell composition first on the Cell composition screen, or untick the adjustment."))
      }
      list(cl = cl, grp = as.character(p$samples[[gname]]), gname = gname,
           clock = k, cf = cf)
    })
    plot.RENDER <- function() {
      d <- plot_data()
      cl <- d$cl; k <- d$clock
      age <- cl$age[[k]]
      ## Acceleration is the residual of DNAm age on chronological age. Without
      ## a chronological age we can only show the raw DNAm age.
      if (!is.null(cl$chron)) {
        df <- data.frame(.a = age, .c = cl$chron, row.names = rownames(cl$age))
        adj <- ""
        if (!is.null(d$cf)) {
          ## Match on sample name, not position - the proportions come from a
          ## separate run and need not be in the clock table's order.
          cf <- d$cf[match(rownames(df), rownames(d$cf)), , drop = FALSE]
          ## One proportion is redundant with the rest, so drop the largest and
          ## keep the design full rank - same rule as the EWAS model.
          cf <- cf[, order(-colMeans(cf, na.rm = TRUE)), drop = FALSE][, -1, drop = FALSE]
          for (cn in colnames(cf)) df[[cn]] <- as.numeric(cf[, cn])
          adj <- " intrinsic"
        }
        y <- stats::resid(stats::lm(.a ~ ., data = df, na.action = stats::na.exclude))
        ylab <- paste0(k, adj, " age acceleration (years)")
      } else {
        y <- age
        ylab <- paste0(k, " DNAm age (years)")
      }
      keep <- !is.na(d$grp) & d$grp != "" & !is.na(y)
      shiny::validate(shiny::need(sum(keep) > 2, "Too few samples with this phenotype."))
      g <- factor(d$grp[keep]); y <- y[keep]
      shiny::validate(shiny::need(nlevels(g) >= 2,
        "This phenotype has a single level across the cohort."))
      op <- graphics::par(mar = c(4.2, 4.4, 1, 1), las = 1); on.exit(graphics::par(op))
      graphics::boxplot(y ~ g, col = "#dce7f2", border = "#444444",
                        xlab = d$gname, ylab = ylab, outline = FALSE)
      graphics::points(jitter(as.numeric(g), 0.6), y, pch = 19,
                       col = grDevices::adjustcolor(MP_PAL$grey, .6), cex = .8)
      if (!is.null(cl$chron)) graphics::abline(h = 0, lty = 3, col = MP_PAL$grey)
      lab <- if (nlevels(g) == 2) {
        sprintf("t-test p = %.3g", stats::t.test(y ~ g)$p.value)
      } else {
        sprintf("Kruskal-Wallis p = %.3g", stats::kruskal.test(y ~ g)$p.value)
      }
      graphics::mtext(lab, side = 3, adj = 1, cex = 0.8, col = MP_PAL$grey)
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive({
        d <- plot_data()
        cbind(sample = rownames(d$cl$age), phenotype = d$grp, d$cl$age)
      }),
      renderFunc = shiny::renderPlot, renderFunc2 = shiny::renderPlot,
      res = c(90, 130), pdf.width = 6, pdf.height = 6, add.watermark = watermark)
  })
}

## ------------------------------------------------------------ coverage table --

methylome_table_coverage_ui <- function(id, title, info.text, caption, height, width) {
  ns <- shiny::NS(id)
  TableModuleUI(ns("tblmod"), title = title, info.text = info.text,
                caption = caption, height = height, width = width, label = "d")
}

methylome_table_coverage_server <- function(id, r.clockset) {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      cl <- r.clockset()
      shiny::validate(shiny::need(!is.null(cl), "Clocks could not be computed."))
      nm <- colnames(cl$all)
      usable <- colnames(cl$age)
      data.frame(
        Clock = nm,
        What = ifelse(nm %in% names(MP_CLOCK_NOTE), MP_CLOCK_NOTE[nm], ""),
        `Median age` = vapply(nm, function(k) {
          v <- cl$all[[k]]
          if (k %in% usable && any(!is.na(v))) sprintf("%.1f", stats::median(v, na.rm = TRUE)) else "withheld"
        }, character(1)),
        Status = ifelse(nm %in% usable, "usable",
                        "below coverage floor - not estimated"),
        check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE
      )
    })
    render <- function() {
      dt <- table_data(); shiny::req(dt)
      DT::datatable(dt, class = "compact hover", rownames = FALSE,
        selection = "none", options = list(dom = "t", paging = FALSE, scrollX = TRUE)) |>
        DT::formatStyle("Status",
          color = DT::styleEqual("usable", MP_PAL$ok, default = MP_PAL$warn))
    }
    TableModuleServer("tblmod", func = render, csvFunc = table_data, selector = "none")
  })
}
