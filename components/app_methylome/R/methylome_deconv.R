##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Reference-based cell-type deconvolution (Houseman projection).
##
## Cell composition is the dominant confounder in bulk-tissue methylation, so
## this screen exists as much to feed covariates into the EWAS model as to be
## read on its own. The estimator ships inside methylclock, which the clocks
## screen already depends on, so this costs no new dependency.

MP_DECONV_REFS <- c(
  "blood gse35069 complete", "blood gse35069", "blood gse35069 chen",
  "combined cord blood", "andrews and bakulski cord blood",
  "cord blood gse68456", "gervin and lyle cord blood",
  "guintivano dlpfc", "saliva gse48472"
)

mp_can_deconv <- function() {
  requireNamespace("methylclock", quietly = TRUE) &&
    "meffilEstimateCellCountsFromBetas" %in% getNamespaceExports("methylclock")
}

#' Estimate cell proportions. Returns a samples x celltype matrix, or NULL.
mp_cell_counts <- function(pgx, reference = MP_DECONV_REFS[1]) {
  if (!mp_can_deconv()) return(NULL)
  X <- mp_beta(pgx)
  cc <- tryCatch(
    methylclock::meffilEstimateCellCountsFromBetas(X, reference),
    error = function(e) {
      warning("[methylome] deconvolution failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(cc)) return(NULL)
  cc <- as.matrix(cc)
  rownames(cc) <- colnames(X)
  cc
}

## ---------------------------------------------------------------- settings --

## wrap = FALSE when the controls are being placed inside a merged screen's
## own tabSettings(), which is the only container the settings drawer reads;
## nesting a second one there renders nothing.
methylome_deconv_inputs <- function(id, wrap = TRUE) {
  ns <- shiny::NS(id)
  wrapper <- if (wrap) bigdash::tabSettings else shiny::tagList
  wrapper(
    withTooltip(
      shiny::selectInput(ns("deconv_ref"), "Reference panel:",
                         choices = MP_DECONV_REFS, selected = MP_DECONV_REFS[1]),
      "Which reference methylome panel to project onto. Pick one matching the tissue: the blood panels are for whole blood, the cord panels for neonatal samples, and there are saliva and brain (DLPFC) panels.",
      placement = "top"
    ),
    withTooltip(
      shiny::selectInput(ns("comp_pheno"), "Compare across:", choices = NULL),
      "Which sample column to compare the estimated proportions against.",
      placement = "top"
    ),
    shiny::div(
      style = "margin-top:6px;",
      shiny::actionButton(ns("run_deconv"), "Estimate composition",
                          class = "btn btn-primary btn-sm", width = "100%"),
      shiny::uiOutput(ns("deconv_stale"))
    )
  )
}

## ---------------------------------------------------------- stacked bars ----

methylome_plot_composition_ui <- function(id, title, caption, info.text, info.methods,
                                          height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    translate = FALSE,
    ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "a"
  )
}

methylome_plot_composition_server <- function(id, pgx, r.cells,
                                              r.pheno = shiny::reactive(NULL),
                                              watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      cc <- r.cells()
      shiny::validate(shiny::need(!is.null(cc),
        "Cell composition is not available. Press Estimate composition in the settings panel."))
      p <- pgx()
      ## Ordered by the phenotype the user picked, not by whichever sample
      ## column happened to come first - on a sheet led by slide or plate that
      ## silently bracketed the bars by batch under a phenotype label.
      ord <- seq_len(nrow(cc)); lab <- NULL; gname <- r.pheno()
      if (!is.null(gname) && nzchar(gname) && gname %in% colnames(p$samples)) {
        v <- as.character(p$samples[rownames(cc), gname])
        if (length(unique(v[!is.na(v) & v != ""])) >= 2) {
          ord <- order(v); lab <- v[ord]
        }
      }
      m <- t(cc[ord, , drop = FALSE])
      op <- graphics::par(mar = c(4.5, 4.2, 1, 7.5), las = 1, xpd = NA)
      on.exit(graphics::par(op))
      cols <- grDevices::hcl.colors(nrow(m), "Spectral")
      graphics::barplot(m, col = cols, border = NA, space = 0,
                        ylab = "estimated proportion", xaxt = "n",
                        xlab = if (is.null(lab)) "samples" else paste("samples, ordered by", gname))
      graphics::legend("topright", inset = c(-0.17, 0), bty = "n",
                       legend = rownames(m), fill = cols, cex = 0.85)
      if (!is.null(lab)) {
        b <- cumsum(rle(lab)$lengths); a <- c(0, utils::head(b, -1))
        graphics::segments(a, -0.03, b, -0.03, lwd = 3, col = "#697586")
        graphics::text((a + b) / 2, -0.075, rle(lab)$values, cex = 0.8, col = "#697586")
      }
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive(as.data.frame(r.cells())),
      renderFunc = shiny::renderPlot, renderFunc2 = shiny::renderPlot,
      res = c(90, 130), pdf.width = 10, pdf.height = 5, add.watermark = watermark)
  })
}

## ------------------------------------------------------ composition by group --

methylome_plot_compgroup_ui <- function(id, title, caption, info.text, info.methods,
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

methylome_plot_compgroup_server <- function(id, pgx, r.cells, r.pheno,
                                            watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      cc <- r.cells()
      shiny::validate(shiny::need(!is.null(cc),
        "Cell composition is not available. Press Estimate composition in the settings panel."))
      p <- pgx(); ph <- r.pheno()
      shiny::validate(shiny::need(!is.null(ph) && ph %in% colnames(p$samples),
        "No phenotype selected to compare composition against."))
      v <- as.character(p$samples[rownames(cc), ph])
      keep <- !is.na(v) & v != ""
      ## The picker is shared with the stacked bars, which order samples by the
      ## phenotype and are happy with a continuous one. This panel is not: it
      ## boxplots each group and tests between them, so age or BMI produced one
      ## box per distinct value at n = 1 and a Kruskal-Wallis across them. The
      ## old message claimed a categorical check that was never performed.
      lv <- unique(v[keep])
      shiny::validate(shiny::need(length(lv) >= 2 && length(lv) <= 8,
        sprintf("'%s' has %d distinct values, so it cannot group the samples for a comparison. Pick a categorical phenotype - two to eight groups - in the settings panel.",
                ph, length(lv))))
      g <- droplevels(factor(v[keep])); cc <- cc[keep, , drop = FALSE]
      k <- ncol(cc)
      op <- graphics::par(mfrow = c(ceiling(k / 4), min(4, k)),
                          mar = c(3.2, 3.6, 2, 0.6), las = 1, mgp = c(2.2, 0.6, 0))
      on.exit(graphics::par(op))
      for (i in seq_len(k)) {
        y <- cc[, i]
        graphics::boxplot(y ~ g, col = "#dce7f2", border = "#444444",
                          xlab = "", ylab = "proportion", outline = FALSE)
        graphics::points(jitter(as.numeric(g), 0.55), y, pch = 19, cex = 0.6,
                         col = grDevices::adjustcolor(MP_PAL$grey, 0.6))
        pv <- tryCatch(
          if (nlevels(g) == 2) stats::t.test(y ~ g)$p.value
          else stats::kruskal.test(y ~ g)$p.value, error = function(e) NA)
        graphics::title(main = sprintf("%s   p = %.3g", colnames(cc)[i], pv),
                        cex.main = 0.95, font.main = 1)
      }
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive(as.data.frame(r.cells())),
      renderFunc = shiny::renderPlot, renderFunc2 = shiny::renderPlot,
      res = c(90, 130), pdf.width = 9, pdf.height = 5, add.watermark = watermark)
  })
}

## --------------------------------------------------------------- table ------

methylome_table_composition_ui <- function(id, title, info.text, caption, height, width) {
  ns <- shiny::NS(id)
  TableModuleUI(
    translate = FALSE,ns("tblmod"), title = title, info.text = info.text,
                caption = caption, height = height, width = width, label = "c")
}

methylome_table_composition_server <- function(id, r.cells, scrollY = "22vh") {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      cc <- r.cells()
      shiny::validate(shiny::need(!is.null(cc),
        "Cell composition is not available. Press Estimate composition in the settings panel."))
      d <- as.data.frame(round(cc, 4))
      cbind(Sample = rownames(cc), d)
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
