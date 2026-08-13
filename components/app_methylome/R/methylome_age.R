##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Epigenetic age: clock agreement, DNAm age by phenotype, and a coverage
## table. A clock fitted on a fraction of its probes returns a confident wrong
## number, so anything with missing probes is withheld rather than caveated.

## Shared across the three modules on this screen.
mp_clock_set <- function(pgx) {
  X <- mp_beta(pgx)
  if (!"package:wateRmelon" %in% search()) suppressPackageStartupMessages(library(wateRmelon))
  raw <- as.data.frame(wateRmelon::agep(X, method = "all"))
  A <- mp_clocks(X)
  miss <- grep("n_missing$", colnames(raw), value = TRUE)
  nmiss <- if (length(miss)) vapply(raw[miss], max, numeric(1)) else stats::setNames(rep(0, ncol(A)), colnames(A))
  names(nmiss) <- sub("\\..*", "", miss)
  list(age = A, raw = raw, nmiss = nmiss, ok = intersect(colnames(A), names(nmiss)[nmiss == 0]))
}

## ------------------------------------------------------------ clock agreement --

methylome_plot_clocks_ui <- function(id, title, caption, info.text, info.methods,
                                     height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "a"
  )
}

methylome_plot_clocks_server <- function(id, pgx, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx())
      mp_clock_set(pgx())
    })
    plot.RENDER <- function() {
      cl <- plot_data()
      shiny::validate(shiny::need(length(cl$ok) >= 2,
        "Need two clocks with full probe coverage to compare."))
      x <- cl$age[[cl$ok[1]]]; y <- cl$age[[cl$ok[2]]]
      op <- graphics::par(mar = c(4, 4, 1, 1), las = 1); on.exit(graphics::par(op))
      plot(x, y, pch = 19, col = grDevices::adjustcolor(MP_PAL$hypo, .7), cex = 1.1,
           xlab = paste(cl$ok[1], "age"), ylab = paste(cl$ok[2], "age"))
      graphics::abline(0, 1, lty = 2, col = MP_PAL$grey)
      graphics::legend("topleft", bty = "n",
                       legend = sprintf("r = %.2f", stats::cor(x, y)))
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive(plot_data()$age),
      renderFunc = shiny::renderPlot, renderFunc2 = shiny::renderPlot,
      res = c(90, 130), pdf.width = 6, pdf.height = 6, add.watermark = watermark)
  })
}

## --------------------------------------------------------- DNAm age by group --

methylome_plot_agegroup_ui <- function(id, title, caption, info.text, info.methods,
                                       height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = "b"
  )
}

methylome_plot_agegroup_server <- function(id, pgx, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx())
      p <- pgx(); cl <- mp_clock_set(p)
      grp <- NULL; gname <- NULL
      for (cn in colnames(p$samples)) {
        v <- as.character(p$samples[[cn]])
        if (length(unique(v[!is.na(v) & v != ""])) == 2) { grp <- v; gname <- cn; break }
      }
      list(cl = cl, grp = grp, gname = gname)
    })
    plot.RENDER <- function() {
      d <- plot_data()
      shiny::validate(shiny::need(length(d$cl$ok) > 0, "No clock has full coverage."))
      shiny::validate(shiny::need(!is.null(d$grp), "No two-level phenotype on this dataset."))
      y <- d$cl$age[[d$cl$ok[1]]]
      op <- graphics::par(mar = c(4, 4, 1, 1), las = 1); on.exit(graphics::par(op))
      graphics::boxplot(y ~ factor(d$grp), col = c("#dce7f2", "#f2dce2"),
        xlab = d$gname, ylab = paste0("DNAm age (", d$cl$ok[1], ")"), border = "#444444")
      graphics::points(jitter(as.numeric(factor(d$grp)), 0.6), y, pch = 19,
        col = grDevices::adjustcolor(MP_PAL$grey, .6), cex = .8)
    }
    PlotModuleServer("pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive(plot_data()$cl$age),
      renderFunc = shiny::renderPlot, renderFunc2 = shiny::renderPlot,
      res = c(90, 130), pdf.width = 6, pdf.height = 6, add.watermark = watermark)
  })
}

## ------------------------------------------------------------ coverage table --

methylome_table_coverage_ui <- function(id, title, info.text, caption, height, width) {
  ns <- shiny::NS(id)
  TableModuleUI(ns("tblmod"), title = title, info.text = info.text,
                caption = caption, height = height, width = width, label = "c")
}

methylome_table_coverage_server <- function(id, pgx) {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      shiny::req(pgx())
      R <- mp_clock_set(pgx())$raw
      ages <- grep("\\.age$", colnames(R), value = TRUE)
      miss <- grep("n_missing$", colnames(R), value = TRUE)
      nm <- if (length(miss)) as.integer(sapply(R[miss], max)) else rep(0L, length(ages))
      data.frame(
        Clock = sub("\\..*", "", ages),
        `Median age` = ifelse(nm == 0, sprintf("%.1f", sapply(R[ages], stats::median)), "withheld"),
        `Probes missing` = nm,
        Status = ifelse(nm == 0, "full coverage",
                        paste0(nm, " probes absent — value not shown")),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })
    render <- function() {
      dt <- table_data(); shiny::req(dt)
      DT::datatable(dt, class = "compact hover", rownames = FALSE,
        selection = "none", options = list(dom = "t", paging = FALSE)) |>
        DT::formatStyle("Status",
          color = DT::styleEqual("full coverage", MP_PAL$ok, default = MP_PAL$warn))
    }
    TableModuleServer("tblmod", func = render, csvFunc = table_data, selector = "none")
  })
}
