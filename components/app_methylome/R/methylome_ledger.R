##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Sample ledger: one TableModule of per-sample checks and one PlotModule of
## the beta distribution.

## ------------------------------------------------------------- ledger table --

methylome_table_ledger_ui <- function(id, title, info.text, caption, height, width) {
  ns <- shiny::NS(id)
  TableModuleUI(
    translate = FALSE,
    ns("tblmod"),
    title = title,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    label = "a"
  )
}

methylome_table_ledger_server <- function(id, pgx, scrollY = "22vh") {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      shiny::req(pgx())
      p <- pgx()
      ## Same Horvath the Epigenetic age screen shows, read off the stored fit
      ## rather than refitted. Without one, sample_qc() falls back to the
      ## wateRmelon single-clock path - fitting ten clocks to print one column
      ## is not worth it here (nor at upload, where there is no stored fit).
      st <- mp_stored_clocks(p)
      ages <- NULL
      if (!is.null(st)) {
        j <- grep("^horvath$", colnames(st$age), ignore.case = TRUE)
        if (length(j)) ages <- st$age[[j[1]]]
      }
      playbase.epigenetics::sample_qc(
        mp_beta(p), annot = p$genes, samples = p$samples, ages = ages
      )
    })

    render <- function(scrollY) {
      dt <- table_data()
      shiny::req(dt)
      DT::datatable(dt,
        class = "compact hover", rownames = TRUE,
        extensions = c("Buttons", "Scroller"), plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = 1),
        options = list(
          dom = "lfrtip", scroller = TRUE, scrollX = TRUE,
          scrollY = scrollY, scrollResize = TRUE, deferRender = TRUE
        )
      ) |>
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") |>
        DT::formatStyle("Verdict",
          color = DT::styleEqual(c("PASS", "CHECK"), c(MP_PAL$ok, MP_PAL$warn)),
          fontWeight = "bold") |>
        ## Only the mismatch is coloured - "ok" and "no record" are both
        ## uneventful, and red on a missing record would cry wolf.
        DT::formatStyle("Sex check",
          color = DT::styleEqual("MISMATCH", MP_PAL$bad, default = MP_PAL$grey))
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

## ------------------------------------------------------- beta distribution --

methylome_plot_betadist_ui <- function(id, title, caption, info.text, info.methods,
                                       height, width) {
  ns <- shiny::NS(id)
  PlotModuleUI(
    translate = FALSE,
    ns("pltmod"),
    plotlib = "base",
    title = title,
    caption = caption,
    info.text = info.text,
    info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width,
    label = "b"
  )
}

methylome_plot_betadist_server <- function(id, pgx, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx())
      X <- mp_beta(pgx())
      set.seed(1)
      X[sample(nrow(X), min(6000, nrow(X))), , drop = FALSE]
    })

    plot.RENDER <- function() {
      X <- plot_data()
      shiny::req(X)
      op <- graphics::par(mar = c(4, 4, 1, 1), las = 1)
      on.exit(graphics::par(op))
      plot(NA, xlim = c(0, 1), ylim = c(0, 4.2), xlab = "Beta value", ylab = "Density")
      for (j in seq_len(min(24, ncol(X)))) {
        graphics::lines(stats::density(X[, j], na.rm = TRUE, from = 0, to = 1),
                        col = grDevices::adjustcolor(MP_PAL$grey, 0.55), lwd = 1)
      }
      graphics::abline(v = c(0.2, 0.8), lty = 2, col = c(MP_PAL$hypo, MP_PAL$hyper))
      graphics::legend("topright", bty = "n", lty = 2,
                       col = c(MP_PAL$hypo, MP_PAL$hyper),
                       legend = c("hypo 0.2", "hyper 0.8"))
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "base",
      func = plot.RENDER,
      csvFunc = plot_data,
      renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot,
      res = c(90, 130),
      pdf.width = 8, pdf.height = 5,
      add.watermark = watermark
    )
  })
}
