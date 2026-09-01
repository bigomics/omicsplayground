##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Methylome character. Both plots read Relation_to_Island and
## genomic_location, which playbase has always produced and which nothing in
## the app read before this board.

## Colour encodes the beta value itself (hypo blue -> hyper red), so the bars
## read on the same scale as the ideogram.
mp_context_bars <- function(v) {
  ramp <- grDevices::colorRamp(c(MP_PAL$hypo, "#dcd6c8", MP_PAL$hyper))
  cols <- apply(ramp(pmin(pmax(v, 0), 1)), 1,
                function(z) grDevices::rgb(z[1], z[2], z[3], maxColorValue = 255))
  op <- graphics::par(mar = c(4, 8, 1, 2), las = 1); on.exit(graphics::par(op))
  bp <- graphics::barplot(rev(v), horiz = TRUE, xlim = c(0, 1), col = rev(cols),
                          border = NA, xlab = "mean beta")
  graphics::text(rev(v) - 0.04, bp, sprintf("%.2f", rev(v)),
                 col = "white", font = 2, cex = .9)
  graphics::abline(v = c(0.2, 0.8), lty = 3, col = MP_PAL$grey)
}

mp_context_profile <- function(pgx, what) {
  ctx <- mp_context_col(pgx, what)
  shiny::validate(shiny::need(!is.null(ctx),
    "This dataset carries no CpG context annotation."))
  ctx <- sub(";.*", "", ctx)
  ctx[is.na(ctx) | ctx == ""] <- "unannotated"
  mb <- rowMeans(mp_beta(pgx), na.rm = TRUE)
  v <- tapply(mb, ctx, mean, na.rm = TRUE)
  ord <- if (what == "island") {
    c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea")
  } else {
    c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR", "unannotated")
  }
  v[intersect(ord, names(v))]
}

methylome_plot_context_ui <- function(id, title, caption, info.text, info.methods,
                                      height, width, label = "a") {
  ns <- shiny::NS(id)
  PlotModuleUI(
    translate = FALSE,
    ns("pltmod"), plotlib = "base", title = title, caption = caption,
    info.text = info.text, info.methods = info.methods,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height, width = width, label = label
  )
}

methylome_plot_context_server <- function(id, pgx, what = c("island", "gene"),
                                          watermark = FALSE) {
  what <- match.arg(what)
  shiny::moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx())
      mp_context_profile(pgx(), what)
    })
    plot.RENDER <- function() {
      v <- plot_data()
      shiny::req(length(v) > 0)
      mp_context_bars(v)
    }
    PlotModuleServer(
      "pltmod", plotlib = "base", func = plot.RENDER,
      csvFunc = shiny::reactive(data.frame(context = names(plot_data()),
                                           mean_beta = as.numeric(plot_data()))),
      renderFunc = shiny::renderPlot, renderFunc2 = shiny::renderPlot,
      res = c(90, 130), pdf.width = 6, pdf.height = 5, add.watermark = watermark
    )
  })
}
