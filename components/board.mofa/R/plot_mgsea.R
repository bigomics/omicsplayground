##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_mgsea_ui <- function(
  id,
  title = "",
  info.text = "",
  info.methods,
  caption = "",
  label = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::radioButtons(ns("size.par"), "size by", c("p-value" = "p", "q-value" = "q"))
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    plotlib = "plotly",
    options = options,
    info.text = info.text,
    info.methods = info.methods,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}


mofa_plot_mgsea_server <- function(id,
                                   mgsea,
                                   input_k = reactive(1),
                                   select = reactive(NULL),
                                   watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      mgsea <- mgsea()
      validate(need(!is.null(mgsea), "missing GSEA data."))
      k <- input_k()
      shiny::req(k)

      df <- mgsea[[k]]
      cn <- colnames(df)
      types <- sub("^pval.", "", grep("^pval", cn, value = TRUE))
      types <- intersect(c("mx", "px", "gx", types), types)
      if (length(types) == 1) types <- rep(types, 2)

      ## scatter coordinates: enrichment (rho) of the two omics layers.
      ## A small jitter separates overlapping genesets.
      R <- df[, grep("^rho", cn), drop = FALSE]
      colnames(R) <- sub("^rho.", "", colnames(R))
      pos <- data.frame(
        x = R[, types[1]] + 2e-3 * rnorm(nrow(R)),
        y = R[, types[2]] + 2e-3 * rnorm(nrow(R))
      )
      rownames(pos) <- rownames(df)

      ## point size by significance (p- or q-value)
      qcomb <- if (input$size.par == "q") df$multi.q else df$multi.p
      cex1 <- 0.1 + (1 - qcomb)**2

      ## table-selected genesets are drawn as red filled dots, the rest as open
      ## grey circles
      selected <- select()
      if (is.null(selected) || length(selected) >= 100) selected <- character(0)

      d <- data.frame(
        pos,
        name = rownames(pos),
        size = 7 * cex1,
        score = df$multi.score,
        selected = rownames(pos) %in% selected
      )
      hover <- ~ paste0("<b>", name, "</b><br>multi.score: ", signif(score, 3))

      plotly::plot_ly() %>%
        plotly::add_markers(
          data = d[!d$selected, , drop = FALSE],
          x = ~x, y = ~y, key = ~name, text = hover, hoverinfo = "text",
          marker = list(
            size = ~size,
            color = "rgba(0,0,0,0)",
            line = list(color = "#999999AA", width = 1.2)
          )
        ) %>%
        plotly::add_markers(
          data = d[d$selected, , drop = FALSE],
          x = ~x, y = ~y, key = ~name, text = hover, hoverinfo = "text",
          marker = list(
            size = ~size,
            color = "#EE0000",
            line = list(color = "#EE0000", width = 1.2)
          )
        ) %>%
        plotly::layout(
          showlegend = FALSE,
          xaxis = list(title = paste(types[1], "enrichment (rho)")),
          yaxis = list(title = paste(types[2], "enrichment (rho)")),
          margin = list(l = 10, r = 10, b = 10, t = 10)
        )
    }

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = plot.RENDER,
      pdf.width = 8, pdf.height = 8,
      res = c(72, 110),
      add.watermark = watermark
    )
  })
}
