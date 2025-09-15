##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_snf_ui <- function(
    id,
    title = "",
    info.text = "",
    info.references,
    info.methods,
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::radioButtons(
      ns("type"),
      label = "Plot type:",
      choices = c("Affinity matrix", "t-SNE"),
      selected = "Affinity matrix"
    ),
    shiny::conditionalPanel(
      condition = "input.type == 't-SNE'",
      ns = ns,
      shiny::radioButtons(
        ns("tsne_colorby"),
        label = "Color by variable:",
        choices = c("Var1", "Var2"),
        selected = "Var1"
      )
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

mofa_plot_snf_server <- function(id,
                                 mofa,
                                 pgx_samples,
                                 watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    observeEvent(mofa(), {
      res <- mofa()
      phenos <- colnames(res$samples)
      updateRadioButtons(
        session,
        "tsne_colorby",
        choices = phenos,
        selected = phenos[1]
      )
    })

    plot.RENDER <- function() {
      res <- mofa()
      pgx_samples <- pgx_samples()
      snf <- res$snf
      validate(need(!is.null(res), "missing MOFA data."))
      type <- input$type

      if (type == "Affinity matrix") {
        par(mfrow = c(2, 2), mar = c(6, 1, 2, 8))
        ndim <- ncol(snf$affinityMatrix[[1]])
        if (ndim > 20) par(mar = c(3, 1, 2, 4))
        nmat <- length(snf$affinityMatrix) + 1
        if (nmat > 4) par(mfrow = c(3, 3))
        playbase::snf.plot_affinity(snf, k = 0.5, par = FALSE)
      }

      if (type == "t-SNE") {
        ph <- input$tsne_colorby
        cc <- factor(pgx_samples[, ph])
        par(mfrow = c(2, 2), mar = c(5, 5, 2, 1))
        i <- 1
        for (i in 1:length(snf$posx)) {
          plot(snf$posx[[i]],
            col = cc, pch = 19, cex = 1,
            xlab = "tSNE 1", ylab = "tSNE 2", las = 1
          )
          title(names(snf$posx)[i], cex.main = 1.4)
        }
      }
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 100),
      add.watermark = watermark
    )
  })
}
