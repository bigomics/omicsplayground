##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_correlation_network_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height,
    width,
    ...) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg"),
    ...
  )
}

wgcna_plot_correlation_network_server <- function(id,
                                                  wgcna,
                                                  pgx,
                                                  selected_module,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    RENDER <- function() {
      out <- wgcna()
      k <- selected_module()
      shiny::req(k, out)
      rownames(out$stats$moduleMembership) <- playbase::probe2symbol(rownames(out$stats$moduleMembership), pgx$genes, "gene_name", fill_na = TRUE)
      colnames(out$datExpr) <- playbase::probe2symbol(colnames(out$datExpr), pgx$genes, "gene_name", fill_na = TRUE)
      playbase::wgcna.plotModuleHubGenes(
        out,
        modules = k, alpha = 0.5, setpar = TRUE
      )
    }

    PlotModuleServer(
      "plot",
      func = RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 120),
      add.watermark = watermark
    )
  })
}
