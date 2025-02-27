##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_module_significance_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height,
    width) {
  ns <- shiny::NS(id)
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_plot_module_significance_server <- function(id,
                                                  wgcna.compute,
                                                  selected_module,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    RENDER <- function() {
      res <- wgcna.compute()
      module <- selected_module()
      shiny::req(module)
      rho <- res$stats$moduleTraitCor[module,]
      ##rho <- rho[order(names(rho))]
      rho <- rho[order(rho, decreasing=TRUE)]
      par(mar=c(6,4,2,0.1))
      barplot( rho,
              ylab = "Trait correlation (rho)",
              main = "",
              width=1, las=3, names.arg='')
      dy <- 0.03*diff(range(rho))
      text(x = (-0.33 + 1:length(rho))*1.2,
           y = par("usr")[3] - dy,
           labels = names(rho),
           xpd = NA,
           srt = 45,
           adj = 0.965,
           cex = 0.9)

    }

    PlotModuleServer(
      "plot",
      func = RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 120),
      add.watermark = watermark
    )
  })
}
