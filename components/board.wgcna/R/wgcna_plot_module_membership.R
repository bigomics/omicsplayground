##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_module_membership_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height,
    width) {
  ns <- shiny::NS(id)

  eigenCorrelation_opts <- shiny::tagList(
    shiny::checkboxInput(ns("eigen_cov"), "covariance", FALSE)
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    options = eigenCorrelation_opts,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_module_membership_server <- function(id,
                                                wgcna,
                                                selected_module,
                                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    RENDER <- function() {
      res <- wgcna()
      module <- selected_module()
      
      rho <- res$stats[['moduleMembership']][,module]
      rho[is.na(rho) | is.infinite(rho)] <- 0
      
      ylab0 <- "Eigengene correlation (rho)"
      if (input$eigen_cov) {
        sdx <- apply(res$datExpr, 2, sd, na.rm = TRUE)
        rho <- (rho * sdx**2)
        ylab0 <- "Eigengene covariance (cov)"
      }

      ii <- unique(c(head(order(rho), ntop), tail(order(rho), ntop)))
      par(mar=c(6,4,2,0.1))
      barplot( sort(rho[ii]), ylab = ylab0, las = 3,
              cex.names = 0.90, main = NULL )
      title(module, line = 1)
    }

    PlotModuleServer(
      "plot",
      func = RENDER,
      pdf.width = 8, pdf.height = 5,
      res = c(80, 120),
      add.watermark = watermark
    )
  })
}
