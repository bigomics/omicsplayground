##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_heatmap_membership_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height,
    width) {
  ns <- shiny::NS(id)

  intraHeatmap_opts <- shiny::tagList(
    shiny::checkboxInput(ns("eigen_cov"), "covariance", FALSE)
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    options = intraHeatmap_opts,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_heatmap_membership_server <- function(id,
                                                 wgcna.compute,
                                                 watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    intraHeatmap.RENDER <- shiny::reactive({
      out <- wgcna.compute()

      MEs <- out$net$MEs
      if (input$eigen_cov) {
        ## rho1 <- cov(MEs, out$datExpr, use = "pairwise")
        rho1 <- cov(MEs, out$datExpr, use = "pairwise.complete.obs")
        rho1[is.na(rho1) | is.infinite(rho1)] <- 0

        ## rho2 <- cov(out$datTraits, out$datExpr, use = "pairwise")
        rho2 <- cov(out$datTraits, out$datExpr, use = "pairwise.complete.obs")
        rho2[is.na(rho2) | is.infinite(rho2)] <- 0

        ## rho3 <- cov(t(rho2), t(rho1), use = "pairwise") ## cor or cov??
        rho3 <- cov(t(rho2), t(rho1), use = "pairwise.complete.obs")
        rho3[is.na(rho3) | is.infinite(rho3)] <- 0
      } else {
        ## rho1 <- cor(MEs, out$datExpr, use = "pairwise")
        rho1 <- cor(MEs, out$datExpr, use = "pairwise.complete.obs")
        rho1[is.na(rho1) | is.infinite(rho1)] <- 0

        ## rho2 <- cor(out$datTraits, out$datExpr, use = "pairwise")
        rho2 <- cor(out$datTraits, out$datExpr, use = "pairwise.complete.obs")
        rho2[is.na(rho2) | is.infinite(rho2)] <- 0

        ## rho3 <- cor(t(rho2), t(rho1), use = "pairwise")
        rho3 <- cor(t(rho2), t(rho1), use = "pairwise.complete.obs")
        rho3[is.na(rho3) | is.infinite(rho3)] <- 0
      }

      playbase::gx.heatmap(rho3,
        nmax = 50, mar = c(5, 10),
        keysize = 0.5, scale = "none", key = FALSE
      )
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = intraHeatmap.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(85, 100),
      add.watermark = watermark
    )
  })
}
