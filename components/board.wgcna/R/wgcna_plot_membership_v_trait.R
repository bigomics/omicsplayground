##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_membership_v_trait_ui <- function(
    id,
    label,
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    height = height,
    caption = caption,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_membership_v_trait_server <- function(id,
                                                 wgcna.compute,
                                                 selected_module,
                                                 watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    intraScatter.RENDER <- shiny::reactive({
      out <- wgcna.compute()

      MEs <- out$net$MEs
      ## rho1 <- cor(MEs, out$datExpr, use = "pairwise")
      rho1 <- cor(MEs, out$datExpr, use = "pairwise.complete.obs")
      rho1[is.na(rho1) | is.infinite(rho1)] <- 0

      ## rho2 <- cor(out$datTraits, out$datExpr, use = "pairwise")
      rho2 <- cor(out$datTraits, out$datExpr, use = "pairwise.complete.obs")
      rho2[is.na(rho2) | is.infinite(rho2)] <- 0

      ## rho3 <- cor(t(rho2), t(rho1), use = "pairwise")
      rho3 <- cor(t(rho2), t(rho1), use = "pairwise.complete.obs")
      rho3[is.na(rho3) | is.infinite(rho3)] <- 0

      k <- selected_module()
      in.mod <- colnames(rho1) %in% out$me.genes[[k]]
      table(in.mod)
      col1 <- c("grey60", out$me.colors[k])[1 + 1 * in.mod]

      ntop <- ifelse(nrow(rho3) >= 20, 20, 12)
      top.px <- head(order(-abs(rho3[, k])), ntop)
      if (ntop == 20) mfrow0 <- c(4, 5)
      if (ntop == 12) mfrow0 <- c(3, 4)

      par(mfrow = mfrow0, mar = c(4, 4, 2, 1), mgp = c(2.0, 0.8, 0))
      i <- top.px[1]
      for (i in top.px) {
        base::plot(rho1[k, ], rho2[i, ],
          pch = 20, cex = 0.7, col = col1,
          xlab = tspan("Module membership (eigengene cor)", js = FALSE),
          ylab = tspan("Gene significance (trait cor)", js = FALSE)
        )
        title(paste(k, "vs.", paste(rownames(rho2)[i])), cex = 1)
      }
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = intraScatter.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(85, 90),
      add.watermark = watermark
    )
  })
}
