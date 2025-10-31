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
  width
) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    shiny::checkboxInput(ns("show_cov"), "covariance", FALSE)
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    options = opts,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_plot_module_membership_server <- function(id,
                                                wgcna,
                                                pgx,
                                                selected_module,
                                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    render_plot <- function(ntop = 30) {
      res <- wgcna()
      module <- selected_module()
      shiny::req(!is.null(module) & module != "")

      rho <- res$stats[["moduleMembership"]][, module]
      rho[is.na(rho) | is.infinite(rho)] <- 0

      ylab0 <- "Eigengene correlation (rho)"
      if (input$show_cov) {
        sdx <- apply(res$datExpr, 2, sd, na.rm = TRUE)
        rho <- (rho * sdx**2)
        ylab0 <- "Eigengene covariance (cov)"
      }

      ## only ME genes
      sel <- which(names(rho) %in% res$me.genes[[module]])
      rho <- rho[sel]

      names(rho) <- playbase::probe2symbol(names(rho), pgx$genes, "gene_name", fill_na = TRUE)
      ## truncate long names, otherwise base plot errors
      names(rho) <- ifelse(nchar(names(rho)) > 15,
        paste0(substr(names(rho), 1, 12), "..."),
        names(rho)
      )

      if (min(rho, na.rm = TRUE) < 0) {
        ii <- unique(c(head(order(rho), ntop / 2), tail(order(rho), ntop / 2)))
      } else {
        ii <- tail(order(rho), ntop)
      }
      len <- max(nchar(names(rho)))
      bmar <- min(max(round(len / 2), 6), 12)
      par(mar = c(bmar, 4, 2, 0.1))
      barplot(sort(rho[ii], decreasing = TRUE),
        ylab = ylab0, las = 3,
        cex.names = 0.90, main = NULL
      )
      title(module, line = 1)
    }

    RENDER <- function() {
      render_plot(ntop = 20)
    }

    RENDER2 <- function() {
      render_plot(ntop = 50)
    }

    PlotModuleServer(
      "plot",
      func = RENDER,
      func2 = RENDER2,
      pdf.width = 8, pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )
  })
}
