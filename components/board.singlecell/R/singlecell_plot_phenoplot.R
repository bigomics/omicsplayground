##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Single cell plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
singlecell_plot_phenoplot_ui <- function(
    id,
    label = "",
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  phenoplot.opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("labelmode"), "Label:", c("groups", "legend"), inline = TRUE),
      "Select whether you want the group labels to be plotted inside the plots or in a seperate legend."
    )
  )

  PlotModuleUI(
    id = ns("plotmodule"),
    plotlib = "plotly",
    label = label,
    info.text = info.text,
    title = title,
    caption = caption,
    options = phenoplot.opts,
    download.fmt = c("png", "pdf"),
    height = height,
    width = width
  )
}

#' Single cell plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @export
singlecell_plot_phenoplot_server <- function(id,
                                             pgx,
                                             pfGetClusterPositions,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      #
      shiny::req(pgx$X)
      clust.pos <- pfGetClusterPositions()
      if (is.null(clust.pos)) {
        return(NULL)
      }

      pos <- pgx$tsne2d
      pos <- clust.pos
      sel <- rownames(pos)
      pheno <- colnames(pgx$Y)
      return(list(
        pos = pos,
        pgx = pgx,
        sel = sel,
        pheno = pheno
      ))
    })

    get_plots <- function(cex = 1) {
      #
      pd <- plot_data()
      shiny::req(pd)

      sel <- pd[["sel"]]
      pheno <- pd[["pheno"]]
      Y <- pd[["pgx"]]$samples
      pos <- pd[["pos"]]

      cex1 <- 1.2 * c(1.8, 1.3, 0.8, 0.5)[cut(nrow(pos), breaks = c(-1, 40, 200, 1000, 1e10))]
      cex1 <- cex1 * ifelse(length(pheno) > 6, 0.8, 1)
      cex1 <- cex1 * ifelse(length(pheno) > 12, 0.8, 1)

      ## is it a float/number???
      is.num <- function(y, fmin = 0.1) {
        suppressWarnings(numy <- as.numeric(as.character(y)))
        t1 <- !all(is.na(numy)) && is.numeric(numy)
        t2 <- (length(unique(y)) / length(y)) > fmin
        (t1 && t2)
      }

      plt <- list()
      i <- 6
      for (i in 1:min(20, length(pheno))) {
        px <- 4
        px <- pheno[i]
        y <- Y[sel, px]
        y[which(y %in% c(NA, "", " ", "NA", "na"))] <- NA
        if (sum(!is.na(y)) == 0) next

        if (playbase::is.num(y)) {
          klrpal <- colorRampPalette(c("grey90", "grey50", "red3"))(16)
          y <- rank(as.numeric(y))
          ny <- round(1 + 15 * (y - min(y)) / (max(y) - min(y)))
          klr0 <- klrpal[ny]
        } else {
          y <- factor(as.character(y))
          klrpal <- playdata::COLORS
          klrpal <- paste0(gplots::col2hex(klrpal), "99")
          klr0 <- klrpal[y]
        }

        p <- playbase::pgx.scatterPlotXY.PLOTLY(
          pos,
          var = y,
          type = "factor", ## always factor?
          col = klrpal,
          cex = 0.5 * cex1 * cex,
          xlab = NA,
          ylab = NA,
          xlim = 1.2 * range(pos[, 1]),
          ylim = 1.2 * range(pos[, 2]),
          axis = FALSE,
          title = tolower(pheno[i]),
          cex.title = 1.2 * cex,
          title.y = 1,
          cex.clust = 1.2 * cex,
          label.clusters = TRUE,
          legend = FALSE,
          gridcolor = "fff"
        ) %>%
          plotly::layout(
            #
            plot_bgcolor = "#f8f8f8",
            margin = list(0, 0, 0, 0)
          )

        plt[[i]] <- p
      }
      return(plt)
    }

    plotly.RENDER <- function() {
      pd <- plot_data()
      plt <- get_plots(cex = 0.9)
      shiny::req(plt)
      ## layout
      nr <- 2
      if (length(plt) > 4) nr <- 3
      if (length(plt) > 6) nr <- 4
      if (length(plt) > 12) nr <- 5
      if (length(plt) > 1) {
        fig <- plotly::subplot(
          plt,
          nrows = nr,
          #
          margin = c(0.01, 0.01, 0.01, 0.045)
        ) %>%
          plotly_default() %>%
          plotly::layout(
            margin = list(l = 10, r = 10, b = 10, t = 20) # lrbt
          )
      } else {
        fig <- plt[[1]] %>%
          plotly_default()
      }
      return(fig)
    }

    plotly_modal.RENDER <- function() {
      pd <- plot_data()
      plt <- get_plots(cex = 1.3)
      ## layout
      nc <- 2
      if (length(plt) > 4) nc <- 3
      if (length(plt) > 6) nc <- 4
      if (length(plt) > 12) nc <- 5
      nr <- ceiling(length(plt) / nc)
      fig <- plotly::subplot(
        plt,
        nrows = nr,
        margin = c(0.01, 0.01, 0.01, 0.045)
      ) %>%
        plotly_modal_default() %>%
        plotly::layout(
          margin = list(l = 20, r = 20, b = 20, t = 40) # lrbt
        )

      return(fig)
    }

    PlotModuleServer(
      "plotmodule",
      func = plotly.RENDER,
      func2 = plotly_modal.RENDER,
      plotlib = "plotly",
      res = c(85, 95),
      pdf.width = 6,
      pdf.height = 10,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
