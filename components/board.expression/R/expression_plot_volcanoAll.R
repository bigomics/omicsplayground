##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
expression_plot_volcanoAll_ui <- function(id,
                                          label = "",
                                          height,
                                          width) {
  ns <- shiny::NS(id)

  info_text <- "Under the <strong>Volcano (all)</strong> tab, the platform simultaneously displays multiple volcano plots for genes across all contrasts. This provides users an overview of the statistics for all comparisons. By comparing multiple volcano plots, the user can immediately see which comparison is statistically weak or strong."

  PlotModuleUI(ns("pltmod"),
    title = "Volcano plots for all contrasts",
    label = label,
    plotlib = "grid",
    info.text = info_text,
    options = NULL,
    download.fmt = c("png", "pdf", "csv"),
    height = height,
    width = width
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
expression_plot_volcanoAll_server <- function(id,
                                              inputData,
                                              getAllContrasts,
                                              features,
                                              fdr,
                                              lfc,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    # reactive function listening for changes in input
    plot_data <- shiny::reactive({
      ngs <- inputData()
      features <- features()


      if (is.null(ngs)) {
        return(NULL)
      }
      ct <- getAllContrasts()
      F <- ct$F
      Q <- ct$Q

      ## comp = names(ngs$gx.meta$meta)
      comp <- names(F)
      if (length(comp) == 0) {
        return(NULL)
      }
      if (is.null(features)) {
        return(NULL)
      }

      fdr <- 1
      lfc <- 0
      fdr <- as.numeric(fdr())
      lfc <- as.numeric(lfc())

      sel.genes <- rownames(ngs$X)
      if (features != "<all>") {
        gset <- getGSETS(features)
        sel.genes <- unique(unlist(gset))
      }


      return(list(
        comp = comp,
        fdr = fdr,
        lfc = lfc,
        sel.genes = sel.genes,
        F = F,
        Q = Q
      ))
    })

    plot.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      shiny::withProgress(message = "rendering volcano plots ...", value = 0, {
        ## plot layout #####
        ng <- length(pd[["comp"]])
        nn <- c(2, max(ceiling(ng / 2), 5))
        ## if(ng>12) nn = c(3,8)
        par(mfrow = nn, mar = c(1, 1, 1, 1) * 0.2, mgp = c(2.6, 1, 0), oma = c(1, 1, 0, 0) * 2)
        nr <- 2
        nc <- ceiling(sqrt(ng))
        if (ng > 24) {
          nc <- max(ceiling(ng / 3), 6)
          nr <- 3
        } else if (TRUE && ng <= 4) {
          nc <- 4
          nr <- 1
        } else {
          nc <- max(ceiling(ng / 2), 6)
          nr <- 2
        }
        nr
        nc
        par(mfrow = c(nr, nc))

        ymax <- 15
        nlq <- -log10(1e-99 + unlist(pd[["Q"]]))
        ymax <- max(1.3, 1.2 * quantile(nlq, probs = 0.999, na.rm = TRUE)[1]) ## y-axis
        xmax <- max(1, 1.2 * quantile(abs(unlist(pd[["F"]])), probs = 0.999, na.rm = TRUE)[1]) ## x-axis


        plt <- list()
        i <- 1
        for (i in 1:length(pd[["comp"]])) {
          qval <- pd[["Q"]][[i]]
          fx <- pd[["F"]][[i]]
          fc.gene <- names(qval)
          is.sig <- (qval <= pd[["fdr"]] & abs(fx) >= pd[["lfc"]])
          sig.genes <- fc.gene[which(is.sig)]
          genes1 <- sig.genes[which(toupper(sig.genes) %in% toupper(pd[["sel.genes"]]))]
          genes2 <- head(genes1[order(-abs(fx[genes1]) * (-log10(qval[genes1])))], 10)
          xy <- data.frame(x = fx, y = -log10(qval))
          is.sig2 <- factor(is.sig, levels = c(FALSE, TRUE))

          plt[[i]] <- pgx.scatterPlotXY.GGPLOT(
            xy,
            title = pd[["comp"]][i], cex.title = 0.85,
            var = is.sig2, type = "factor",
            col = c("#bbbbbb", "#1e60bb"),
            legend.pos = "none", ## plotlib="ggplot",
            hilight = NULL, hilight2 = genes2,
            xlim = xmax * c(-1, 1), ylim = c(0, ymax),
            xlab = "difference  (log2FC)",
            ylab = "significance  (-log10q)",
            hilight.lwd = 0, hilight.col = "#1e60bb", hilight.cex = 1.5,
            cex = 0.45, cex.lab = 0.62
          )
          ## ggplot2::theme(legend.position='none')
          ## ggplot2::theme_bw(base_size=11)

          if (!interactive()) shiny::incProgress(1 / length(comp))
        }
      }) ## progress

      gridExtra::grid.arrange(grobs = plt, nrow = nr, ncol = nc)
    }

    # modal_plot.RENDER <- function() {
    #   plot.RENDER() %>%
    #     plotly::layout(
    #       ## showlegend = TRUE,
    #       font = list(
    #         size = 16
    #       )
    #     )
    # }


    PlotModuleServer(
      "pltmod",
      plotlib = "grid",
      func = plot.RENDER,
      # func2 = modal_plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(70, 90), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
