##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
expression_plot_volcanoAll_ui <- function(
  id,
  title,
  caption,
  info.text,
  label = "",
  height,
  width) {
  ns <- shiny::NS(id)

  PlotModuleUI(ns("pltmod"),
    title = title,
    label = label,
    plotlib = "grid",
    info.text = info.text,
    caption = caption,
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
                                              pgx,
                                              getAllContrasts,
                                              features,
                                              fdr,
                                              lfc,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    ## reactive function listening for changes in input
    plot_data <- shiny::reactive({

      if (is.null(pgx)) {
        return(NULL)
      }
      
      features <- features()
      ct <- getAllContrasts()
      F <- ct$F
      Q <- ct$Q

      
      ## comp = names(pgx$gx.meta$meta)
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

      sel.genes <- rownames(pgx$X)
      if (features != "<all>") {
        gset <- getGSETS(features)
        sel.genes <- unique(unlist(gset))
      }

      ## combined matrix for output
      matF <- do.call(cbind,F)
      colnames(matF) <- paste0("fc.",names(F))
      matQ <- do.call(cbind,Q)
      colnames(matQ) <- paste0("q.",names(Q))
      FQ <- cbind(matF, matQ)
      
      pd <- list(
        FQ = FQ,   ## Remember: the first element is returned as downloadable CSV
        comp = comp,
        fdr = fdr,
        lfc = lfc,
        sel.genes = sel.genes,
        F = F,
        Q = Q
      )

      return(pd)
    })

    render_plots <- function(cex=0.45, base_size=11) {
      pd <- plot_data()
      shiny::req(pd)

      ymax <- 15
      nlq <- -log10(1e-99 + unlist(pd[["Q"]]))
      ymax <- max(1.3, 1.2 * quantile(nlq, probs = 0.999, na.rm = TRUE)[1]) ## y-axis
      xmax <- max(1, 1.2 * quantile(abs(unlist(pd[["F"]])), probs = 0.999, na.rm = TRUE)[1]) ## x-axis
      
      plt <- list()

      shiny::withProgress(message = "rendering volcano plots ...", value = 0, {
        i <- 1
        for (i in 1:length(pd[["comp"]])) {

          qval <- pd[["Q"]][[i]]
          fx <- pd[["F"]][[i]]
          fc.genes <- names(qval)
          is.sig <- (qval <= pd[["fdr"]] & abs(fx) >= pd[["lfc"]])
          sig.genes <- fc.genes[which(is.sig)]
          genes1 <- sig.genes[which(toupper(sig.genes) %in% toupper(pd[["sel.genes"]]))]
          genes2 <- head(genes1[order(-abs(fx[genes1]) * (-log10(qval[genes1])))], 10)
          xy <- data.frame(x = fx, y = -log10(qval))
          is.sig1 <- factor(is.sig, levels = c(FALSE, TRUE))

          plt[[i]] <- playbase::pgx.scatterPlotXY.GGPLOT(
            xy,
            title = pd[["comp"]][i],
            cex.title = 0.85,
            var = is.sig1,
            type = "factor",
            col = c("#bbbbbb", "#1e60bb"),
            legend.pos = "none", ## plotlib="ggplot",
            hilight = NULL,
            hilight2 = genes2,
            xlim = xmax * c(-1, 1),
            ylim = c(0, ymax),
            xlab = "difference  (log2FC)",
            ylab = "significance  (-log10q)",
            hilight.lwd = 0,
            hilight.col = "#1e60bb",
            hilight.cex = 1.5,
            cex = cex,
            cex.lab = 1.8*cex,
            base_size = base_size
          ) + ggplot2::theme_bw(base_size = base_size)
          
          if (!interactive()) shiny::incProgress(1 / length(pd[["comp"]]))
        }
      }) ## progress

      return(plt)
    }


    plot.RENDER <- function() {
      plt <- render_plots(cex=0.5, base_size=11)
      nplots <- length(plt)
      ## plot layout #####
      ## layout
      nr = 1
      nc = max(4,nplots)
      if(nplots > 6) {
        nr = 2
        nc = ceiling(nplots/nr)
      }
      if (nplots > 12) {
        nr = 3
        nc = ceiling(nplots/nr)
      }
      gridExtra::grid.arrange(grobs = plt, nrow = nr, ncol = nc)
    }


    modal_plot.RENDER <- function() {      
      plt <- render_plots(cex=0.9, base_size=15)
      nplots <- length(plt)
      ## layout
      nr = 1
      nc = max(2,nplots)
      if(nplots > 3) {
        nr = 2
        nc = ceiling(nplots/nr)
      }
      if (nplots > 8) {
        nr = 3
        nc = ceiling(nplots/nr)
      }
      gridExtra::grid.arrange(grobs = plt, nrow = nr, ncol = nc)
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "grid",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      ## csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(70, 90), ## resolution of plots
      pdf.width = 12, pdf.height = 5,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
