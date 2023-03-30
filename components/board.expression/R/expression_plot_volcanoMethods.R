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
#' @param width
#'
#' @export
expression_plot_volcanoMethods_ui <- function(id,
                                              label = "",
                                              height,
                                              width) {
  ns <- shiny::NS(id)

  info_text <- "<b>Volcano plot for all statistical methods.</b> Simultaneous visualisation of volcano plots of genes by multiple differential expression methods for the selected contrast. Methods showing better statistical significance will show volcano plots with 'higher' wings. This provides a comparative overview of the relative statistical power between all methods."

  PlotModuleUI(ns("pltmod"),
    title = "Volcano plots for all methods",
    label = label,
    plotlib = "ggplot",
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
expression_plot_volcanoMethods_server <- function(id,
                                                  pgx,
                                                  comp, # input$gx_contrast
                                                  features, # input$gx_features
                                                  fdr, # input$gx_fdr
                                                  lfc, # input$gx_lfc
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    ## reactive function listening for changes in input
    plot_data <- shiny::reactive({
      comp <- comp()
      features <- features()

      if (is.null(comp)) {
        return(NULL)
      }
      shiny::req(pgx)
      if (is.null(features)) {
        return(NULL)
      }

      fdr <- as.numeric(fdr()) 
      lfc <- as.numeric(lfc()) 
      genes <- NULL

      gset <- getGSETS(features)
      sel.genes <- unique(unlist(gset))

      pd <- list(
          pgx = pgx,
          fdr = fdr,
          lfc = lfc,
          comp = comp,
          sel.genes = sel.genes
      )
        
      return(pd)
    })

    plot.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      ## meta tables
      mx <- pd[["pgx"]]$gx.meta$meta[[pd[["comp"]]]]
      fc <- unclass(mx$fc)
      ## pv = unclass(mx$p)
      qv <- unclass(mx$q)
      nlq <- -log10(1e-99 + qv)
      ymax <- max(3, 1.2 * quantile(nlq, probs = 0.999, na.rm = TRUE)[1]) ## y-axis
      xlim <- c(-1.1, 1.1) * max(abs(fc))
      xlim <- 1.3 * c(-1, 1) * quantile(abs(fc), probs = 0.999)
      fc.genes <- pd[["pgx"]]$genes[rownames(mx), "gene_name"]
      nplots <- min(24, ncol(qv))

      ## methods = names(pgx$gx.meta$output)
      methods <- colnames(pd[["pgx"]]$gx.meta$meta[[1]]$fc)
      nc <- 6
      par(mfrow = c(2, 6), mar = c(4, 4, 2, 2) * 0, oma = c(1, 1, 0, 0) * 2)
      if (nplots > 12) {
        nplots <- min(nplots, 24)
        par(mfrow = c(3, 8), mar = c(4, 4, 2, 2) * 0)
        nc <- 8
      }

      shiny::withProgress(message = "computing volcano plots ...", value = 0, {
        i <- 1
        for (i in 1:nplots) {
          fx <- fc[, i]
          ## pval = pv[,i]
          qval <- qv[, i]
          sig.genes <- fc.genes[which(qval <= pd[["fdr"]] & abs(fx) >= pd[["lfc"]])]
          ## genes1 = intersect(sig.genes, sel.genes)
          genes1 <- sig.genes[which(toupper(sig.genes) %in% toupper(pd[["sel.genes"]]))]
          gx.volcanoPlot.XY(
            x = fx, pv = qval, gene = fc.genes,
            render = "canvas", n = 5000, nlab = 5,
            xlim = xlim, ylim = c(0, ymax), axes = FALSE,
            use.fdr = TRUE, p.sig = pd[["fdr"]], lfc = pd[["lfc"]],
            ## main=comp[i],
            ## ma.plot=TRUE, use.rpkm=TRUE,
            cex = 0.6, lab.cex = 1.5, highlight = genes1
          )

          is.first <- (i %% nc == 1)
          last.row <- ((i - 1) %/% nc == (nplots - 1) %/% nc)
          is.first
          last.row
          if (is.first) axis(2, mgp = c(2, 0.7, 0), cex.axis = 0.8)
          if (last.row) axis(1, mgp = c(2, 0.7, 0), cex.axis = 0.8)
          graphics::box(lwd = 1, col = "black", lty = "solid")
          legend("top",
            legend = colnames(fc)[i], cex = 1.2,
            bg = "white", box.lty = 0, inset = c(0, 0.01),
            x.intersp = 0.1, y.intersp = 0.1
          )
          shiny::incProgress(1 / length(nplots))
        }
      })
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
      plotlib = "ggplot",
      func = plot.RENDER,
      # func2 = modal_plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(80, 170), ## resolution of plots
      pdf.width = 12, pdf.height = 5,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
