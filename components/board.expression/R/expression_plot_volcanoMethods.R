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

    render_plots <- function(base_size=11) {
      pd <- plot_data()
      shiny::req(pd)

      ## meta tables
      comp <- pd[["comp"]]
      mx <- pd[["pgx"]]$gx.meta$meta[[comp]]
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
      plt <- list()
      
      shiny::withProgress(message = "computing volcano plots ...", value = 0, {
        i <- 1
        for (i in 1:nplots) {
          fx <- fc[, i]
          ## pval = pv[,i]
          qval <- qv[, i]
          sig.genes <- fc.genes[which(qval <= pd[["fdr"]] & abs(fx) >= pd[["lfc"]])]
          ## genes1 = intersect(sig.genes, sel.genes)
          genes2 <- sig.genes[which(toupper(sig.genes) %in% toupper(pd[["sel.genes"]]))]

          ## gx.volcanoPlot.XY(
          ##   x = fx, pv = qval, gene = fc.genes,
          ##   render = "canvas", n = 5000, nlab = 5,
          ##   xlim = xlim, ylim = c(0, ymax), axes = FALSE,
          ##   use.fdr = TRUE, p.sig = pd[["fdr"]], lfc = pd[["lfc"]],
          ##   ## main=comp[i],
          ##   ## ma.plot=TRUE, use.rpkm=TRUE,
          ##   cex = 0.6, lab.cex = 1.5, highlight = genes1
          ## )

          xy <- data.frame(x = fx, y = -log10(qval))
          is.sig1 <- fc.genes %in% sig.genes
          is.sig2 <- fc.genes %in% genes2
          
          plt[[i]] <- pgx.scatterPlotXY.GGPLOT(
            xy,
            title = methods[i],
            cex.title = 0.85,
            var = is.sig1,
            type = "factor",
            col = c("#bbbbbb", "#1e60bb"),
            legend.pos = "none", ## plotlib="ggplot",
            hilight = NULL,
            hilight2 = genes2,
            xlim = xlim,
            ylim = c(0,ymax),
            xlab = "difference  (log2FC)",
            ylab = "significance  (-log10q)",
            hilight.lwd = 0,
            hilight.col = "#1e60bb",
            hilight.cex = 1.5,
            cex = 0.45,
            cex.lab = 0.62,
            base_size = base_size
          )

          if (!interactive()) shiny::incProgress(1 / length(pd[["comp"]]))
        }
      })

      return(plt)
    }
      
    plot.RENDER <- function() {      
      plt <- render_plots(base_size=12)
      nplots <- length(plt)
      
      ## layout
      nr = 1
      nc = 5
      if(nplots > 5) {
        nr = 2
        nc = 6
      }
      if (nplots > 12) {
        nr = 3
        nc = 8
      }
      gridExtra::grid.arrange(grobs = plt, nrow = nr, ncol = nc)
    }

    modal_plot.RENDER <- function() {      
      plt <- render_plots(base_size=18)
      nplots <- length(plt)
      
      ## layout
      nr = 1
      nc = 2
      if(nplots > 3) {
        nr = 2
        nc = 4
      }
      if (nplots > 8) {
        nr = 3
        nc = 6
      }
      gridExtra::grid.arrange(grobs = plt, nrow = nr, ncol = nc)
    }
    
    PlotModuleServer(
      "pltmod",
##      plotlib = "ggplot",
      plotlib = "grid",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      ## csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(80, 170), ## resolution of plots
      pdf.width = 12, pdf.height = 5,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
