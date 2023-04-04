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
singlecell_plot_icpplot_ui <- function(id,
                                       label = "",
                                       height,
                                       width,
                                       parent) {
  ns <- shiny::NS(id)

  icp.opts <- shiny::tagList(
    withTooltip(shiny::selectInput(parent("refset"), "Reference:", choices = NULL),
      "Select a reference dataset for the cell type prediction.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(shiny::selectInput(parent("dcmethod"), "Method:", choices = NULL),
      "Choose a method for the cell type prediction.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(
      shiny::radioButtons(parent("sortby"), "Sort by:",
        choices = c("probability", "name"), inline = TRUE
      ),
      "Sort by name or probability.",
      placement = "top",
      options = list(container = "body")
    ),
    withTooltip(
      shiny::radioButtons(parent("layout"), "Layout:",
        choices = c("4x4", "6x6"),
        ## selected="6x6",
        inline = TRUE
      ),
      "Choose layout.",
      placement = "top", options = list(container = "body")
    )
  )

  icp_info <- "<strong>Cell type profiling</strong> infers the type of cells using computational deconvolution methods and reference datasets from the literature. Currently, we have implemented a total of 8 methods and 9 reference datasets to predict immune cell types (4 datasets), tissue types (2 datasets), cell lines (2 datasets) and cancer types (1 dataset). However, we plan to expand the collection of methods and databases and to infer other cell types."

  PlotModuleUI(
    id = ns("plot"),
    ##    plotlib = "plotly",
    plotlib = "ggplot",      
    label = label,
    info.text = icp_info,
    title = "Cell type profiling",
    options = icp.opts,
    download.fmt = c("png", "pdf", "csv"),
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
singlecell_plot_icpplot_server <- function(id,
                                           pgx,
                                           pfGetClusterPositions,
                                           method, # input$dcmethod
                                           refset, # input$refset
                                           lyo, # input$layout
                                           sortby, # input$sortby
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({

      shiny::req(pgx)
        
      method <- "meta"
      refset <- "LM22"
      method <- method() # input$dcmethod
      if (is.null(method)) {
        return(NULL)
      }
      refset <- refset() #
      lyo <- lyo()
      sortby <- sortby()

      if (!("deconv" %in% names(pgx))) {
        return(NULL)
      }
      results <- pgx$deconv[[refset]][[method]]
      ## threshold everything (because DCQ can be negative!!!)
      results <- pmax(results, 0)

      clust.pos <- pfGetClusterPositions()

      # alertDataLoaded(session,pgx)

      if (is.null(clust.pos)) {
        return(NULL)
      }
      pos <- pgx$tsne2d
      pos <- clust.pos
      score <- pgx$deconv[[1]][["meta"]]
      score <- results
      if (is.null(score) || length(score) == 0) {
        return(NULL)
      }

      ## normalize
      score <- score[rownames(pos), , drop = FALSE]
      score[is.na(score)] <- 0
      score <- pmax(score, 0)
      score <- score / (1e-20 + rowSums(score))
      score <- tanh(score / mean(abs(score)))
      score <- score / max(score, na.rm = TRUE)
      summary(as.vector(score))

      ## take top10 features
      jj.top <- unique(as.vector(apply(score, 1, function(x) head(order(-x), 10))))
      score <- score[, jj.top]
      score <- score[, order(-colMeans(score**2))]
      score <- score[, 1:min(50, ncol(score))]
      ii <- hclust(dist(score))$order
      jj <- hclust(dist(t(score)))$order
      score <- score[ii, jj]

      score0 <- score
      pos <- pos[rownames(score), ]
      b0 <- 1 + 0.85 * pmax(30 - ncol(score), 0)

      pd <- list(
          score = score,
          pos = pos,
          lyo = lyo,
          refset = refset,
          sortby = sortby
      )
      return(pd)
    })
    
    base.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)
      
      cex1 <- 1.2
      cex.bin <- cut(nrow(pd[["pos"]]), breaks = c(-1, 40, 200, 1000, 1e10))
      cex1 <- 0.9 * c(2.2, 1.1, 0.6, 0.3)[cex.bin]
      klrpal <- colorRampPalette(c("grey90", "grey50", "red3"))(16)
      klrpal <- paste0(gplots::col2hex(klrpal), "66")

      ntop <- 25
      par(mfrow = c(5, 5), mar = c(0.2, 0.2, 1.8, 0.2), oma = c(1, 1, 1, 1) * 0.8)
      par(mfrow = c(5, 5), mar = c(0, 0.2, 0.5, 0.2), oma = c(1, 1, 6, 1) * 0.5)
      if (ncol(pd[["score"]]) > 25) par(mfrow = c(6, 6), mar = c(0, 0.2, 0.5, 0.2) * 0.6)
      if (pd[["lyo"]] == "4x4") {
        par(mfrow = c(4, 4), mar = c(0, 0.2, 0.5, 0.2) * 0.6)
        ntop <- 16
      }
      if (pd[["lyo"]] == "6x6") {
        par(mfrow = c(6, 6), mar = c(0, 0.2, 0.5, 0.2) * 0.6)
        ntop <- 36
      }

      i <- 1
      sel <- NULL
      sel <- head(order(-colMeans(pd[["score"]]**2)), ntop)
      if (pd[["sortby"]] == "name") {
        sel <- sel[order(colnames(pd[["score"]])[sel])]
      }
      plt <- list()
      
      for (i in 1:length(sel)) {
        j <- sel[i]  
        gx <- pmax(pd[["score"]][, j], 0)
        gx <- 1 + round(15 * gx / (1e-8 + max(pd[["score"]])))
        klr0 <- klrpal[gx]
        ii <- order(gx)
        pos <- pd[["pos"]][ii,]
        tt <- colnames(pd[["score"]])[j]
        
        ii <- sample(nrow(pos))
        base::plot(
           pd[["pos"]][ii, ],
           pch = 19,
           cex = 1 * cex1,
           col = klr0[ii],
           xlim = 1.2 * range(pd[["pos"]][, 1]),
           ylim = 1.2 * range(pd[["pos"]][, 2]),
           fg = gray(0.8),
           bty = "o",
           xaxt = "n",
           yaxt = "n",
           xlab = "",
           ylab = ""
        )
        legend("topleft",
          legend = colnames(pd[["score"]])[j], bg = "#AAAAAA88",
          cex = 1.2, text.font = 1, y.intersp = 0.8, bty = "n",
          inset = c(-0.05, -0.0)
          )
      }
      refset <- input$refset
      mtext(refset, outer = TRUE, line = 0.5, cex = 1.0)     
    }

    get_ggplots <- function() {
      pd <- plot_data()
      shiny::req(pd)
      
      cex1 <- 1.2
      cex.bin <- cut(nrow(pd[["pos"]]), breaks = c(-1, 40, 200, 1000, 1e10))
      cex1 <- 0.9 * c(2.2, 1.1, 0.6, 0.3)[cex.bin]
      klrpal <- colorRampPalette(c("grey90", "grey50", "red3"))(16)
      klrpal <- paste0(gplots::col2hex(klrpal), "66")

      ntop <- 25
      if (pd[["lyo"]] == "4x4") ntop <- 16
      if (pd[["lyo"]] == "6x6") ntop <- 36

      i <- 1
      sel <- NULL
      sel <- head(order(-colMeans(pd[["score"]]**2)), ntop)
      if (pd[["sortby"]] == "name") {
        sel <- sel[order(colnames(pd[["score"]])[sel])]
      }
      
      plt <- list()      
      for (i in 1:length(sel)) {
        j <- sel[i]  
        gx <- pmax(pd[["score"]][, j], 0)
        gx <- 1 + round(15 * gx / (1e-8 + max(pd[["score"]])))
        klr0 <- klrpal[gx]
        ii <- order(gx)
        pos <- pd[["pos"]][ii,]
        tt <- colnames(pd[["score"]])[j]
        ## ------- start plot ----------       
        p <- pgx.scatterPlotXY.GGPLOT(
          pos,
          var = gx,
          col = klrpal,
          cex = 0.8*cex1,
          xlab = "",
          ylab = "",
          xlim = 1.2*range(pd[["pos"]][, 1]),
          ylim = 1.2*range(pd[["pos"]][, 2]),
          ##axis = FALSE,
          title = tt,
          cex.title = cex*0.85,
          ##title.y = 0.9,
          ##cex.clust = cex*0.8,
          label.clusters = FALSE,
          legend = FALSE,
          barscale = 0.3
        ) 

        plt[[i]] <- p
      }
      return(plt)
    }

    get_plotly <- function() {
      pd <- plot_data()
      shiny::req(pd)
      
      cex1 <- 1.2
      cex.bin <- cut(nrow(pd[["pos"]]), breaks = c(-1, 40, 200, 1000, 1e10))
      cex1 <- 0.6 * c(2.2, 1.1, 0.6, 0.3)[cex.bin]
      klrpal <- colorRampPalette(c("grey90", "grey50", "red3"))(16)
      klrpal <- paste0(gplots::col2hex(klrpal), "66")

      ntop <- 25
      if (pd[["lyo"]] == "4x4") ntop <- 16
      if (pd[["lyo"]] == "6x6") ntop <- 36

      i <- 1
      sel <- NULL
      sel <- head(order(-colMeans(pd[["score"]]**2)), ntop)
      if (pd[["sortby"]] == "name") {
        sel <- sel[order(colnames(pd[["score"]])[sel])]
      }
      plt <- list()
      
      for (i in 1:length(sel)) {
        j <- sel[i]  
        gx <- pmax(pd[["score"]][, j], 0)
        gx <- 1 + round(15 * gx / (1e-8 + max(pd[["score"]])))
        klr0 <- klrpal[gx]
        ii <- order(gx)
        pos <- pd[["pos"]][ii,]
        tt <- colnames(pd[["score"]])[j]
        ## ------- start plot ----------       
        p <- pgx.scatterPlotXY.PLOTLY(
          pos,
          var = gx,
          col = klrpal,
          cex = 0.6*cex1,
          xlab = "",
          ylab = "",
          xlim = 1.2*range(pd[["pos"]][, 1]),
          ylim = 1.2*range(pd[["pos"]][, 2]),
          axis = FALSE,
          title = tt,
          cex.title = cex1*0.5,
          title.y = 0.9,
#         cex.clust = cex1*0.8,
          label.clusters = FALSE,
          legend = FALSE,
          gridcolor = "fff",
          bgcolor = "#f8f8f8"          
        ) %>% plotly::layout(
          ## showlegend = TRUE,
          plot_bgcolor = "#f8f8f8"
        )
        
        plt[[i]] <- p
      }
      return(plt)
    }
    
    plotly.RENDER <- function() {
      pd <- plot_data()  
      plt <- get_plotly()       
      nr <- 5
      if (pd[["lyo"]] == "4x4") nr <- 4
      if (pd[["lyo"]] == "6x6") nr <- 6      
      fig <- plotly::subplot(
        plt,
        nrows = nr,
        margin = 0.01
      ) %>% plotly::layout(
        title = list(text=pd$refset, size=16)
##        margin = c(l=0,r=0,b=0,t=30) # lrbt
      ) ## %>% plotly_default()
      return(fig)
    }

    plotly_modal.RENDER <- function() {
      fig <- plotly.RENDER() %>%
        plotly::layout(
          margin = list(l=0,r=0,b=0,t=50) # lfbt  
        ) %>%
          plotly_modal_default()      
      return(fig)
    }

    ggplot.RENDER <- function() {
      pd <- plot_data()  
      plt <- get_ggplots()       
      nr <- 5
      if (pd[["lyo"]] == "4x4") nr <- 4
      if (pd[["lyo"]] == "6x6") nr <- 6      
      fig <- gridExtra::grid.arrange(
        grobs = plt,
        nrow = nr,
        ncol = nr
      )
      return(fig)
    }
    
    PlotModuleServer(
      id = "plot",        
      func = ggplot.RENDER,
      ##func = plotly.RENDER,
      ##func2 = plotly_modal.RENDER,
      ##plotlib = "plotly",
      plotlib = "ggplot",      
      res = c(85, 95),
      pdf.width = 12, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
