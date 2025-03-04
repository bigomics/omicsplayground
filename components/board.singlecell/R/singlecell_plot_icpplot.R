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
singlecell_plot_icpplot_ui <- function(
    id,
    title,
    info.text,
    caption,
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
        inline = TRUE
      ),
      "Choose layout.",
      placement = "top", options = list(container = "body")
    )
  )

  PlotModuleUI(
    id = ns("plot"),
    plotlib = "plotly",
    ## plotlib = "ggplot",
    label = label,
    info.text = info.text,
    title = title,
    caption = caption,
    options = icp.opts,
    download.fmt = c("png", "pdf", "svg"),
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
                                           layout, # input$layout
                                           sortby, # input$sortby
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      method <- method()
      refset <- refset()
      layout <- layout()
      sortby <- sortby()

      if (!("deconv" %in% names(pgx)) || length(pgx$deconv) == 0) {
        shiny::validate(shiny::need(FALSE, "Cell type mapping requires deconvolution"))
      }
      results <- pgx$deconv[[refset]][[method]]
      ## threshold everything (because DCQ can be negative!!!)
      results <- pmax(results, 0)

      shiny::req(pfGetClusterPositions())
      clust.pos <- pfGetClusterPositions()
      pos <- clust.pos
      score <- results
      if (is.null(score) || length(score) == 0) {
        return(NULL)
      }

      ## normalize
      score <- score[rownames(pos), ]
      score[is.na(score)] <- 0
      score <- pmax(score, 0)
      score <- score / (1e-20 + rowSums(score, na.rm = TRUE))
      score <- tanh(score / mean(abs(score)))
      score <- score / max(score, na.rm = TRUE)

      ## take top10 features
      jj.top <- unique(as.vector(apply(score, 1, function(x) head(order(-x), 10))))
      score <- score[, jj.top]
      score <- score[, order(-colMeans(score**2))]
      score <- score[, 1:min(50, ncol(score))]
      ii <- hclust(dist(score))$order
      jj <- hclust(dist(t(score)))$order
      score <- score[ii, jj]
      pos <- pos[rownames(score), ]

      # Return list
      pd <- list(
        score = score,
        pos = pos,
        layout = layout,
        refset = refset,
        sortby = sortby
      )
      return(pd)
    })

    get_ggplots <- function(cex = 1) {
      pd <- plot_data()
      shiny::req(pd)

      cex1 <- 1.2
      cex.bin <- cut(nrow(pd[["pos"]]), breaks = c(-1, 40, 200, 1000, 1e10))
      cex1 <- cex * c(2.2, 1.1, 0.6, 0.3)[cex.bin]
      klrpal <- colorRampPalette(c("grey95", "grey65", "red3"))(16)
      klrpal <- paste0(gplots::col2hex(klrpal), "66") ## add opacity...

      ntop <- 25
      if (pd[["layout"]] == "4x4") ntop <- 16
      if (pd[["layout"]] == "6x6") ntop <- 36

      i <- 1
      sel <- NULL
      sel <- head(order(-colMeans(pd[["score"]]**2)), ntop)
      if (pd[["sortby"]] == "name") {
        sel <- sel[order(colnames(pd[["score"]])[sel])]
      }

      cmin <- min(pd[["score"]], na.rm = TRUE)
      cmax <- max(pd[["score"]], na.rm = TRUE)

      plt <- list()
      for (i in 1:length(sel)) {
        j <- sel[i]
        gx <- pmax(pd[["score"]][, j], 0)
        pos <- pd[["pos"]]
        tt <- colnames(pd[["score"]])[j]

        if (i == 1) {
          legend <- TRUE
        } else {
          legend <- FALSE
        }
        ## ------- start plot ----------
        p <- playbase::pgx.scatterPlotXY.GGPLOT(
          pos,
          var = gx,
          col = klrpal,
          cex = 0.6 * cex1,
          xlab = "",
          ylab = "",
          zlim = c(cmin, cmax),
          #          cmin = cmin,
          #          cmax = cmax,
          xlim = 1.2 * range(pd[["pos"]][, 1]),
          ylim = 1.2 * range(pd[["pos"]][, 2]),
          axis = FALSE,
          title = tt,
          cex.title = 0.55,
          label.clusters = FALSE,
          legend = legend,
          gridcolor = "#ffffff",
          bgcolor = "#f8f8f8",
          box = TRUE,
          guide = "legend"
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
      klrpal <- colorRampPalette(c("grey95", "grey65", "red3"))(16)
      klrpal <- paste0(gplots::col2hex(klrpal), "66")

      cmin <- min(pd[["score"]], na.rm = TRUE)
      cmax <- max(pd[["score"]], na.rm = TRUE)

      ntop <- 25
      if (pd[["layout"]] == "4x4") ntop <- 16
      if (pd[["layout"]] == "6x6") ntop <- 36

      i <- 1
      sel <- NULL
      sel <- head(order(-colMeans(pd[["score"]]**2)), ntop)
      if (pd[["sortby"]] == "name") {
        sel <- sel[order(colnames(pd[["score"]])[sel])]
      }

      plt <- vector("list", length(sel))
      for (i in 1:length(sel)) {
        j <- sel[i]
        gx <- pmax(pd[["score"]][, j], 0)
        gx <- 1 + round(15 * gx / (1e-8 + max(pd[["score"]])))
        if (length(unique(gx)) == 1) {
          col <- klrpal[gx[1]]
          col <- substr(col, 1, 7)
        } else {
          col <- klrpal
        }
        ii <- order(gx)
        pos <- pd[["pos"]][ii, ]
        tt <- colnames(pd[["score"]])[j]
        ## ------- start plot ----------
        p <- playbase::pgx.scatterPlotXY.PLOTLY(
          pos,
          var = gx,
          col = col,
          zlim = c(cmin, cmax),
          cex = 1 * cex1,
          xlab = "",
          ylab = "",
          xlim = 1.25 * range(pd[["pos"]][, 1]),
          ylim = 1.25 * range(pd[["pos"]][, 2]),
          axis = FALSE,
          title = tt,
          cex.title = 0.9,
          title.y = 0.85,
          label.clusters = FALSE,
          legend = FALSE,
          box = TRUE,
          gridcolor = "#ffffff",
          bgcolor = "#f8f8f8",
        )
        plt[[i]] <- p
      }
      return(plt)
    }


    plotly.RENDER <- function() {
      pd <- plot_data()
      plt <- get_plotly()
      nr <- 2
      if (pd[["layout"]] == "4x4") nr <- 4
      if (pd[["layout"]] == "6x6") nr <- 6
      plt <- head(plt, nr * nr)

      fig <- plotly::subplot(
        plt,
        nrows = nr,
        margin = 0.01
      ) %>%
        plotly::layout(
          title = list(text = pd$refset, size = 14),
          margin = list(l = 0, r = 0, b = 0, t = 30) # lrbt
        )
      return(fig)
    }

    plotly_modal.RENDER <- function() {
      fig <- plotly.RENDER() %>%
        plotly::layout(
          margin = list(l = 0, r = 0, b = 0, t = 50) # lfbt
        ) %>%
        plotly_modal_default()
      return(fig)
    }


    ggplot.RENDER <- function() {
      pd <- plot_data()
      plt <- get_ggplots()
      shiny::req(plt)
      nr <- ceiling(sqrt(length(plt)))
      title <- pd$refset
      fig <- gridExtra::grid.arrange(
        grobs = plt,
        nrow = nr,
        ncol = nr,
        padding = unit(0.01, "line"),
        top = textGrob(title, gp = gpar(fontsize = 15))
      )
      return(fig)
    }

    PlotModuleServer(
      id = "plot",
      func = plotly.RENDER,
      func2 = plotly_modal.RENDER,
      plotlib = "plotly",
      ##      func = ggplot.RENDER,
      ##      plotlib = "ggplot",
      res = c(85, 95),
      pdf.width = 12, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
