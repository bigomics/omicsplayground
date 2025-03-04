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
singlecell_plot_mappingplot_ui <- function(
    id,
    title,
    info.text,
    caption,
    label = "",
    height,
    width,
    parent) {
  ns <- shiny::NS(id)

  VIEWTYPES2 <- c("dotmap" = "dotmap", "heatmap (by method)" = "heatmap")

  mapping.opts <- shiny::tagList(
    withTooltip(shiny::selectInput(parent("view2"), "plot type:", VIEWTYPES2),
      "Specify the plot type: dotmap, or heatmap.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(shiny::selectInput(parent("refset2"), "reference:", choices = NULL),
      "Select a reference dataset for the cell type prediction.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(shiny::selectInput(parent("dcmethod2"), "method:", choices = NULL),
      "Choose a method for the cell type prediction.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(shiny::selectInput(parent("group2"), "group by:", "group", selected = NULL),
      "Group the samples/cells by grouping factor.",
      placement = "top", options = list(container = "body")
    )
  )

  PlotModuleUI(
    id = ns("plot"),
    label = label,
    info.text = info.text,
    options = mapping.opts,
    title = title,
    caption = caption,
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
singlecell_plot_mappingplot_server <- function(id,
                                               pgx,
                                               pfGetClusterPositions,
                                               getDeconvResults2,
                                               grpvar, # input$group2
                                               refset, # input$refset2
                                               group, # input$group2
                                               view, # input$view2
                                               watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      shiny::req(pfGetClusterPositions())
      if (!("deconv" %in% names(pgx)) || length(pgx$deconv) == 0) {
        shiny::validate(shiny::need(FALSE, "Cell type profile requires deconvolution"))
      }
      clust.pos <- pfGetClusterPositions()
      pos <- pgx$tsne2d
      pos <- clust.pos

      score <- pgx$deconv[["LM22"]][["meta"]]
      score <- getDeconvResults2()
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

      grpvar <- grpvar() # input$group2
      refset <- refset() # input$refset2
      view <- view()

      return(list(
        grpvar = grpvar,
        score = score,
        refset = refset,
        pgx = pgx,
        pos = pos,
        view = view
      ))
    })

    plot.render <- function() {
      pd <- plot_data()

      if (pd[["grpvar"]] != "<ungrouped>" && pd[["grpvar"]] %in% colnames(pd[["pgx"]]$samples)) {
        grp <- pd[["pgx"]]$samples[rownames(pd[["score"]]), pd[["grpvar"]]]
        shiny::validate(shiny::need(length(unique(grp)) > 1, "Filter is too restrictive, two samples at least are required."))
        pd[["pos"]] <- apply(pd[["pos"]], 2, function(x) tapply(x, grp, median))
        pd[["score"]] <- apply(pd[["score"]], 2, function(x) tapply(x, grp, mean))
        ii <- hclust(dist(pd[["score"]]))$order
        jj <- hclust(dist(t(pd[["score"]])))$order
        pd[["score"]] <- pd[["score"]][ii, jj]
      }
      b0 <- 0.1 + 0.70 * pmax(30 - ncol(pd[["score"]]), 0)

      if (pd[["view"]] == "dotmap") {
        #
        par(mfrow = c(1, 1), mar = c(0, 0, 8, 1), oma = c(1, 1, 1, 1) * 0.25)
        score3 <- pd[["score"]]**1.5
        rownames(score3) <- paste("", rownames(score3), "  ")
        tl.srt <- 90
        tl.cex <- ifelse(nrow(pd[["score"]]) > 60, 0.7, 0.85)
        if (max(sapply(rownames(score3), nchar)) > 30) tl.srt <- 45
        corrplot::corrplot(t(score3),
          mar = c(b0, 1, 4, 0.5),
          cl.lim = c(0, max(score3)), cl.pos = "n",
          tl.cex = tl.cex, tl.col = "grey20",
          tl.srt = tl.srt
        )
      }

      if (pd[["view"]] == "heatmap") {
        usermode <- "PRO"
        if (!is.null(usermode) && usermode >= "PRO") {
          kk <- head(colnames(pd[["score"]])[order(-colMeans(pd[["score"]]**2))], 18)
          kk <- intersect(colnames(pd[["score"]]), kk)
          all.scores <- pd[["pgx"]]$deconv[["LM22"]]
          all.scores <- pd[["pgx"]]$deconv[[pd[["refset"]]]]
          if (pd[["grpvar"]] != "<ungrouped>" && pd[["grpvar"]] %in% colnames(pd[["pgx"]]$samples)) {
            grp <- pd[["pgx"]]$samples[rownames(all.scores[[1]]), pd[["grpvar"]]]
            for (i in 1:length(all.scores)) {
              all.scores[[i]] <- apply(
                all.scores[[i]], 2,
                function(x) tapply(x, grp, mean)
              )
              ii <- rownames(pd[["score"]])
              all.scores[[i]] <- all.scores[[i]][ii, kk]
            }
          }

          nm <- length(all.scores)
          m <- 3
          n <- 2
          if (nm > 6) {
            m <- 3
            n <- 3
          }
          if (nm > 9) {
            m <- 4
            n <- 3
          }
          rr <- 2 + max(nchar(colnames(pd[["score"]]))) / 2
          par(mfrow = c(m, n), mar = c(0, 0.3, 2, 0.3), oma = c(10, 0, 0, rr), xpd = TRUE)
          k <- 1
          for (k in 1:length(all.scores)) {
            ii <- rownames(pd[["score"]])
            score1 <- all.scores[[k]][ii, kk]
            if (k %% n != 0) colnames(score1) <- rep("", ncol(score1))
            if ((k - 1) %/% n != (nm - 1) %/% n) rownames(score1) <- rep("", nrow(score1))
            score1 <- score1 / (1e-8 + rowSums(score1))
            if (nrow(score1) > 100) rownames(score1) <- rep("", nrow(score1))
            playbase::gx.imagemap(t(score1**1), cex = 0.85, main = "", clust = FALSE)
            title(main = names(all.scores)[k], cex.main = 1.1, line = 0.4, font.main = 1)
          }
        } else {
          score1 <- pd[["score"]]
          score1 <- score1 / (1e-8 + rowSums(score1))
          if (nrow(score1) > 100) rownames(score1) <- rep("", nrow(pd[["score"]]))
          playbase::gx.heatmap(t(score1**2),
            scale = "none",
            cexRow = 1, cexCol = 0.6, col = heat.colors(16),
            mar = c(b0, 15), key = FALSE, keysize = 0.5
          )
        }
      }
    }

    PlotModuleServer(
      "plot",
      func = plot.render,
      res = c(85, 95),
      pdf.width = 8, pdf.height = 8,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
