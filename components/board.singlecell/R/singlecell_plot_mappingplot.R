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
  parent
) {
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
    outputFunc = shiny::plotOutput,
    outputFunc2 = shiny::plotOutput,
    download.fmt = c("png", "pdf", "svg", "csv"),
    height = height,
    width = width,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "gradient"
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

      # Apply grouping aggregation if needed
      if (grpvar != "<ungrouped>" && grpvar %in% colnames(pgx$samples)) {
        grp <- pgx$samples[rownames(score), grpvar]
        shiny::validate(shiny::need(length(unique(grp)) > 1, "Filter is too restrictive, two samples at least are required."))
        pos <- apply(pos, 2, function(x) tapply(x, grp, median))
        score <- apply(score, 2, function(x) tapply(x, grp, mean))
        ii <- hclust(dist(score))$order
        jj <- hclust(dist(t(score)))$order
        score <- score[ii, jj]
      }

      # Prepare final plot data based on view type
      final_data <- list()
      if (view == "dotmap") {
        # For dotmap: transform the score data
        score_transformed <- score**1.5
        rownames(score_transformed) <- paste("", rownames(score_transformed), "  ")
        final_data$plot_matrix <- t(score_transformed)
        final_data$plot_type <- "dotmap"
      }

      if (view == "heatmap") {
        usermode <- "PRO"

        if (!is.null(usermode) && usermode >= "PRO") {
          # PRO mode: prepare all.scores with all methods
          kk <- head(colnames(score)[order(-colMeans(score**2))], 18)
          kk <- intersect(colnames(score), kk)
          all.scores <- pgx$deconv[[refset]]

          if (grpvar != "<ungrouped>" && grpvar %in% colnames(pgx$samples)) {
            grp <- pgx$samples[rownames(all.scores[[1]]), grpvar]
            for (i in 1:length(all.scores)) {
              all.scores[[i]] <- apply(
                all.scores[[i]], 2,
                function(x) tapply(x, grp, mean)
              )
              ii <- rownames(score)
              all.scores[[i]] <- all.scores[[i]][ii, kk]
            }
          }

          # Process each method's scores for plotting
          processed_scores <- list()
          for (k in 1:length(all.scores)) {
            ii <- rownames(score)
            score1 <- all.scores[[k]][ii, kk]
            score1 <- score1 / (1e-8 + rowSums(score1))
            processed_scores[[names(all.scores)[k]]] <- t(score1**1)
          }

          final_data$plot_matrices <- processed_scores
          final_data$plot_type <- "heatmap_pro"
          final_data$selected_features <- kk
        } else {
          # Non-PRO mode: single heatmap
          score1 <- score
          score1 <- score1 / (1e-8 + rowSums(score1))
          final_data$plot_matrix <- t(score1**2)
          final_data$plot_type <- "heatmap_basic"
        }
      }

      return(list(
        grpvar = grpvar,
        score = score,
        refset = refset,
        pgx = pgx,
        pos = pos,
        view = view,
        final_data = final_data
      ))
    })

    plot.render <- function() {
      pd <- plot_data()
      if (is.null(pd) || is.null(pd$final_data)) {
        return(NULL)
      }

      ## Editor: color palette (white → user color)
      color_low <- if (!is.null(input$color_low)) input$color_low else get_color_theme()$secondary
      col_pal <- grDevices::colorRampPalette(c("white", "grey90", color_low))

      final_data <- pd$final_data
      b0 <- 0.1 + 0.70 * pmax(30 - ncol(pd[["score"]]), 0)

      if (final_data$plot_type == "dotmap") {
        # Dotmap visualization
        par(mfrow = c(1, 1), mar = c(0, 0, 8, 1), oma = c(1, 1, 1, 1) * 0.25)
        score3 <- final_data$plot_matrix
        tl.srt <- 90
        tl.cex <- ifelse(nrow(score3) > 60, 0.7, 0.85)
        if (max(sapply(rownames(score3), nchar)) > 30) tl.srt <- 45
        corrplot::corrplot(score3,
          col = col_pal(100),
          mar = c(b0, 1, 4, 0.5),
          cl.lim = c(0, max(score3)), cl.pos = "n",
          tl.cex = tl.cex, tl.col = "grey20",
          tl.srt = tl.srt
        )
      }

      if (final_data$plot_type == "heatmap_pro") {
        # PRO mode heatmap with multiple methods
        processed_scores <- final_data$plot_matrices
        nm <- length(processed_scores)
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
        rr <- 2 + max(nchar(final_data$selected_features)) / 2
        par(mfrow = c(m, n), mar = c(0, 0.3, 2, 0.3), oma = c(10, 0, 0, rr), xpd = TRUE)

        k <- 1
        for (k in 1:length(processed_scores)) {
          method_name <- names(processed_scores)[k]
          score_matrix <- processed_scores[[k]]
          # Handle subplot labeling
          if ((k - 1) %/% n != (nm - 1) %/% n) colnames(score_matrix) <- rep("", ncol(score_matrix))
          if (k %% n != 0) rownames(score_matrix) <- rep("", nrow(score_matrix))
          if (nrow(score_matrix) > 100) rownames(score_matrix) <- rep("", nrow(score_matrix))

          playbase::gx.imagemap(score_matrix, cex = 0.85, main = "", clust = FALSE,
            col = col_pal(64))
          title(main = method_name, cex.main = 1.1, line = 0.4, font.main = 1)
          k <- k + 1
        }
      }

      if (final_data$plot_type == "heatmap_basic") {
        # Basic heatmap
        score_matrix <- final_data$plot_matrix
        if (nrow(pd[["score"]]) > 100) rownames(score_matrix) <- rep("", nrow(pd[["score"]]))

        playbase::gx.heatmap(score_matrix,
          scale = "none",
          cexRow = 1, cexCol = 0.6, col = col_pal(16),
          mar = c(b0, 15), key = FALSE, keysize = 0.5
        )
      }
    }

    plot_data_csv <- function() {
      plot_result <- plot_data()
      final_plot_data <- plot_result$final_data

      if (final_plot_data$plot_type == "dotmap") {
        final_plot_data$plot_matrix
      } else if (final_plot_data$plot_type == "heatmap_pro") {
        do.call(rbind, final_plot_data$plot_matrices)
      } else if (final_plot_data$plot_type == "heatmap_basic") {
        final_plot_data$plot_matrix
      }
    }

    PlotModuleServer(
      "plot",
      func = plot.render,
      res = c(85, 95),
      pdf.width = 8, pdf.height = 8,
      add.watermark = watermark,
      csvFunc = plot_data_csv,
      parent_session = session
    )
  }) ## end of moduleServer
}
