##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' Clustering plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
clustering_plot_splitmap_ui <- function(
  id,
  label = "",
  title,
  caption,
  info.text,
  info.methods,
  info.references,
  info.extra_link,
  height,
  width
) {
  ns <- shiny::NS(id)

  splitmap_opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(
        ns("plot_type"), "Plot type:",
        choices = c("heatmap", "sample correlation"),
        inline = TRUE
      ),
      "Show gene-sample heatmap or sample-sample correlation plot."
    ),
    shiny::hr(),
    withTooltip(
      shiny::radioButtons(
        ns("hm_scale"), "Row scaling:",
        choices = c("center" = "row.center", "absolute" = "none", "z-score" = "row"),
        inline = TRUE
      ),
      "Show relative (i.e. mean-centered), absolute expression values or batch-mean-centered.",
      placement = "right", options = list(container = "body")
    ),
    shiny::hr(),
    withTooltip(
      shiny::checkboxInput(
        ns("show_legend"), "show legend",
        value = TRUE
      ), "Show or hide the legend."
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("show_rownames"), "show row names",
        value = TRUE
      ), "Show or hide the row names (features)."
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("show_colnames"), "show column names",
        value = TRUE
      ), "Show or hide the column names (samples)."
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    plotlib = c("plotly", "base"),
    info.text = info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption = caption,
    options = splitmap_opts,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height,
    cards = TRUE,
    card_names = c("dynamic", "static"),
    editor = TRUE,
    ns_parent = ns,
    plot_type = "heatmap"
  )
}

#' Clustering plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#' @param watermark
#'
#'
#'
#' @export
clustering_plot_splitmap_server <- function(id,
                                            pgx,
                                            getTopMatrix,
                                            selected_phenotypes,
                                            hm_level,
                                            hm_ntop,
                                            hm_topmode,
                                            hm_clustk,
                                            hm_average_group,
                                            watermark = FALSE,
                                            labeltype) {
  moduleServer(id, function(input, output, session) {
    fullH <- 850

    ns <- session$ns

    shiny::observeEvent(pgx$Y, {
      if (nrow(pgx$Y) > 100) { # Put cexCol (heatmap) to 0 if more than 100 samples
        shiny::updateNumericInput(session, "hm_cexCol", value = 0)
      } else {
        shiny::updateNumericInput(session, "hm_cexCol", value = 1)
      }
    })

    ## Editor: get custom group order for split heatmaps
    get_splitx_order <- function() {
      custom_order <- input$hm_group_order
      if (is.null(custom_order) || length(custom_order) == 0) {
        return(NULL)
      }
      custom_order
    }

    ## Editor: drag-and-drop group ordering for split heatmaps
    output$hm_group_order_ui <- shiny::renderUI({
      pd <- plot_data()
      shiny::req(pd)
      splitx <- pd$filt$grp
      if (is.null(splitx) || length(unique(splitx)) <= 1) {
        return(shiny::tags$em("Split the heatmap by a variable to reorder groups."))
      }
      grp_levels <- unique(splitx)
      rank_list_ui(grp_levels, ns, input_id = "hm_group_order")
    })

    ## Editor: column clustering for the (non-split) heatmap. The same
    ## hclust object seeds the sample-order list AND is handed to the static
    ## plot via cluster_columns, so the list and the plot show the exact same
    ## order while keeping the column dendrogram. NULL when the heatmap is
    ## split (columns are then clustered within each group) or too small.
    hm_col_hclust <- shiny::reactive({
      filt <- getTopMatrix()
      shiny::req(filt)
      splitx <- filt$grp
      if (!is.null(splitx) && length(unique(splitx)) > 1) {
        return(NULL)
      }
      m <- filt$mat
      if (is.null(m) || ncol(m) < 2) {
        return(NULL)
      }
      ## scale rows the same way the plot does, so the dendrogram matches
      scale.mode <- input$hm_scale
      if (!is.null(scale.mode)) {
        if (scale.mode == "row.center") m <- m - rowMeans(m, na.rm = TRUE)
        if (scale.mode == "row") m <- t(scale(t(m)))
      }
      d <- stats::dist(t(m))
      if (anyNA(d)) d[is.na(d)] <- max(d, na.rm = TRUE)
      fastcluster::hclust(d, method = "ward.D2")
    })

    ## Editor: default (clustered) display order of the samples. Used both to
    ## seed the drag-and-drop list and as the baseline for detecting whether
    ## the user has actually rearranged it.
    default_sample_order <- function() {
      filt <- getTopMatrix()
      if (is.null(filt) || is.null(filt$mat)) {
        return(NULL)
      }
      cn <- colnames(filt$mat)
      hc <- hm_col_hclust()
      if (is.null(hc)) {
        return(cn)
      }
      cn[hc$order]
    }

    ## Editor: get custom individual-sample order. Returns the reordered
    ## sample vector only when the user has actually rearranged the list
    ## (i.e. it differs from the default clustered order shown in the plot);
    ## otherwise NULL so the default column clustering is preserved. Disabled
    ## while the heatmap is split into groups (use Group Order instead).
    get_sample_order <- function() {
      custom_order <- input$hm_sample_order
      if (is.null(custom_order) || length(custom_order) == 0) {
        return(NULL)
      }
      filt <- getTopMatrix()
      if (is.null(filt)) {
        return(NULL)
      }
      splitx <- filt$grp
      if (!is.null(splitx) && length(unique(splitx)) > 1) {
        return(NULL)
      }
      custom_order <- as.character(custom_order)
      if (identical(custom_order, default_sample_order())) {
        return(NULL)
      }
      custom_order
    }

    ## Editor: reorder heatmap columns by the custom sample order. Permutes
    ## the data matrix and the sample-split vector together so they stay
    ## aligned. 'custom' flags whether a reorder was applied, so the caller
    ## can switch off column clustering to honour the manual order.
    apply_sample_order <- function(zx, splitx) {
      ord <- get_sample_order()
      if (is.null(ord)) {
        return(list(zx = zx, splitx = splitx, custom = FALSE))
      }
      cn <- colnames(zx)
      neworder <- c(intersect(ord, cn), setdiff(cn, ord))
      perm <- match(neworder, cn)
      zx <- zx[, perm, drop = FALSE]
      if (!is.null(splitx)) splitx <- splitx[perm]
      list(zx = zx, splitx = splitx, custom = TRUE)
    }

    ## Editor: drag-and-drop individual-sample ordering. Only offered when
    ## the heatmap is not split; when split, columns are ordered via groups.
    output$hm_sample_order_ui <- shiny::renderUI({
      filt <- getTopMatrix()
      shiny::req(filt)
      splitx <- filt$grp
      if (!is.null(splitx) && length(unique(splitx)) > 1) {
        return(shiny::tags$em("Sample ordering is unavailable while the heatmap is split. Remove the split to reorder individual samples, or use Group Order."))
      }
      ## seed the list in the same (clustered) order the plot displays
      samples <- default_sample_order()
      if (is.null(samples) || length(samples) <= 1) {
        return(shiny::tags$em("Not enough samples to reorder."))
      }
      rank_list_ui(samples, ns, input_id = "hm_sample_order")
    })

    plot_data <- shiny::reactive({
      filt <- getTopMatrix()
      shiny::req(filt)

      zx <- filt$mat
      annot <- filt$annot
      zx.idx <- filt$idx

      ## For large data (eg scRNAseq)
      if (pgx$datatype == "scRNAseq" && ncol(zx) > 3000) {
        dbg("[clustering_plot_splitmap] More than 3K samples detected in data matrix.")
        dbg("[clustering_plot_splitmap] Balanced downsampling of 40 cells per celltype.")
        grp <- annot$celltype
        sel <- tapply(1:ncol(zx), grp, function(ii) head(sample(ii), 40))
        ii <- unname(unlist(sel))
        zx <- zx[, ii]
        annot <- annot[colnames(zx), , drop = FALSE]
        filt$mat <- zx
        filt$annot <- annot
        filt$samples <- colnames(zx)
        filt$grp <- annot$celltype
      }

      shiny::validate(shiny::need(
        ncol(zx) > 1, "Filtering too restrictive. Please change 'Filter samples' settings."
      ))

      sample_cor <- FALSE
      if (input$plot_type == "sample correlation") {
        sample_cor <- TRUE
        X <- zx
        if (input$hm_scale == "row.center") X <- X - rowMeans(X, na.rm = TRUE)
        if (input$hm_scale == "row") X <- t(scale(t(X)))
        zx <- cor(X, method = "pearson", use = "pairwise.complete.obs")
        D <- as.dist(1 - zx)
        D[which(is.nan(D) | is.na(D))] <- 1
        hc <- fastcluster::hclust(D, method = "ward.D2")
        zx.idx <- try(paste0("S", cutree(hc, hm_clustk())), silent = TRUE)
        if (inherits(zx.idx, "try-error")) zx.idx <- paste0("S", cutree(hc, 1))
      }

      return(list(zx = zx, annot = annot, zx.idx = zx.idx, filt = filt, sample_cor = sample_cor))
    })

    base_splitmap.RENDER <- function() {
      pd <- plot_data()
      zx <- pd[["zx"]]
      annot <- pd[["annot"]]
      zx.idx <- pd[["zx.idx"]]
      filt <- pd[["filt"]]
      sample_cor <- pd[["sample_cor"]]

      if (nrow(zx) <= 1) {
        return(NULL)
      }
      show_rownames <- input$show_rownames

      # Use editor settings for text sizes. normalize to base size 10
      label_cex <- ifelse(!is.null(input$label_size), input$label_size / 10, 1)

      ## split genes dimension in 5 groups
      splity <- 5
      splity <- 6
      if (!is.null(zx.idx)) {
        splity <- zx.idx
      }

      ## split samples
      splitx <- filt$grp

      ## Editor: apply custom individual-sample order (heatmap mode only;
      ## skipped for sample-correlation where columns and rows are samples)
      sample_order_custom <- FALSE
      if (!sample_cor) {
        sample_ord <- apply_sample_order(zx, splitx)
        zx <- sample_ord$zx
        splitx <- sample_ord$splitx
        sample_order_custom <- sample_ord$custom
      }

      # Calculate text sizes
      cex0 <- ifelse(!is.null(splitx) && length(splitx) <= 10, 1.05, 0.85) # fixed title size scaling
      cex1 <- label_cex * ifelse(ncol(zx) > 200, 0, ifelse(ncol(zx) > 100, 0.5, ifelse(ncol(zx) > 50, 0.75, 1)))
      cex2 <- label_cex * ifelse(nrow(zx) > 60, 0.8, 1)

      shiny::validate(shiny::need(
        !is.null(cex1) && !is.null(cex2),
        "Text size settings must not be empty."
      ))

      # Set column text size to 0 if show_colnames is FALSE
      if (!input$show_colnames) cex1 <- 0
      show_colnames <- (cex1 > 0)

      annotF <- data.frame(as.list(annot), stringsAsFactors = TRUE, check.names = FALSE)
      rownames(annotF) <- rownames(annot)

      # Select annot to display (user input)
      sel <- selected_phenotypes()
      sel <- intersect(sel, colnames(annotF))
      if (length(sel) == 0) {
        annot <- NULL
      } else {
        annot <- annotF[, sel, drop = FALSE]
      }

      if (hm_level() == "gene") {
        rownames(zx) <- sub(".*:", "", rownames(zx))
        rownames(zx) <- playbase::probe2symbol(rownames(zx), pgx$genes, labeltype(), fill_na = TRUE)
      }

      if (hm_level() == "geneset") {
        rownames(zx) <- sub("HALLMARK:HALLMARK_", "HALLMARK:", rownames(zx))
        rownames(zx) <- gsub(playdata::GSET_PREFIX_REGEX, "", rownames(zx))
        rownames(zx) <- tolower(rownames(zx))
      }
      rownames(zx) <- substring(rownames(zx), 1, 50) ## cut long names...

      crot <- 0
      totnchar <- nchar(paste0(unique(splitx), collapse = ""))
      nx <- length(unique(splitx))
      if (!is.null(splitx) & (totnchar > 44 || nx >= 6)) crot <- 90

      nrownames <- 9999
      if (cex2 == 0) nrownames <- 0

      shiny::showNotification("Rendering heatmap...")

      margin_top <- ifelse(input$margin_checkbox && !is.na(input$margin_top), input$margin_top, 5)
      margin_right <- ifelse(input$margin_checkbox && !is.na(input$margin_right), input$margin_right, 5)
      margin_bottom <- ifelse(input$margin_checkbox && !is.na(input$margin_bottom), input$margin_bottom, 5)
      margin_left <- ifelse(input$margin_checkbox && !is.na(input$margin_left), input$margin_left, 5)

      scale.mode <- input$hm_scale
      if (sample_cor) {
        splity <- NULL
        splity <- splitx
        scale.mode <- "none"
        cluster_rows <- FALSE
        cluster_columns <- FALSE
        zlim <- c(min(zx, na.rm = TRUE), max(zx, na.rm = TRUE))
      } else {
        cluster_rows <- TRUE
        if (sample_order_custom) {
          ## manual order: columns already permuted, don't recluster
          cluster_columns <- FALSE
        } else {
          ## force the same column clustering that seeds the editor list so
          ## the plot and the list agree; falls back to default clustering
          ## (TRUE) when split. A precomputed hclust keeps the dendrogram.
          col_hc <- hm_col_hclust()
          cluster_columns <- if (!is.null(col_hc)) col_hc else TRUE
        }
        zlim <- NULL
      }

      playbase::gx.splitmap(
        zx,
        split = splity,
        splitx = splitx,
        splitx_order = get_splitx_order(),
        scale = scale.mode,
        show_legend = input$show_legend,
        show_colnames = show_colnames,
        column_title_rot = crot,
        column_names_rot = input$column_names_rot,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        color_low = input$color_low,
        color_mid = input$color_mid,
        color_high = input$color_high,
        zlim = zlim,
        show_rownames = show_rownames,
        rownames_width = input$rownames_width,
        softmax = 0,
        na_col = "green",
        na_text = "x",
        title_cex = 1,
        cexCol = cex1,
        cexRow = cex2,
        annot.cex = input$annot_cex / 10,
        col.annot = annot,
        row.annot = NULL,
        annot.ht = 2.3,
        key.offset = c(0.89, 1.01),
        main = " ",
        nmax = -1,
        mar = c(margin_bottom, margin_left, margin_top, margin_right)
      )
      p <- grDevices::recordPlot()
      p
    }

    plotly_splitmap.RENDER_get <- function() {
      shiny::req(pgx$genes)

      ## -------------- variable to split samples
      plt <- NULL
      pd <- plot_data()
      filt <- pd[["filt"]]
      X <- pd[["zx"]]
      annot <- pd[["annot"]]
      splity <- pd[["zx.idx"]]
      sample_cor <- pd[["sample_cor"]]

      ## sample clustering index
      splitx <- filt$grp

      ## Editor: apply custom individual-sample order (heatmap mode only;
      ## skipped for sample-correlation where columns and rows are samples)
      sample_order_custom <- FALSE
      if (!sample_cor) {
        sample_ord <- apply_sample_order(X, splitx)
        X <- sample_ord$zx
        splitx <- sample_ord$splitx
        sample_order_custom <- sample_ord$custom
      }

      ## iheatmapr needs factors for sharing between groups
      annotF <- data.frame(as.list(annot), stringsAsFactors = TRUE, check.names = FALSE)
      rownames(annotF) <- rownames(annot)

      sel <- selected_phenotypes()
      sel <- intersect(sel, colnames(annotF))
      if (length(sel) == 0) {
        annotF <- NULL
      } else {
        annotF <- annotF[, sel, drop = FALSE]
      }

      # Use editor settings for text sizes
      label_cex <- ifelse(!is.null(input$label_size), input$label_size / 10, 1) # normalize to base size 10

      # Adjust label size based on data dimensions
      col_cex <- label_cex
      if (ncol(X) < 50) col_cex <- label_cex * 1.2
      if (ncol(X) > 50) col_cex <- label_cex * 0.75
      if (ncol(X) > 100) col_cex <- label_cex * 0.5
      if (ncol(X) > 200) col_cex <- 0

      # Set column text size to 0 if show_colnames is FALSE
      if (!input$show_colnames) col_cex <- 0

      row_cex <- label_cex
      if (nrow(X) < 60) row_cex <- label_cex * 1.2
      if (nrow(X) > 60) row_cex <- label_cex * 0.8

      tooltips <- NULL
      if (hm_level() == "gene") {
        getInfo <- function(g) {
          aa <- paste0(
            "<b>", pgx$genes[g, "feature"], "</b> ",
            "(", pgx$genes[g, "symbol"], "): ",
            pgx$genes[g, "gene_title"], "."
          )
          playbase::breakstring2(aa, 50, brk = "<br>")
        }
        tooltips <- sapply(rownames(X), getInfo)
        labeled_features <- NULL
        symbol <- playbase::probe2symbol(rownames(X), pgx$genes, labeltype(), fill_na = TRUE)
        rownames(X) <- symbol
      } else {
        aa <- gsub("_", " ", rownames(X)) ## just geneset names
        tooltips <- sapply(aa, function(x) {
          playbase::breakstring2(x, 50, brk = "<br>")
        })
      }
      shiny::showNotification("Rendering iHeatmap...")

      if (sample_cor) {
        idx <- splitx
        if (hm_average_group()) {
          idx <- NULL
          splitx <- NULL
        }
        scale.mode <- "none"
        zlim <- c(min(X, na.rm = TRUE), max(X, na.rm = TRUE))
        symm <- TRUE
      } else {
        idx <- splity
        scale.mode <- input$hm_scale
        zlim <- NULL
        symm <- FALSE
      }
      plt <- playbase::pgx.splitHeatmapFromMatrix(
        X = X,
        annot = annotF,
        ytips = tooltips,
        idx = idx,
        splitx = splitx,
        splitx_order = get_splitx_order(),
        col_clust = !sample_order_custom,
        scale = scale.mode,
        zlim = zlim,
        symm = symm,
        row_annot_width = 0.025,
        rowcex = ifelse(input$show_rownames, row_cex, 0),
        colcex = col_cex,
        show_legend = input$show_legend,
        return_x_matrix = TRUE,
        heatmap_colors = extract_heatmap_lmh_colors(input)
      )
      return(plt)
    }

    plotly_splitmap.RENDER <- function() {
      plt <- plotly_splitmap.RENDER_get()$plt
      if (any(grepl("Iheatmap", class(plt)))) {
        plt <- plt %>% iheatmapr::to_plotly_list()
        plt <- plotly::as_widget(plt)
      }
      plt <- plt %>%
        plotly::layout(margin = list(l = 10, r = 5, t = 5, b = 5))
      return(plt)
    }

    plot_grid <- list(
      list(plotlib = "plotly", func = plotly_splitmap.RENDER, card = 1),
      list(plotlib = "base", func = base_splitmap.RENDER, card = 2)
    )

    plot_data_csv <- function() {
      plotly_splitmap.RENDER_get()$X
    }

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "pltmod",
        plotlib = x$plotlib,
        func = x$func,
        csvFunc = plot_data_csv,
        res = c(80, 95), # resolution of plots
        pdf.width = 10,
        pdf.height = 8,
        add.watermark = watermark,
        card = x$card,
        parent_session = session
      )
    })
  }) ## end of moduleServer
}
