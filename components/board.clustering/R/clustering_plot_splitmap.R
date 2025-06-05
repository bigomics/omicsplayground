##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
    width) {
  ns <- shiny::NS(id)

  splitmap_opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("show_legend"), "show legend", value = TRUE),
      "Show or hide the legend.",
      placement = "right",
      options = list(container = "body")
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("sample_cor"), "show sample-sample correlation", value = FALSE),
      "Show sample-sample correlation heatmap.",
      placement = "right",
      options = list(container = "body")
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
                                            hm_scale,
                                            hm_topmode,
                                            hm_clustk,
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

      sel <- selected_phenotypes()

      sample_cor <- FALSE
      if (input$sample_cor) {
        sample_cor <- TRUE
        zx <- cor(pgx$X, method = "pearson")
        D <- as.dist(1 - zx)
        D[which(is.nan(D) | is.na(D))] <- 1
        hc <- fastcluster::hclust(D, method = "ward.D2")
        zx.idx <- paste0("S", cutree(hc, hm_clustk()))        
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

      if (nrow(zx) <= 1) { return(NULL) }
      show_rownames <- TRUE
      if (nrow(zx) > input$num_rownames) show_rownames <- FALSE

      # Use editor settings for text sizes
      label_cex <- ifelse(!is.null(input$label_size), input$label_size/10, 1)  # normalize to base size 10
      
      scale.mode <- "none"
      if (hm_scale() == "relative") { scale.mode <- "row.center" }
      if (hm_scale() == "BMC") { scale.mode <- "row.bmc" }
      scale.mode

      ## split genes dimension in 5 groups
      splity <- 5
      splity <- 6
      if (!is.null(zx.idx)) { splity <- zx.idx }

      ## split samples
      splitx <- NULL
      splitx <- filt$grp

      # Calculate text sizes
      cex0 <- ifelse(!is.null(splitx) && length(splitx) <= 10, 1.05, 0.85)  # fixed title size scaling
      cex1 <- label_cex * ifelse(ncol(zx) > 200, 0, ifelse(ncol(zx) > 100, 0.5, ifelse(ncol(zx) > 50, 0.75, 1)))
      cex2 <- label_cex * ifelse(nrow(zx) > 60, 0.8, 1)

      shiny::validate(shiny::need(
        !is.null(cex1) && !is.null(cex2),
        "Text size settings must not be empty."
      ))

      # Set column text size to 0 if show_colnames is FALSE
      if (!input$show_colnames) cex1 <- 0
      show_colnames <- (cex1 > 0)

      # Select annot to display (user input)
      sel <- selected_phenotypes()
      if (length(sel) == 0) {
        annot <- NULL
      } else {
        annot <- annot[, sel, drop = FALSE]
      }

      rownames(zx) <- sub("HALLMARK:HALLMARK_", "HALLMARK:", rownames(zx))
      rownames(zx) <- gsub(playdata::GSET_PREFIX_REGEX, "", rownames(zx))
      rownames(zx) <- substring(rownames(zx), 1, 50) ## cut long names...
      if (hm_level() == "gene") {
        rownames(zx) <- sub(".*:", "", rownames(zx))
        rownames(zx) <- playbase::probe2symbol(rownames(zx), pgx$genes, labeltype(), fill_na = TRUE)
      }
      
      if (hm_level() == "geneset") { rownames(zx) <- tolower(rownames(zx)) }

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

      if (sample_cor) {
        splity <- NULL
        scale.mode <- "none"
        cluster_rows <- FALSE
        cluster_columns <- FALSE
        zlim <-  c(min(zx),max(zx))
      } else {
        cluster_rows <- TRUE
        cluster_columns <- TRUE
        zlim <- NULL
      }
      
      playbase::gx.splitmap(
        zx,
        split = splity,
        splitx = splitx,
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
        annot.cex = input$annot_cex/10,
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
      scale <- "none"
      if (hm_scale() == "relative") { scale <- "row.center" }
      if (hm_scale() == "BMC") { scale <- "row.bmc" }

      plt <- NULL
      pd <- plot_data()
      filt <- pd[["filt"]]
      X <- pd[["zx"]]
      annot <- pd[["annot"]]
      splity <- pd[["zx.idx"]]
      sample_cor <- pd[["sample_cor"]]      

      ## sample clustering index
      splitx <- NULL
      splitx <- filt$grp

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
      label_cex <- ifelse(!is.null(input$label_size), input$label_size/10, 1)  # normalize to base size 10
      
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
            "<b>", pgx$genes[g, "gene_name"], "</b>. ",
            pgx$genes[g, "gene_title"], "."
          )
          playbase::breakstring2(aa, 50, brk = "<br>")
        }
        tooltips <- sapply(rownames(X), getInfo)
        labeled_features <- NULL
        rownames(X) <- playbase::probe2symbol(rownames(X), pgx$genes, labeltype(), fill_na = TRUE)
      } else {
        aa <- gsub("_", " ", rownames(X)) ## just geneset names
        tooltips <- sapply(aa, function(x) {
          playbase::breakstring2(x, 50, brk = "<br>")
        })
        names(tooltips) <- rownames(X)
      }
      shiny::showNotification("Rendering iHeatmap...")

      if (sample_cor) idx = NULL else idx = splity
      plt <- playbase::pgx.splitHeatmapFromMatrix(
        X = X,
        annot = annotF,
        ytips = tooltips,
        idx = idx,
        splitx = splitx,
        scale = scale,
        row_annot_width = 0.025,
        rowcex = row_cex,
        colcex = col_cex,
        show_legend = input$show_legend,
        return_x_matrix = TRUE
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
    
    plot_data_csv <- function() { plotly_splitmap.RENDER_get()$X }

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
