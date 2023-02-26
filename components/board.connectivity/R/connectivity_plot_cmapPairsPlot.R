##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' Importance plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
connectivity_plot_cmapPairsPlot_ui <- function(id,
                                               label = "",
                                               fullH = 750) {
  ns <- shiny::NS(id)
  info_text <- strwrap(
    "The <strong>Pairs</strong> panel provides pairwise scatterplots of
    differential expression profiles for the selected contrasts. The main
    purpose of this panel is to identify similarity or dissimilarity between
    selected contrasts. The pairs plot is interactive and shows information of
    each gene with a mouse hover-over."
  )
  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::selectInput(
        ns("cmap_logFC"),
        "logFC threshold:",
        c(0, 0.5, 1, 2, 3, 4),
        selected = 1
      ),
      "Threshold for (log) foldchange to highlight in plot.",
      placement = "right",
      options = list(container = "body")
    )
  )
  PlotModuleUI(ns("plot"),
    title = "Scatterplot matrix (pairs)",
    label = label,
    plotlib = "plotly",
    info.text = info_text,
    options = plot_opts,
    download.fmt = c("html"),
    height = c(fullH - 80, 700),
    width = c("auto", 1000),
  )
}

#' Importance plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
connectivity_plot_cmapPairsPlot_server <- function(id,
                                                   inputData,
                                                   cmap_contrast,
                                                   cmap_sigdb,
                                                   getConnectivityContrasts,
                                                   getCurrentContrast,
                                                   connectivityScoreTable,
                                                   getConnectivityScores,
                                                   getConnectivityMatrix,
                                                   watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
        res <- list(
          pgx = inputData(),
          cmap_contrast = cmap_contrast(),
          cmap_sigdb = cmap_sigdb()
        )
        return(res)
      })

      plot_RENDER <- shiny::reactive({
        res <- plot_data()
        pgx <- res$pgx
        cmap_contrast <- res$cmap_contrast
        cmap_sigdb <- res$cmap_sigdb

        # pgx <- inputData()
        shiny::req(pgx, cmap_contrast)
        sigdb <- cmap_sigdb

        all.ct <- getConnectivityContrasts(sigdb)
        ct1 <- all.ct[1]
        fc1 <- pgx$gx.meta$meta[[1]]$meta.fx
        names(fc1) <- rownames(pgx$gx.meta$meta[[1]])

        fc1 <- getCurrentContrast()$fc
        ct1 <- getCurrentContrast()$name

        sigdb <- cmap_sigdb
        shiny::req(sigdb)
        ct2 <- all.ct[1]
        sel.row <- connectivityScoreTable$rows_selected()
        if (is.null(sel.row)) {
          return(
            plotly::plotly_empty(type = "scatter", mode = "markers") %>%
              plotly::config(
                displayModeBar = FALSE
              ) %>%
              plotly::layout(
                title = list(
                  text = "Select dataset/contrast",
                  yref = "paper",
                  y = 0.5
                )
              )
          )
        }
        df <- getConnectivityScores()
        df <- df[abs(df$score) > 0, , drop = FALSE]
        ct2 <- rownames(df)[sel.row]
        fc2 <- getConnectivityMatrix(sigdb, select = ct2)[, 1]

        ## match with selection filter
        gg <- unique(c(names(fc1), names(fc2))) ## union or intersection??
        fc1 <- fc1[match(gg, names(fc1))]
        fc2 <- fc2[match(gg, names(fc2))]
        df <- data.frame(fc2, fc1)
        rownames(df) <- gg
        colnames(df) <- c(ct2, ct1)
        df[is.na(df)] <- 0 ## missing as zero???
        df <- df[order(-rowMeans(abs(df**2), na.rm = TRUE)), ]

        ## Number of selected genes
        sel.genes <- grep("^CD", rownames(df), value = TRUE)
        logfc <- as.numeric(input$cmap_logFC)
        sel.genes <- rownames(df)[rowSums(abs(df) > logfc) >= 1] ## minimum FC

        not.sel <- setdiff(rownames(df), sel.genes)
        if (length(not.sel) > 0) {
          nr <- 1000
          ii <- c(sel.genes, sample(not.sel, nr, replace = TRUE))
          df <- df[unique(ii), , drop = FALSE]
        }

        ## remove NA genes from selection
        na.fc <- rownames(df)[rowSums(is.na(df)) > 0] ## probably was missing
        na.zero <- rownames(df)[rowSums(df == 0) > 0] ## probably was missing
        sel.genes <- setdiff(sel.genes, c(na.fc, na.zero))

        ## resort selection so that selected genes are drawn last to avoid
        ## covering them up.
        is.sel <- (rownames(df) %in% sel.genes)
        highlight <- TRUE
        if (highlight) {
          df.color <- c("#00000033", "#0066FF")[1 + is.sel]
          df.color <- c("#AAAAAA55", "#1e60BB88")[1 + is.sel]
        } else {
          df.color <- rep("#77777788", nrow(df))
          df.color <- rep("#1e60BB88", nrow(df))
        }

        ## Labels for top 50
        label.text0 <- head(rownames(df)[which(is.sel)], 50)
        label.text <- shortstring(label.text0, 30)
        if (sum(is.na(label.text))) label.text[is.na(label.text)] <- ""

        ## reorder so the selected genes don't get overlapped
        jj <- order(is.sel)
        df <- df[jj, ]
        df.color <- df.color[jj]
        sel1 <- match(label.text0, rownames(df)) ## index for labeled

        ## Tooltip text for all
        gg <- rownames(df)
        tt <- paste0("<b>", gg, "</b> ", pgx$genes[gg, "gene_title"])
        tt <- gsub("_", " ", tt)
        tt <- sapply(tt, breakstring2, 50, brk = "<br>")

        ## plotly
        ##
        axis <- list(
          showline = TRUE,
          zeroline = TRUE,
          gridcolor = "#dddf",
          ticklen = 4
        )

        if (ncol(df) >= 3) {
          dimensions <- lapply(colnames(df), function(a) list(label = a, values = df[, a]))

          ## compute correlations
          rho <- cor(df, use = "pairwise")
          rho.text <- paste("r=", as.vector(round(rho, digits = 3)))
          n <- ncol(df)

          ## annotation positions (approximated by eye...)
          xann <- 1.02 * (as.vector(mapply(rep, seq(0, 0.98, 1 / n), n)) + 0.05 * 1 / n)
          yann <- 1.08 * (as.vector(rep(seq(1, 0.02, -1 / n), n)) - 0.15 * 1 / n - 0.04)

          p <- plotly::plot_ly(df, source = "cmapSPLOM", key = rownames(df)) %>%
            plotly::add_trace(
              type = "splom",
              dimensions = dimensions,
              text = tt,
              hovertemplate = paste0("<br>%{text}<br>x: %{x}<br>y: %{y}<extra></extra>"),
              marker = list(
                color = df.color,
                size = 5,
                line = list(
                  width = 0.3,
                  color = "rgb(0,0,0)"
                )
              )
            ) %>%
            plotly::add_annotations(
              x = xann,
              y = yann,
              text = rho.text,
              font = list(size = 11),
              xanchor = "left",
              align = "left",
              showarrow = FALSE,
              xref = "paper",
              yref = "paper",
              borderpad = 3,
              bordercolor = "black",
              borderwidth = 0.6
            ) %>%
            plotly::layout(
              hovermode = "closest",
              dragmode = "select",
              xaxis = c(domain = NULL, axis),
              yaxis = c(domain = NULL, axis),
              xaxis2 = axis, xaxis3 = axis, xaxis4 = axis, xaxis5 = axis, xaxis6 = axis, xaxis7 = axis,
              yaxis2 = axis, yaxis3 = axis, yaxis4 = axis, yaxis5 = axis, yaxis6 = axis, yaxis7 = axis
            )
        } else {
          rho <- cor(df[, 1], df[, 2], use = "pairwise")
          annot.rho <- list(
            text = paste("r=", round(rho, 4)),
            font = list(size = 14),
            align = "left",
            showarrow = FALSE,
            xref = "paper",
            yref = "paper",
            x = 0.03,
            y = 0.97,
            borderpad = 8,
            bordercolor = "black",
            borderwidth = 0.6
          )

          p <- plotly::plot_ly(
            data = df[, 1:2], x = df[, 1], y = df[, 2],
            type = "scattergl", mode = "markers",
            source = "cmapSPLOM", key = rownames(df),
            text = tt,
            hovertemplate = paste0("<br>%{text}<br>x: %{x}<br>y: %{y}<extra></extra>"),
            marker = list(
              color = df.color,
              size = 8,
              line = list(
                width = 0.3,
                color = "rgb(0,0,0)"
              )
            )
          )
          if (length(sel1) > 0) {
            p <- p %>%
              plotly::add_annotations(
                x = df[sel1, 1],
                y = df[sel1, 2],
                text = as.character(label.text),
                xanchor = "center",
                yanchor = "top",
                font = list(size = 14),
                xref = "x", yref = "y",
                showarrow = FALSE,
                ax = 20, ay = -40
              )
          }

          p <- p %>%
            plotly::layout(
              annotations = annot.rho,
              hovermode = "closest",
              xaxis = c(title = paste(colnames(df)[1], "   (logFC)"), axis),
              yaxis = c(title = paste(colnames(df)[2], "   (logFC)"), axis)
            )
        }

        p <- p %>%
          plotly::config(toImageButtonOptions = list(format = "svg", height = 800, width = 800, scale = 1.1)) %>%
          plotly::event_register("plotly_selected")

        p
      })


      PlotModuleServer(
        "plot",
        plotlib = "plotly",
        func = plot_RENDER,
        csvFunc = plot_data,
        pdf.width = 8, pdf.height = 8,
        res = 95,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
