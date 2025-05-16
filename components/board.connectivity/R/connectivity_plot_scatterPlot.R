##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
connectivity_plot_scatterPlot_ui <- function(
    id,
    title,
    info.text,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::selectInput(
        ns("logFC"),
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
    title = title,
    label = label,
    plotlib = "plotly",
    caption = caption,
    info.text = info.text,
    options = plot_opts,
    download.fmt = c("pdf", "png", "svg"),
    height = height,
    width = width,
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
connectivity_plot_scatterPlot_server <- function(id,
                                                 pgx,
                                                 r_sigdb,
                                                 getConnectivityContrasts,
                                                 getCurrentContrast,
                                                 connectivityScoreTable,
                                                 getConnectivityScores,
                                                 getConnectivityMatrix,
                                                 watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      ## this returns data to plot. all reactivity should be resolved
      ## and isolated here before doing the plotting
      plot_data <- shiny::reactive({
        sigdb <- r_sigdb()
        logfc <- as.numeric(input$logFC)

        shiny::req(pgx$X)
        shiny::req(sigdb)
        shiny::req(logfc)
        shiny::req(connectivityScoreTable$rows_selected())

        all.ct <- getConnectivityContrasts(sigdb)
        fc1 <- getCurrentContrast()$fc
        ct1 <- getCurrentContrast()$name
        sel.row <- connectivityScoreTable$rows_selected()

        df <- sel.genes <- NULL
        if (!is.null(sel.row)) {
          df <- getConnectivityScores()
          shiny::req(df)
          if (is.null(df)) {
            return(NULL)
          }

          df <- df[abs(df$score) > 0, , drop = FALSE]
          ct2 <- rownames(df)[sel.row]
          fc2 <- getConnectivityMatrix(sigdb, select = ct2)[, 1]

          ## match with selection filter
          gg <- unique(c(names(fc1), names(fc2))) ## union or intersection??
          fc1 <- fc1[match(gg, names(fc1))]
          fc2 <- fc2[match(gg, names(fc2))]
          names(fc1) <- gg
          names(fc2) <- gg
          df <- data.frame(fc2, fc1)
          rownames(df) <- gg
          colnames(df) <- c(ct2, ct1)
          df[is.na(df)] <- 0 ## missing as zero???
          df <- df[order(-rowMeans(abs(df**2), na.rm = TRUE)), ]

          ## Number of selected genes
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
        }

        ## Result object. Remember the first element of the list must
        ## be a dataframe for CSV download from the plotmodule UI.
        res <- list(
          df = df,
          sel.row = sel.row,
          sel.genes = sel.genes,
          pgx = pgx
        )
        return(res)
      })


      ## Plots get data from plot_data. This should be a simple
      ## function, not reactive. Any slow pre-plotting should be in
      #
      plot_RENDER <- function() {
        res <- plot_data()
        pgx <- res$pgx
        sel.row <- res$sel.row
        sel.genes <- res$sel.genes
        df <- res$df

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
        rownames_df <- playbase::probe2symbol(rownames(df), pgx$genes, "gene_name", fill_na = TRUE)
        label.text0 <- head(rownames_df[which(is.sel)], 50)
        label.text <- playbase::shortstring(label.text0, 30)
        if (sum(is.na(label.text))) label.text[is.na(label.text)] <- ""

        ## reorder so the selected genes don't get overlapped
        jj <- order(is.sel)
        df <- df[jj, ]
        df.color <- df.color[jj]
        sel1 <- match(label.text0, playbase::probe2symbol(rownames(df), pgx$genes, "gene_name", fill_na = TRUE)) ## index for labeled

        ## Tooltip text for all
        gg <- playbase::probe2symbol(rownames(df), pgx$genes, "gene_name", fill_na = TRUE)
        tt <- paste0("<b>", gg, "</b> ", playbase::probe2symbol(gg, pgx$genes, "gene_title", fill_na = TRUE))
        tt <- gsub("_", " ", tt)
        tt <- sapply(tt, playbase::breakstring2, 50, brk = "<br>")

        ## plotly
        ##
        axis <- list(
          showline = TRUE,
          zeroline = TRUE,
          gridcolor = "#dddf",
          ticklen = 4
        )

        p <- plotly::plot_ly(
          data = df[, 1:2],
          x = df[, 1],
          y = df[, 2],
          type = "scattergl",
          mode = "markers",
          source = "cmapSPLOM",
          key = rownames(df),
          text = tt,
          hovertemplate = paste0("<br>%{text}<br>x: %{x}<br>y: %{y}<extra></extra>"),
          marker = list(
            color = df.color,
            size = 7,
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

        p <- p %>%
          plotly::layout(
            annotations = annot.rho,
            hovermode = "closest",
            xaxis = c(title = paste(colnames(df)[1], "   (logFC)"), axis),
            yaxis = c(title = paste(colnames(df)[2], "   (logFC)"), axis)
          )

        p <- p %>%
          plotly::config(toImageButtonOptions = list(format = "svg", height = 800, width = 800, scale = 1.1)) %>%
          plotly::event_register("plotly_selected")

        p
      }

      PlotModuleServer(
        "plot",
        plotlib = "plotly",
        func = plot_RENDER,
        csvFunc = plot_data,
        pdf.width = 8,
        pdf.height = 8,
        res = 95,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
