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
connectivity_plot_scatterPlot_ui <- function(id,
                                               label = "",
                                               fullH = 750) {
  ns <- shiny::NS(id)
  info_text <- strwrap(
    "The <strong>FC-FC scatter plot</strong> provides a pairwise scatterplot of
    logFC fold-change profiles for the selected contrasts. The main
    purpose of this panel is to identify similarity or dissimilarity between
    selected contrasts. The scatter plot is interactive and shows information of
    each gene with a mouse hover-over."
  )
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
    title = "FC-FC scatterplot",
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
connectivity_plot_scatterPlot_server <- function(id,
                                               pgx,
                                               r_contrast,
                                               r_sigdb,
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
          pgx = pgx,
          contrast = r_contrast(),
          sigdb = r_sigdb()
        )
        return(res)
      })

      plot_RENDER <- shiny::reactive({
        res <- plot_data()
        pgx <- res$pgx
        contrast <- res$contrast
        sigdb <- res$sigdb
        
        shiny::req(pgx, contrast)
        shiny::req(sigdb)

        all.ct <- getConnectivityContrasts(sigdb)
        ct1 <- all.ct[1]
        fc1 <- pgx$gx.meta$meta[[1]]$meta.fx
        names(fc1) <- rownames(pgx$gx.meta$meta[[1]])

        fc1 <- getCurrentContrast()$fc
        ct1 <- getCurrentContrast()$name

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
        logfc <- as.numeric(input$logFC)
        shiny::req(logfc)
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
        label.text <- playbase::shortstring(label.text0, 30)
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
          data = df[, 1:2], x = df[, 1], y = df[, 2],
          type = "scattergl", mode = "markers",
          source = "cmapSPLOM", key = rownames(df),
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
