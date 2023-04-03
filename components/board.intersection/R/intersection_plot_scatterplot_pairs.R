##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

intersection_scatterplot_pairs_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)

  info_text <- "For the selected contrasts, the <strong>Pairs</strong> panel provides pairwise scatterplots for the differential expression profiles corresponding to multiple contrasts. The main purpose of this panel is to identify similarity or dissimilarity between selected contrasts. When K >= 3 contrasts are selected, the figure shows a KxK scatterplot matrix. When K <= 2, the Pairs panel provides an interactive pairwise scatterplots for the differential expression profiles of the two selected contrasts. The pairs plot is interactive and shows information of each gene with a mouse hover-over. Users can also select a number points by selecting points with the mouse, using the box selection or the lasso selection tool. Note that the selected genes will appear in input panel on the left sidebar as '<custom>' selection."

  scatterplot_pairs.opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(ns("splom_highlight"), "Highlight genes", TRUE),
      "Enable highlighting genes on the plots. Users can highlight points by selecting them with the mouse, using the box selection or the lasso selection tool."
    )
  )

  PlotModuleUI(
    ns("scatterplot"),
    plotlib = "plotly",
    title = "Scatterplot pairs",
    label = "a",
    info.text = info_text,
    options = scatterplot_pairs.opts,
    download.fmt = c("png", "pdf", "csv"),
    height = c(740, 750),
    width = c("100%", 1000)
  )
}


intersection_scatterplot_pairs_server <- function(id,
                                                  getActiveFoldChangeMatrix,
                                                  level,
                                                  pgx,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      res <- getActiveFoldChangeMatrix()
      fc0 <- res$fc.full
      fc1 <- res$fc

      if (is.null(res)) {
        return(NULL)
      }

      ## fc0 <- fc0[order(-rowSums(fc0**2)),]
      fc0 <- fc0[order(-apply(abs(fc0), 1, max)), ]
      fc0 <- fc0[order(-rowMeans(abs(fc0**2))), ]

      ## selected genes
      sel.genes <- grep("^CD", rownames(fc0), value = TRUE)
      sel.genes <- head(rownames(fc0), 100) ## top100
      sel.genes <- intersect(rownames(fc0), rownames(fc1))
      sel.genes <- intersect(sel.genes, rownames(fc0))
      head(sel.genes)

      ## subsample for speed: take top1000 + 1000
      df <- data.frame(fc0)
      if (1) {
        ntop <- 99999
        ## ntop <- input$splom_ntop
        jj <- match(sel.genes, rownames(df))
        jj <- c(jj, 1:min(ntop, nrow(df)))
        if (nrow(df) > ntop) {
          nremain <- setdiff(1:nrow(df), jj)
          jj <- c(jj, sample(nremain, min(1000, length(nremain)))) ## add 1000 random
        }
        jj <- unique(jj)
        ## df <- data.frame(head(fc0,ntop))
        df <- data.frame(df[jj, ])
        list(df, sel.genes)
      }
    })

    scatterPlotMatrix.PLOT <- function() {
      data <- plot_data()
      df <- data[[1]]
      sel.genes <- data[[2]]
      ## resort selection so that selected genes are drawn last to avoid
      ## covering them up.
      is.sel <- (rownames(df) %in% sel.genes)
      if (input$splom_highlight) {
        df.color <- c("#00000033", "#0066FF")[1 + is.sel]
        df.color <- c("#AAAAAA", "#1e60BB")[1 + is.sel]
        df.color <- c("#CCCCCC22", "#1e60BB88")[1 + is.sel]
      } else {
        df.color <- rep("#00000088", nrow(df))
        df.color <- rep("#1e60BB88", nrow(df))
      }

      ## Labels for top 50
      label.text <- label.text0 <- head(rownames(df)[which(is.sel)], 50)
      label.text <- sub(".*[:]", "", label.text) ## strip prefix??
      label.text <- shortstring(label.text, 30)
      if (sum(is.na(label.text))) label.text[is.na(label.text)] <- ""

      ## reorder so the selected genes don't get overlapped
      jj <- order(is.sel)
      df <- df[jj, ]
      df.color <- df.color[jj]
      sel1 <- match(label.text0, rownames(df)) ## index for labeled

      ## Tooltip text for all
      tt <- rownames(df) ## strip prefix
      ## tt <- sub("","",tt)  ## strip prefix??
      # if(input$level == "gene") {
      if (level == "gene") {
        g <- rownames(df)
        tt <- paste0("<b>", g, "</b> ", pgx$genes[g, "gene_title"])
      }
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

      if (ncol(df) <= 2) {
        ## ----------------------------------------------------
        ## Single pairs plot
        ## ----------------------------------------------------

        rho <- cor(df[, 1], df[, 2])
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
          source = "splom", key = rownames(df),
          ## type = 'scatter', mode = 'markers',
          text = tt,
          hovertemplate = paste0("<br>%{text}<br>x: %{x}<br>y: %{y}<extra></extra>"),
          marker = list(
            color = df.color,
            size = 8,
            line = list(
              width = 0.3,
              ## color = 'rgb(230,230,230)'
              color = "rgb(0,0,0)"
            )
          )
        ) %>%
          plotly::add_annotations(
            x = df[sel1, 1],
            y = df[sel1, 2],
            text = as.character(label.text),
            ## text = rep("mylabel",length(sel1)),
            ## xanchor = 'left',
            xanchor = "center",
            yanchor = "top",
            font = list(size = 14),
            xref = "x",
            yref = "y",
            showarrow = FALSE,
            ax = 20,
            ay = -40
          ) %>%
          plotly::layout(
            ## title= 'Scatterplot',
            annotations = annot.rho,
            hovermode = "closest",
            dragmode = "select",
            ## plot_bgcolor='rgba(240,240,240, 0.95)',
            ## template = "plotly_dark",
            xaxis = c(title = paste(colnames(df)[1], "   (logFC)"), axis),
            yaxis = c(title = paste(colnames(df)[2], "   (logFC)"), axis)
          )
      } else {
        ## ----------------------------------------------------
        ## Scatter pairs matrix plot
        ## ----------------------------------------------------

        dimensions <- lapply(colnames(df), function(a) list(label = a, values = df[, a]))

        ## compute correlations
        rho <- cor(df)
        rho.text <- paste("r=", as.vector(round(rho, digits = 3)))
        n <- ncol(df)

        ## annotation positions (approximated by eye...)
        xann <- 1.02 * (as.vector(mapply(rep, seq(0, 0.98, 1 / n), n)) + 0.05 * 1 / n)
        ## xann = as.vector(mapply(rep,seq(0,1,1/(n-1)),n))
        yann <- 1.08 * (as.vector(rep(seq(1, 0.02, -1 / n), n)) - 0.15 * 1 / n - 0.04)
        ## yann = as.vector(rep(seq(1,0.0,-1/(n-1)),n))

        p <- plotly::plot_ly(df, source = "splom", key = rownames(df)) %>%
          plotly::add_trace(
            type = "splom",
            dimensions = dimensions,
            text = tt,
            hovertemplate = paste0("<br>%{text}<br>x: %{x}<br>y: %{y}<extra></extra>"),
            marker = list(
              color = df.color,
              ## colorscale = pl_colorscale,
              size = 5,
              line = list(
                width = 0.3,
                ## color = 'rgb(230,230,230)'
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
            ## title= 'Scatterplot matrix',
            hovermode = "closest",
            dragmode = "select",
            ## annotations = annot,
            ## plot_bgcolor='rgba(240,240,240, 0.95)',
            ## template = "plotly_dark",
            xaxis = c(domain = NULL, axis),
            yaxis = c(domain = NULL, axis),
            xaxis2 = axis, xaxis3 = axis, xaxis4 = axis, xaxis5 = axis, xaxis6 = axis, xaxis7 = axis,
            yaxis2 = axis, yaxis3 = axis, yaxis4 = axis, yaxis5 = axis, yaxis6 = axis, yaxis7 = axis
          )
        ## %>% plotly::style(diagonal = list(visible = F))
      }

      p <- p %>%
        plotly::layout(margin = list(80, 80, 80, 80)) ## l,r,b,t

      p <- p %>%
        ## config(displayModeBar = FALSE) %>% ## disable buttons
        plotly::config(modeBarButtonsToRemove = setdiff(all.plotly.buttons, "toImage")) %>%
        plotly::config(toImageButtonOptions = list(
          format = "svg",
          height = 800, width = 800, scale = 1.1
        )) %>%
        plotly::config(displaylogo = FALSE) %>%
        plotly::event_register("plotly_selected")

      p
    }

    # observeEvent(input$test, {browser()})

    PlotModuleServer(
      "scatterplot",
      plotlib = "plotly",
      func = scatterPlotMatrix.PLOT,
      csvFunc = plot_data,
      res = 95,
      pdf.width = 5, pdf.height = 5,
      add.watermark = watermark
    )
  })
}
