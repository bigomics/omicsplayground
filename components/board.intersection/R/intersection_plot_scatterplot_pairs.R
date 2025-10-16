##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

intersection_scatterplot_pairs_ui <- function(
    id,
    title,
    label = "",
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  scatterplot_pairs.opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(ns("splom_highlight"), tspan("Highlight genes"), TRUE),
      "Enable highlighting genes on the plots. Users can highlight points by selecting them with the mouse, using the box selection or the lasso selection tool."
    )
  )

  PlotModuleUI(
    ns("scatterplot"),
    plotlib = "plotly",
    title = title,
    label = "a",
    caption = caption,
    info.text = info.text,
    options = scatterplot_pairs.opts,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width
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
      shiny::req(res)
      fc0 <- res$fc.full
      fc1 <- res$fc
      
      if (is.null(res)) return(NULL)
      
      fc0 <- fc0[order(-apply(abs(fc0), 1, max, na.rm = TRUE)), ]
      fc0 <- fc0[order(-rowMeans(abs(fc0**2), na.rm = TRUE)), ]
      sel.genes <- intersect(rownames(fc0), rownames(fc1))

      ## subsample for speed: take top1000 + 1000
      df <- data.frame(fc0)
      ntop <- 99999
      jj <- match(sel.genes, rownames(df))
      jj <- c(jj, 1:min(ntop, nrow(df)))
      if (nrow(df) > ntop) {
        nremain <- setdiff(1:nrow(df), jj)
        jj <- c(jj, sample(nremain, min(1000, length(nremain)))) ## add 1000 random
      }
      jj <- unique(jj)
      df <- data.frame(df[jj, ])
      qv <- data.frame(res$qv.full[rownames(df), ])
      list(df, qv, sel.genes)
    })

    scatterPlotMatrix.PLOT <- function() {
      data <- plot_data()
      df <- data[[1]]
      qv <- data[[2]]
      sel.genes <- data[[3]]
      shiny::req(sel.genes)

      is.sel <- (rownames(df) %in% sel.genes)
      df.color <- rep(omics_colors("grey"), nrow(df))
      if (input$splom_highlight)
        df.color <- c("#CCCCCC22", omics_colors("grey"))[1 + is.sel]  

      ## Labels for top 50
      label.text <- label.text0 <- head(rownames(df)[which(is.sel)], 50)
      label.text <- sub(".*[:]", "", label.text)
      label.text <- playbase::shortstring(label.text, 30)
      if (sum(is.na(label.text))) label.text[is.na(label.text)] <- ""

      ## reorder so the selected genes don't get overlapped
      jj <- order(is.sel)
      df <- df[jj, ]
      df.color <- df.color[jj]
      sel1 <- match(label.text0, rownames(df)) ## index for labeled

      ## Tooltip text for all
      tt <- rownames(df) ## strip prefix
      if (level() == "gene") {
        g <- rownames(df)
        tt <- paste0("<b>", g, "</b> ", pgx$genes[g, "gene_title"])
      }
      tt <- gsub("_", " ", tt)
      tt <- sapply(tt, playbase::breakstring2, 50, brk = "<br>")

      ## plotly
      axis <- list(
        showline = TRUE,
        zeroline = TRUE,
        gridcolor = "#dddf",
        ticklen = 4
      )

      if (ncol(df) <= 2) {
        
        rho <- cor.test(df[, 1], df[, 2], use = "pairwise")
        rho.coeff <- round(rho$estimate, 2)
        rho.pv <- paste0("\np = ", round(rho$p.value, 2))
        if (is.na(rho$p.value)) rho.pv = ""
        rho.text <- paste0("r = ", rho.coeff, rho.pv)

        df.color1 <- df.color
        sig.fc <- apply(df, 1, function(x) sum(abs(x)>=1) == 2)
        sig.qv <- apply(qv, 1, function(x) sum(x<=0.05) == 2)
        jj <- which(sig.fc & sig.qv)
        if (any(jj)) df.color1[jj] <- omics_colors("green")

        jj1 <- abs(df[,1])>=1 & qv[,1]<=0.05
        jj2 <- abs(df[,2])>=1 & qv[,2]<=0.05
        jj3 <- unique(c(which(jj1 & !jj2), which(!jj1 & jj2)))
        if (any(jj3)) df.color1[jj3] <- omics_colors("orange")
        
        annot.rho <- list(
          text = rho.text,
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
          text = tt,
          hovertemplate = paste0("<br>%{text}<br>x: %{x}<br>y: %{y}<extra></extra>"),
          marker = list(
            color = df.color1,
            size = 8,
            line = list(
              width = 0.3,
              color = "rgb(0,0,0)"
            )
          )
        ) %>%
          plotly::add_annotations(
            x = df[sel1, 1],
            y = df[sel1, 2],
            text = as.character(label.text),
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
            annotations = annot.rho,
            hovermode = "closest",
            dragmode = "select",
            xaxis = c(title = paste(colnames(df)[1], " (log2FC)"), axis),
            yaxis = c(title = paste(colnames(df)[2], " (log2FC)"), axis)
          )
      } else {

        dimensions <- lapply(colnames(df), function(a) list(label = a, values = df[, a]))

        rho <- Hmisc::rcorr(as.matrix(df))
        rho.coeff <- as.vector(round(rho$r, 2))
        rho.pv <- as.vector(round(rho$P, 2))
        jj <- is.na(rho.pv)
        rho.pv <- paste0("\np = ", rho.pv) 
        rho.pv[jj] <- ""
        rho.text <- paste0("r = ", rho.coeff, rho.pv)
        n <- ncol(df)
        
        xann <- 1.02 * (as.vector(mapply(rep, seq(0, 0.98, 1 / n), n)) + 0.05 * 1 / n)
        yann <- 1.08 * (as.vector(rep(seq(1, 0.02, -1 / n), n)) - 0.15 * 1 / n - 0.04)

        i=1; cols_list=list()
        for(i in 1:ncol(df)) {
          t=1;
          for(t in 1:ncol(df)) {
            comp <- paste(colnames(df)[i], colnames(df)[t], sep = "--VS--")
            df1 <- df[, c(colnames(df)[i], colnames(df)[t]), drop = FALSE]
            qv1 <- qv[, c(colnames(df)[i], colnames(df)[t]), drop = FALSE]
            sig.fc <- apply(df1, 1, function(x) sum(abs(x)>=1) == 2)
            sig.qv <- apply(qv1, 1, function(x) sum(x<=0.05) == 2)
            cols_list[[comp]] <- df.color
            jj <- which(sig.fc & sig.qv)
            if (any(jj)) cols_list[[comp]][jj] <- omics_colors("green")
            jj1 <- abs(df1[,1])>=1 & qv1[,1]<=0.05
            jj2 <- abs(df1[,2])>=1 & qv1[,2]<=0.05
            jj3 <- unique(c(which(jj1 & !jj2), which(!jj1 & jj2)))
            if (any(jj3)) cols_list[[comp]][jj3] <- omics_colors("orange")
          }
        }

        cols <- unname(do.call(c, cols_list))
        
        p <- plotly::plot_ly(df, source = "splom", key = rownames(df)) %>%
          plotly::add_trace(
            type = "splom",
            dimensions = dimensions,
            text = tt,
            hovertemplate = paste0("<br>%{text}<br>x: %{x}<br>y: %{y}<extra></extra>"),
            marker = list(
              color = cols,
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
      }

      p <- p %>%
        plotly::layout(margin = list(80, 80, 80, 80)) ## l,r,b,t

      p <- p %>%
        plotly::config(modeBarButtonsToRemove = setdiff(all.plotly.buttons, "toImage")) %>%
        plotly::config(toImageButtonOptions = list(
          format = "svg",
          height = 800, width = 800, scale = 1.1
        )) %>%
        plotly::config(displaylogo = FALSE) %>%
        plotly::event_register("plotly_selected")

      p
    }

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
