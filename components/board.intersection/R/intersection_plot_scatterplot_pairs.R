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
  width
) {
  ns <- shiny::NS(id)

  scatterplot_pairs.opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("annotate"),
        tspan("Annotate top features"),
        TRUE
      ),
      "Annotate top 50 features"
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("corr_line"),
        tspan("Show correlation (r=1) line"),
        FALSE
      ),
      "Show correlation (r=1) line."
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
    width = width,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "significance"
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

      if (is.null(res)) {
        return(NULL)
      }

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
      qv <- data.frame(res$qv.full)
      cm <- intersect(rownames(df), rownames(qv))
      df <- df[cm, , drop = FALSE]
      qv <- qv[cm, , drop = FALSE]
      list(df, qv, sel.genes)
    })

    scatterPlotMatrix.PLOT <- function() {
      data <- plot_data()
      df <- data[[1]]
      qv <- data[[2]]
      sel.genes <- data[[3]]
      shiny::req(sel.genes)

      ## Editor: significance colors
      clr_ns <- if (!is.null(input$color_ns)) input$color_ns else omics_colors("grey")
      clr_both <- if (!is.null(input$color_both)) input$color_both else omics_colors("green")
      clr_one <- if (!is.null(input$color_one)) input$color_one else omics_colors("orange")

      is.sel <- (rownames(df) %in% sel.genes)
      df.color <- rep(clr_ns, nrow(df))
      ## if (input$splom_highlight)
      ##  df.color <- c("#CCCCCC22", omics_colors("grey"))[1 + is.sel]

      ## Labels for top 50 (or custom from editor)
      if (isTRUE(input$custom_labels) && !is.null(input$label_features) && input$label_features != "") {
        custom_genes <- strsplit(input$label_features, "\\s+")[[1]]
        ## match custom genes to rownames (by symbol or rowname)
        all_symbols <- playbase::probe2symbol(rownames(df), pgx$genes, "gene_name", fill_na = TRUE)
        idx <- which(rownames(df) %in% custom_genes | all_symbols %in% custom_genes)
        label.text0 <- rownames(df)[idx]
        label.text <- all_symbols[idx]
        label.text <- playbase::shortstring(label.text, 30)
      } else {
        label.text <- label.text0 <- head(rownames(df)[which(is.sel)], 50)
        label.text <- sub(".*[:]", "", label.text)
        label.text <- playbase::probe2symbol(label.text, pgx$genes, "gene_name", fill_na = TRUE)
        label.text <- playbase::shortstring(label.text, 30)
      }
      if (sum(is.na(label.text))) label.text[is.na(label.text)] <- ""

      ## reorder so the selected genes don't get overlapped
      jj <- order(is.sel)
      df <- df[jj, ]
      qv <- qv[jj, ]
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
      axis <- list(showline = TRUE, zeroline = TRUE, gridcolor = "#dddf", ticklen = 4)

      if (ncol(df) <= 2) {
        rho <- cor.test(df[, 1], df[, 2], use = "pairwise")
        rho.coeff <- round(rho$estimate, 2)
        rho.pv <- paste0("\np = ", round(rho$p.value, 2))
        if (is.na(rho$p.value)) rho.pv <- ""
        rho.text <- paste0("r = ", rho.coeff, rho.pv)

        df.color1 <- df.color
        sig.fc <- apply(df, 1, function(x) sum(abs(x) >= 1) == 2)
        sig.qv <- apply(qv, 1, function(x) sum(x <= 0.05) == 2)
        jj <- which(sig.fc & sig.qv)
        if (any(jj)) df.color1[jj] <- clr_both

        jj1 <- abs(df[, 1]) >= 1 & qv[, 1] <= 0.05
        jj2 <- abs(df[, 2]) >= 1 & qv[, 2] <= 0.05
        jj3 <- unique(c(which(jj1 & !jj2), which(!jj1 & jj2)))
        if (any(jj3)) df.color1[jj3] <- clr_one

        ## color just selected: dim non-labeled points
        if (isTRUE(input$color_selection) && length(label.text0) > 0) {
          is.labeled <- rownames(df) %in% label.text0
          df.color1[!is.labeled] <- "#DDDDDD"
        }

        ## make non-selected genes transparent
        opacity <- rep(1, nrow(df))
        if (isTRUE(input$color_selection) && length(label.text0) > 0) {
          is.labeled <- rownames(df) %in% label.text0
          opacity[!is.labeled] <- 0.15
        } else if (sum(is.sel) > 0) {
          no.sel <- !rownames(df) %in% sel.genes
          opacity[no.sel] <- 0.1
        }

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
            opacity = opacity,
            line = list(
              width = 0.3,
              color = "rgb(0,0,0)"
            )
          )
        )
        if (input$annotate) {
          p <- p %>%
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
            )
        }
        if (input$corr_line) {
          rng <- range(c(df[, 1], df[, 2]), na.rm = TRUE)
          p <- p %>% plotly::add_lines(x = rng, y = rng,
            line = list(color = "black", dash = "dash", width = 1),
            showlegend = FALSE, inherit = FALSE)
        }
        p <- p %>%
          plotly::layout(
            annotations = annot.rho,
            hovermode = "closest",
            dragmode = "select",
            xaxis = c(title = paste(colnames(df)[1], " (log2FC)"), axis),
            yaxis = c(title = paste(colnames(df)[2], " (log2FC)"), axis)
          )
      } else {
        ctx.comp <- unique(paste0(rep(colnames(df), each = ncol(df)), "--VS--", rep(colnames(df))))
        scale_factor <- 1 / length(ctx.comp) / 3
        scale_factor <- max(min(scale_factor, 1), 0.85)

        plot_list <- list()
        for (i in 1:length(ctx.comp)) {
          c1 <- strsplit(ctx.comp[i], "--VS--")[[1]][1]
          c2 <- strsplit(ctx.comp[i], "--VS--")[[1]][2]
          if (c1 == c2) next
          cc <- unique(c(paste0(c1, "--VS--", c2), paste0(c2, "--VS--", c1)))
          if (any(cc %in% names(plot_list))) next

          df1 <- df[, c(c1, c2), drop = FALSE]
          qv1 <- qv[, c(c1, c2), drop = FALSE]

          qv1 <- qv1[rownames(df1), , drop = FALSE]
          ff <- rownames(df1)
          ff <- paste0("<b>", ff, "</b> ", pgx$genes[ff, "gene_title"])
          ff <- sapply(gsub("_", " ", ff), playbase::breakstring2, 50, brk = "<br>")
          hovertext <- paste0(
            ff, "<br>",
            "x: ", round(df1[, 1], 2), "<br>",
            "y: ", round(df1[, 2], 2)
          )

          rho <- cor.test(df1[, 1], df1[, 2], use = "pairwise")
          rho.coeff <- round(rho$estimate, 2)
          rho.pv <- paste0("\np = ", round(rho$p.value, 2))
          if (is.na(rho$p.value)) rho.pv <- ""
          rho.text <- paste0("r = ", rho.coeff, rho.pv)

          df.color1 <- df.color
          sig.fc <- apply(df1, 1, function(x) sum(abs(x) >= 1) == 2)
          sig.qv <- apply(qv1, 1, function(x) sum(x <= 0.05) == 2)
          jj <- which(sig.fc & sig.qv)
          if (any(jj)) df.color1[jj] <- clr_both

          jj1 <- abs(df1[, 1]) >= 1 & qv1[, 1] <= 0.05
          jj2 <- abs(df1[, 2]) >= 1 & qv1[, 2] <= 0.05
          jj3 <- unique(c(which(jj1 & !jj2), which(!jj1 & jj2)))
          if (any(jj3)) df.color1[jj3] <- clr_one

          ## color just selected: dim non-labeled points
          if (isTRUE(input$color_selection) && length(label.text0) > 0) {
            is.labeled <- rownames(df) %in% label.text0
            df.color1[!is.labeled] <- "#DDDDDD"
          }

          ## make non-selected genes transparent
          opacity <- rep(1, nrow(df))
          if (isTRUE(input$color_selection) && length(label.text0) > 0) {
            is.labeled <- rownames(df) %in% label.text0
            opacity[!is.labeled] <- 0.15
          } else if (sum(is.sel) > 0) {
            no.sel <- !rownames(df) %in% sel.genes
            opacity[no.sel] <- 0.1
          }

          annot.rho <- list(
            text = rho.text,
            font = list(size = 13 * scale_factor),
            align = "left",
            showarrow = FALSE,
            xref = "paper",
            yref = "paper",
            x = 0.02,
            y = 0.98,
            xanchor = "left",
            yanchor = "top"
          )

          p <- plotly::plot_ly(
            data = df1, x = df1[, c1], y = df1[, c2],
            type = "scattergl", mode = "markers",
            marker = list(
              color = df.color1, size = 8 * scale_factor, opacity = opacity,
              line = list(width = 0.3, color = "rgb(0,0,0)")
            ),
            text = hovertext, hoverinfo = "text",
            hovertemplate = "%{text}<extra></extra>"
          )
          if (input$annotate) {
            p <- p %>%
              plotly::add_annotations(
                x = df1[sel1, 1],
                y = df1[sel1, 2],
                text = as.character(label.text),
                xanchor = "center",
                yanchor = "top",
                font = list(size = 14 * scale_factor),
                xref = "x",
                yref = "y",
                showarrow = FALSE,
                ax = 20,
                ay = -40
              )
          }
          if (input$corr_line) {
            rng <- range(c(df1[, 1], df1[, 2]), na.rm = TRUE)
            p <- p %>% plotly::add_lines(
              x = rng, y = rng,
              line = list(color = "black", dash = "dash", width = 2),
              showlegend = FALSE, inherit = FALSE
            )
          }
          p <- p %>%
            plotly::layout(
              annotations = annot.rho,
              hovermode = "closest", dragmode = "select",
              xaxis = list(
                title = list(
                  text = paste(colnames(df1)[1], " (log2FC)"),
                  font = list(size = 14 * scale_factor)
                ),
                showline = TRUE, ticklen = 4
              ),
              yaxis = list(
                title = list(
                  text = paste(colnames(df1)[2], " (log2FC)"),
                  font = list(size = 14 * scale_factor)
                ),
                showline = TRUE, ticklen = 4
              ),
              showlegend = FALSE
            ) %>%
            plotly::layout(margin = list(80, 40, 100, 60)) %>%
            plotly::config(modeBarButtonsToRemove = setdiff(all.plotly.buttons, "toImage")) %>%
            plotly::config(toImageButtonOptions = list(
              format = "svg", height = 800, width = 800, scale = 1.1
            )) %>%
            plotly::config(displaylogo = FALSE) %>%
            plotly::event_register("plotly_selected")

          plot_list[[ctx.comp[i]]] <- p
        }

        nr <- ceiling(length(plot_list) / 2)
        fig <- plotly::subplot(plot_list,
          nrows = nr, shareX = FALSE, shareY = FALSE,
          titleX = TRUE, titleY = TRUE, margin = 0.05
        )
        fig
      }
    }

    PlotModuleServer(
      "scatterplot",
      plotlib = "plotly",
      func = scatterPlotMatrix.PLOT,
      csvFunc = plot_data,
      res = 95,
      pdf.width = 5, pdf.height = 5,
      add.watermark = watermark,
      parent_session = session
    )
  })
}
