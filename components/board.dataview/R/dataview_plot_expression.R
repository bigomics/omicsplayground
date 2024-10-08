##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

dataview_plot_expression_ui <- function(
    id,
    label = "",
    height,
    title,
    info
#    caption,
#    info.text
    ) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::radioButtons(
      inputId = ns("geneplot_type"),
      label = "Plot type (only available when {Group by} setting is in use)",
      choices = c("barplot", "box", "violin")
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    info = info,
#    caption = caption,
#    info.text = info.text,
    options = options,
    outputFunc = plotly::plotlyOutput,
    outputFunc2 = plotly::plotlyOutput,
    download.fmt = c("png", "pdf", "csv"),
    height = height,
    show.ai = TRUE
  )
}

dataview_plot_expression_server <- function(id,
                                            pgx,
                                            r.gene = reactive(""),
                                            r.samples = reactive(""),
                                            r.data_type = reactive("counts"),
                                            r.data_groupby = reactive("<ungrouped>"),
                                            info = info,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(r.gene(), r.data_type())
      shiny::req(all(colnames(pgx$X) == rownames(pgx$samples)))
      ## dereference reactives
      gene <- r.gene()
      samples <- r.samples()
      data_type <- r.data_type()
      groupby <- r.data_groupby()
      xgenes <- rownames(pgx$X)

      if (samples[1] == "") {
        return(NULL)
      }
      if (gene == "") {
        return(NULL)
      }
      if (!gene %in% xgenes) {
        return(NULL)
      }

      grpvar <- 1
      grp <- rep(NA, length(samples))
      if (groupby != "<ungrouped>") {
        grp <- factor(as.character(pgx$Y[samples, groupby]))
      }
      # Req to avoid error on dataset change
      shiny::req(length(grp) == length(samples))

      pp <- rownames(pgx$genes)[match(gene, rownames(pgx$genes))]
      if (data_type %in% c("counts", "abundance")) {
        gx <- pgx$counts[pp, samples]
        ylab <- tspan("Counts", js = FALSE)
      } else if (data_type %in% c("logCPM", "log2")) {
        gx <- pgx$X[pp, samples]
        ylab <- tspan("Counts (log2)", js = FALSE)
      }

      # geneplot_type <- "barplot"
      geneplot_type <- input$geneplot_type
      pd <- list(
        df = data.frame(
          x = gx,
          samples = samples,
          group = grp
        ),
        geneplot_type = geneplot_type,
        groupby = groupby,
        ylab = ylab,
        gene = gene
      )
      return(pd)
    })


    plot.RENDER.SAVE <- function() {
      pd <- plot_data()

      shiny::req(pd)

      df <- pd[["df"]]

      par(mar = c(7, 3.5, 2, 1), mgp = c(2.1, 0.8, 0))

      BLUE <- rgb(0.2, 0.5, 0.8, 0.8)
      bee.cex <- ifelse(length(df$x) > 500, 0.1, 0.2)
      bee.cex <- c(0.3, 0.1, 0.05)[cut(length(df$x), c(0, 100, 500, 99999))]

      if (pd$groupby != "<ungrouped>") {
        nnchar <- nchar(paste(unique(df$group), collapse = ""))
        srt <- ifelse(nnchar < 20, 0, 35)
        ngrp <- length(unique(df$group))
        cx1 <- ifelse(ngrp < 10, 1, 0.8)
        cx1 <- ifelse(ngrp > 20, 0.6, cx1)
        if (pd$geneplot_type == "barplot") {
          playbase::gx.b3plot(
            df$x,
            df$group,
            las = 3,
            main = pd$gene,
            ylab = pd$ylab,
            cex.main = 1,
            col.main = "#7f7f7f",
            bar = TRUE,
            border = NA,
            bee.cex = bee.cex,
            xlab = "",
            names.cex = cx1,
            srt = srt,
            col = rgb(0.4, 0.6, 0.85, 0.85)
          )
        } else if (pd$geneplot_type == "violin") {
          playbase::pgx.violinPlot(
            df$x,
            df$group,
            main = pd$gene,
            cex.main = 1,
            xlab = "",
            ylab = ylab,
            vcol = rgb(0.4, 0.6, 0.85, 0.85),
            srt = srt
          )
        } else {
          boxplot(
            df$x ~ df$group,
            main = pd$gene,
            cex.main = 1.0,
            ylab = pd$ylab,
            xlab = "",
            xaxt = "n",
            col = rgb(0.4, 0.6, 0.85, 0.85)
          )
          yy <- sort(unique(df$group))
          text(
            x = 1:length(yy),
            y = par("usr")[3] - 0.03 * diff(range(df$x)),
            labels = yy,
            xpd = NA,
            srt = srt,
            adj = ifelse(srt == 0, 0.5, 0.965),
            cex = cx1
          )
        }
      } else {
        ## plot as bars
        barplot(
          df$x,
          col = BLUE,
          las = 3,
          cex.names = 0.8,
          ylab = pd$ylab,
          xlab = "",
          main = pd$gene,
          cex.main = 1,
          col.main = "#7f7f7f",
          border = NA,
          names.arg = rep(NA, length(df$x))
        )

        ## add labels if needed
        nx <- length(df$x)
        if (nx < 100) {
          cx1 <- ifelse(nx > 20, 0.8, 0.9)
          cx1 <- ifelse(nx > 40, 0.6, cx1)
          cx1 <- ifelse(nx < 10, 1, cx1)
          text(
            x = (1:nx - 0.5) * 1.2,
            y = -0.04 * max(df$x),
            labels = names(df$x),
            las = 3,
            cex = cx1,
            pos = 2,
            adj = 0,
            offset = 0,
            srt = 45,
            xpd = TRUE
          )
        }
      }
    }

    plotly.RENDER <- function() {
      pd <- plot_data()

      shiny::req(pd)

      df <- pd[["df"]]

      BLUE <- omics_colors("brand_blue")
      bee.cex <- ifelse(length(df$x) > 500, 0.1, 0.2)
      bee.cex <- c(0.3, 0.1, 0.05)[cut(length(df$x), c(0, 100, 500, 99999))]

      if (pd$groupby != "<ungrouped>") {
        nnchar <- nchar(paste(unique(df$group), collapse = ""))
        srt <- ifelse(nnchar < 20, 0, 35)
        ngrp <- length(unique(df$group))
        cx1 <- ifelse(ngrp < 10, 1, 0.8)
        cx1 <- ifelse(ngrp > 20, 0.6, cx1)

        if (pd$geneplot_type == "barplot") {
          data_mean <- tapply(df$x, df$group, mean)
          data_sd <- tapply(df$x, df$group, sd)
          data <- data.frame(group = names(data_mean), mean = data_mean, sd = data_sd)

          fig <- plotly::plot_ly(
            data = data,
            x = ~group, y = ~mean, type = "bar", name = pd$gene,
            error_y = ~ list(array = sd, color = "#000000")
          )

          fig <- fig %>% plotly::add_markers(
            x = df$group, y = df$x,
            type = "scatter", showlegend = FALSE,
            marker = list(color = "black", size = 8)
          )
          fig
          ## fig
        } else if (pd$geneplot_type == "violin") {
          fig <- df %>%
            plotly::plot_ly(
              x = ~group,
              y = ~x,
              split = ~group,
              type = "violin",
              box = list(
                visible = TRUE
              ),
              meanline = list(
                visible = TRUE
              ),
              x0 = "",
              color = ~group,
              colors = omics_pal_d()(length(unique(df$group)))
            ) %>%
            plotly::layout(
              yaxis = list(
                zeroline = FALSE
              )
            )
          ## fig
        } else {
          ## boxplot
          fig <- plotly::plot_ly(
            df,
            y = ~x,
            split = ~group,
            boxpoints = "all",
            jitter = 0.3,
            pointpos = 0.0,
            type = "box",
            color = ~group,
            colors = omics_pal_d()(length(unique(df$group)))
          )
          ## fig
        }
      } else {
        ## plot as regular bar plot
        fig <- plotly::plot_ly(
          df,
          x = ~samples,
          y = ~x,
          type = "bar",
          name = pd$gene,
          marker = list(color = BLUE),
          hovertemplate = "<b>Sample: </b>%{x}<br><b>%{yaxis.title.text}:</b> %{y:.2f}<extra></extra>"
        )
        pd$groupby <- ""
        ## fig
      }

      fig <- fig %>%
        plotly::layout(
          xaxis = list(title = "", fixedrange = TRUE),
          yaxis = list(title = pd$ylab, fixedrange = TRUE),
          font = list(family = "Lato"),
          showlegend = FALSE
        ) %>%
        plotly_default()
      fig
    }

    modal_plotly.RENDER <- function() {
      fig <- plotly.RENDER() %>%
        plotly_modal_default()
      fig
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      download.fmt = c("png", "pdf", "csv", "obj"),
      res = c(90, 170) * 1, ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark,
      info = info,
      show.ai = TRUE
    )
  }) ## end of moduleServer
}
