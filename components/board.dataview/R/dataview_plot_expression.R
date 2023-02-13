##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

dataview_plot_expression_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)

  info_text <- paste0("Expression barplot of grouped samples (or cells) for the gene selected in the <code>Search gene</code> Samples (or cells) in the barplot can be ungrouped by setting the <code>grouped</code> under the main <i>Options</i>.")

  PlotModuleUI(
    ns("pltmod"),
    title = "Abundance/expression",
    label = label,
    outputFunc = plotly::plotlyOutput,
    outputFunc2 = plotly::plotlyOutput,
    info.text = info_text,
    download.fmt = c("png", "pdf", "csv", "obj"),
    ## width = c("auto","100%"),
    height = height
  )
}

dataview_plot_expression_server <- function(id,
                                            pgx,
                                            r.gene = reactive(""),
                                            r.samples = reactive(""),
                                            r.data_type = reactive("counts"),
                                            r.data_groupby = reactive("<ungrouped>"),
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    dbg("[dataview_expressionplot_server] created!")

    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(r.gene(), r.data_type())

      ## dereference reactives
      gene <- r.gene()
      samples <- r.samples()
      data_type <- r.data_type()
      groupby <- r.data_groupby()

      if (samples[1] == "") samples <- colnames(pgx$X)
      if (gene == "") genes <- rownames(pgx$X)[1]

      grpvar <- 1
      grp <- rep(NA, length(samples))
      if (groupby != "<ungrouped>") {
        ## grp  = factor(as.character(pgx$samples[,3]))
        grp <- factor(as.character(pgx$Y[samples, groupby]))
      }

      pp <- rownames(pgx$genes)[match(gene, pgx$genes$gene_name)]
      gx <- NULL
      ylab <- NULL
      if (data_type == "counts") {
        gx <- pgx$counts[pp, samples]
        ylab <- "expression (counts)"
      } else if (data_type == "CPM") {
        gx <- 2**pgx$X[pp, samples]
        ylab <- "expression (CPM)"
      } else if (data_type == "logCPM") {
        gx <- pgx$X[pp, samples]
        ylab <- "expression (log2CPM)"
      }

      geneplot_type <- "barplot"
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
        if (pd$geneplot_type == "bar") {
          gx.b3plot(
            df$x,
            df$group,
            las = 3,
            main = pd$gene,
            ylab = pd$ylab,
            cex.main = 1,
            col.main = "#7f7f7f",
            bar = TRUE,
            border = NA,
            ## bee  =  ifelse(length(df$x) < 500,TRUE,FALSE),
            bee.cex = bee.cex,
            ## sig.stars = TRUE,
            ## max.stars = 5,
            xlab = "",
            names.cex = cx1,
            srt = srt,
            col = rgb(0.4, 0.6, 0.85, 0.85)
          )
        } else if (pd$geneplot_type == "violin") {
          pgx.violinPlot(
            df$x,
            df$group,
            main = pd$gene,
            cex.main = 1,
            xlab = "",
            ylab = ylab,
            ## vcol = rgb(0.2,0.5,0.8,0.8),
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

    plot.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      df <- pd[["df"]]
      ## par(mar=c(7,3.5,2,1), mgp=c(2.1,0.8,0))

      BLUE <- rgb(0.2, 0.5, 0.8, 0.8)
      bee.cex <- ifelse(length(df$x) > 500, 0.1, 0.2)
      bee.cex <- c(0.3, 0.1, 0.05)[cut(length(df$x), c(0, 100, 500, 99999))]

      if (pd$groupby != "<ungrouped>") {
        nnchar <- nchar(paste(unique(df$group), collapse = ""))
        srt <- ifelse(nnchar < 20, 0, 35)
        ngrp <- length(unique(df$group))
        cx1 <- ifelse(ngrp < 10, 1, 0.8)
        cx1 <- ifelse(ngrp > 20, 0.6, cx1)

        if (pd$geneplot_type == "bar") {
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
              x0 = ""
            ) %>%
            plotly::layout(
              yaxis = list(
                ## title = "",
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
            type = "box"
          )
          ## fig
        }
      } else {
        ## plot as regular bar plot
        fig <- plotly::plot_ly(df, x = ~samples, y = ~x, type = "bar", name = pd$gene)
        # fig <- plotly::plot_ly(x = df$samples, y = df$x, type = 'bar', name = pd$gene)
        pd$groupby <- ""
        ## fig
      }

      ## DEBUG
      if (0) {
        message("[DEBUG] saving fig.rda ...")
        y <- array(c(1, 2, 3, 4, 5, 6), dimnames = list(letters[1:6]))
        ## fig2 <- plotly::plot_ly(x = names(y), y = y, type = "bar")
        # fig2 <- plotly::plot_ly(pd$df,  x = ~samples, y = ~x, type = "bar")
        ## fig2 <- plotly::plot_ly(x = pd$df$samples, y = pd$df$x, type = "bar")
        fig2 <- plotly::plot_ly(x = df$samples, y = df$x, type = "bar")
        saveRDS(fig2, file = "/tmp/figplotexpression-render-debug.rds")
        message("[DEBUG] done saving!")
      }

      fig <- fig %>%
        plotly::layout(
          xaxis = list(title = pd$groupby, fixedrange = TRUE),
          yaxis = list(title = pd$ylab, fixedrange = TRUE),
          font = list(family = "Lato"),
          showlegend = FALSE
          ## title = pd$gene
        ) %>%
        ## plotly::config(displayModeBar = FALSE) %>%
        plotly::config(displaylogo = FALSE) %>%
        plotly::config(
          modeBarButtons = list(list("toImage")),
          toImageButtonOptions = list(format = "svg", height = 500, width = 900)
        )
      fig
    }

    modal_plot.RENDER <- function() {
      fig <- plot.RENDER() %>%
        plotly::layout(
          showlegend = TRUE,
          font = list(
            size = 18
          )
        )
      fig <- plotly::style(fig, marker.size = 14)
      fig
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      plotlib2 = "plotly",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      download.fmt = c("png", "pdf", "csv", "obj"),
      renderFunc = plotly::renderPlotly,
      renderFunc2 = plotly::renderPlotly,
      ## renderFunc = shiny::renderPlot,
      ## renderFunc2 = shiny::renderPlot,
      res = c(90, 170) * 1, ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
