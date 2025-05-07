##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

dataview_plot_expression_ui <- function(
    id,
    label = "",
    height,
    title,
    caption,
    info.text) {
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
    caption = caption,
    options = options,
    outputFunc = plotly::plotlyOutput,
    outputFunc2 = plotly::plotlyOutput,
    info.text = info.text,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    ns_parent = ns,
    editor = TRUE,
    plot_type = "barplot"
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

    plot_data <- shiny::reactive({
      
      shiny::req(pgx$X)
      shiny::req(r.gene(), r.data_type())
      shiny::req(all(colnames(pgx$X) == rownames(pgx$samples)))
      
      gene <- r.gene()
      samples <- r.samples()
      data_type <- r.data_type()
      groupby <- r.data_groupby()
      xgenes <- rownames(pgx$X)

      ##--------------------------------
      ## metadata_vars <- colnames(pgx$samples)
      ## new.options <- tagList(
      ##   shiny::radioButtons(
      ##     inputId = ns("colorby_var"),
      ##     label = "Annotate by:",
      ##     choices = metadata_vars,
      ##     selected = metadata_vars[1],
      ##     inline = FALSE
      ##   )
      ## )
      ##---------------------------------
      
      if (samples[1] == "") return(NULL)
      if (gene == "") return(NULL)
      if (!gene %in% xgenes) return(NULL)

      grpvar <- 1
      grp <- rep(NA, length(samples))
      if (groupby != "<ungrouped>")
        grp <- factor(as.character(pgx$samples[samples, groupby]))
 
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
      df <- data.frame(x = gx, samples = samples, group = grp)
      kk <- which(colnames(pgx$samples) != groupby)
      df <- cbind(df, pgx$samples[samples, kk, drop = FALSE])
      df <- df[, unique(colnames(df)), drop = FALSE]
      pd <- list(
        df = df,
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

    output$rank_list <- renderUI({
      pd <- plot_data()
      shiny::req(pd)
      if (pd$groupby != "<ungrouped>") {
        sortable::bucket_list(
          header = NULL,
          class = "default-sortable custom-sortable",
          sortable::add_rank_list(
            input_id = session$ns("rank_list_basic"),
            text = NULL,
            labels = unique(pd[["df"]]$group),
          )
        )
      } else {
        sortable::bucket_list(
          header = NULL,
          class = "default-sortable custom-sortable",
          sortable::add_rank_list(
            input_id = session$ns("rank_list_basic"),
            text = NULL,
            labels = unique(pd[["df"]]$samples),
          )
        )
      }
    })

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

        data_mean <- tapply(df$x, df$group, mean)
        data_sd <- tapply(df$x, df$group, sd)
        data <- data.frame(group = names(data_mean), mean = data_mean, sd = data_sd)
        
        if (!is.null(input$bars_order)) {
          if (input$bars_order == "ascending") {
            data$group <- reorder(data$group, data$mean)
            df$group <- reorder(df$group, df$x)
          } else if (input$bars_order == "descending") {
            data$group <- reorder(data$group, -data$mean)
            df$group <- reorder(df$group, -df$x)
          } else if (input$bars_order == "custom" && !is.null(input$rank_list_basic) && 
                    all(input$rank_list_basic %in% unique(data$group))) {
            data$group <- factor(data$group, levels = valid_ranks)
            df$group <- factor(df$group, levels = valid_ranks) 
          }
        }

        if (pd$geneplot_type == "barplot") {
          fig <- plotly::plot_ly(
            data = data,
            x = ~group, y = ~mean, type = "bar", name = pd$gene,
            error_y = ~ list(array = sd, color = "#000000"),
            marker = list(color = input$bar_color)
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
        if (!is.null(input$bars_order)) {
          if (input$bars_order == "ascending") {
            df <- df[order(df$x), ]
          } else if (input$bars_order == "descending") {
            df <- df[order(-df$x), ]
          } else if (input$bars_order == "custom" && !is.null(input$rank_list_basic) && 
                    all(input$rank_list_basic %in% unique(df$samples))) {
            df$samples <- factor(df$samples, levels = input$rank_list_basic)
            df <- df[order(df$samples), ]
          }
          df$samples <- factor(df$samples, levels = df$samples)
        }

        metadata <- df[, !colnames(df) %in% c("x","samples","group")]
        df$metadata <- apply(metadata, 1, function(x) {
          paste0(paste0("<b>", colnames(metadata), ": </b>", x, "<br>", collapse = ""))
        })

        fig <- plotly::plot_ly(
          df,
          x = ~samples,
          y = ~x,
          type = "bar",
          name = pd$gene,
          marker = list(color = input$bar_color),
          text = ~metadata,
          textposition = "none",
          hovertemplate = paste0(
            "<b>Sample: </b>%{x}<br>",
            "<b>%{yaxis.title.text}: </b>%{y:.2f}<br>",
            "%{text}<extra></extra>"
          )
        )
        pd$groupby <- ""
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
      download.fmt = c("png", "pdf", "csv", "obj", "svg"),
      res = c(90, 170) * 1, ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
