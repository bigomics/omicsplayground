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
    ),
    shiny::tags$div(
      style = "margin-top: 15px;",
      withTooltip(
        shiny::checkboxInput(
          inputId = ns("show_imputed_values"),
          label = "Show imputed sample values",
          value = FALSE
        ),
        "Show samples (if any) in which the selected feature was imputed.",
        placement = "top"
      )
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
        gx[which(is.na(gx))] <- 0
        ylab <- tspan("Counts", js = FALSE)
      } else if (data_type %in% c("logCPM", "log2")) {
        gx <- pgx$X[pp, samples]
        gx[which(is.na(gx))] <- 0
        ylab <- tspan("Counts (log2)", js = FALSE)
      }
      
      pd <- list(
        df = data.frame(x = gx, samples = samples, group = grp),
        geneplot_type = input$geneplot_type,
        groupby = groupby,
        ylab = ylab,
        gene = gene
      )

      return(pd)

    })

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

      points.color <- rep("black", nrow(df))
      names(points.color) <- df$sample
      if (input$show_imputed_values) {
        counts.values <- pgx$counts[unique(pd$gene), df$samples]
        x.values <- pgx$X[unique(pd$gene), df$samples]
        if (any(is.na(counts.values)) && !any(is.na(x.values)))
          points.color[which(is.na(counts.values))] <- "#C1C1C1"
      }

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

        points.color <- points.color[match(df$samples, names(points.color))]
        df$points.color <- unname(points.color)

        if (pd$geneplot_type == "barplot") {
          fig <- plotly::plot_ly(
            data = data, x = ~group, y = ~mean, type = "bar",
            name = pd$gene, error_y = ~ list(array = sd, color = "#000000"),
            marker = list(color = input$bar_color)
          )          
          fig <- fig %>% plotly::add_markers(
            x = df$group, y = df$x,
            type = "scatter", showlegend = FALSE,
            marker = list(color = ~points.color, size = 8)
          )
          fig
        } else if (pd$geneplot_type == "violin") {
          fig <- df %>%
            plotly::plot_ly(
              x = ~group, y = ~x, split = ~group, type = "violin",
              box = list(visible = TRUE), meanline = list(visible = TRUE),
              points = FALSE, x0 = "", color = ~group,
              colors = omics_pal_d()(length(unique(df$group)))
            )
          fig <- fig %>%
            plotly::add_trace(
              data = df, x = ~group, y = ~x,
              type = 'scatter', mode = 'markers',
              marker = list(color = ~points.color, size = 8),
              showlegend = FALSE, inherit = FALSE
            ) %>%
            plotly::layout(yaxis = list(zeroline = FALSE))
        } else {
          ## boxplot
          fig <- plotly::plot_ly(df, y = ~x, x = ~group, type = "box", boxpoints = FALSE)
          fig <- fig %>% plotly::add_trace(
            data = df, y = ~x, x = ~group,
            type = "scatter", mode = "markers",
            marker = list(color = ~points.color, size = 8),
            showlegend = FALSE
          )
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

        points.color[which(points.color=="black")] <- input$bar_color 
        fig <- plotly::plot_ly(
          df, x = ~samples, y = ~x,
          type = "bar", name = pd$gene,
          marker = list(color = points.color),
          hovertemplate = "<b>Sample: </b>%{x}<br><b>%{yaxis.title.text}:</b> %{y:.2f}<extra></extra>"
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
      fig <- plotly.RENDER() %>% plotly_modal_default()
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
