##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' @export
across_plot_barplot_ui <- function(id,
                                   title,
                                   info.text,
                                   info.methods = NULL,
                                   info.extra_link = NULL,
                                   width,
                                   height) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(
        inputId = ns("multigene_display"),
        label = "Multi-gene display",
        choices = c("Bar plots" = "barplot", "Heatmap" = "heatmap"),
        selected = "barplot",
        inline = TRUE
      ),
      "Choose how to display multiple genes: individual bar plots or a combined heatmap.",
      placement = "top"
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    plotlib = "plotly",
    label = "a",
    options = options,
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width
  )
}

#' @export
across_plot_barplot_server <- function(id,
                                       getPlotData,
                                       watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- reactive({
      df <- getPlotData()
      shiny::validate(shiny::need(
        !is.null(df) && nrow(df) > 0,
        "No data available. Please select genes and click 'Query'."
      ))
      return(df)
    })

    barplot.RENDER <- function() {
      df <- plot_data()

      color_by <- attr(df, "color_by")
      if (is.null(color_by)) color_by <- "dataset"
      legend_title <- tools::toTitleCase(color_by)

      genes <- unique(df$gene)
      n_genes <- length(genes)
      display_mode <- if (is.null(input$multigene_display)) "barplot" else input$multigene_display

      ## Use heatmap for multiple genes when heatmap mode is selected
      if (n_genes > 1 && display_mode == "heatmap") {
        p <- render_expression_heatmap(df, genes, color_by, legend_title)
        return(p)
      }

      if (n_genes == 1) {
        p <- plotly::plot_ly(
          data = df,
          x = ~reorder(sample, count),
          y = ~count,
          color = ~color_group,
          type = "bar",
          colors = "Set2",
          hovertemplate = paste(
            "<b>%{x}</b><br>",
            "Count: %{y:.2f}<br>",
            "<extra></extra>"
          )
        ) %>%
          plotly::layout(
            title = list(text = genes[1], font = list(size = 14)),
            xaxis = list(
              title = "Sample",
              showticklabels = FALSE
            ),
            yaxis = list(title = "Expression"),
            legend = list(orientation = "h", y = -0.15, title = list(text = legend_title)),
            margin = list(l = 60, r = 20, t = 40, b = 40)
          )
      } else {
        plots <- lapply(genes, function(gene) {
          gene_df <- df[df$gene == gene, ]
          plotly::plot_ly(
            data = gene_df,
            x = ~reorder(sample, count),
            y = ~count,
            color = ~color_group,
            type = "bar",
            colors = "Set2",
            showlegend = (gene == genes[1]),
            hovertemplate = paste(
              "<b>%{x}</b><br>",
              "Count: %{y:.2f}<br>",
              "<extra></extra>"
            )
          ) %>%
            plotly::layout(
              annotations = list(
                list(
                  text = gene,
                  xref = "paper", yref = "paper",
                  x = 0.5, y = 1.05,
                  showarrow = FALSE,
                  font = list(size = 12, weight = "bold")
                )
              ),
              xaxis = list(showticklabels = FALSE)
            )
        })

        p <- plotly::subplot(
          plots,
          nrows = n_genes,
          shareX = FALSE,
          shareY = FALSE,
          titleY = TRUE
        ) %>%
          plotly::layout(
            legend = list(orientation = "h", y = -0.05, title = list(text = legend_title)),
            margin = list(l = 60, r = 20, t = 20, b = 40)
          )
      }

      p <- p %>%
        plotly::config(
          modeBarButtonsToRemove = setdiff(all.plotly.buttons, "toImage"),
          displaylogo = FALSE
        )

      return(p)
    }

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = barplot.RENDER,
      csvFunc = plot_data,
      res = c(85, 100),
      pdf.width = 12, pdf.height = 6,
      add.watermark = watermark
    )
  })
}

#' Render heatmap for multi-gene expression data using playbase::pgx.splitHeatmapFromMatrix
#' @keywords internal
render_expression_heatmap <- function(df, genes, color_by, legend_title) {
  ## Create expression matrix (genes x samples)
  samples <- unique(df$sample)
  expr_matrix <- matrix(NA, nrow = length(genes), ncol = length(samples),
                        dimnames = list(genes, samples))

  for (i in seq_len(nrow(df))) {
    gene <- df$gene[i]
    sample <- df$sample[i]
    if (gene %in% genes && sample %in% samples) {
      expr_matrix[gene, sample] <- df$count[i]
    }
  }

  ## Get sample annotations for color grouping
  sample_info <- df[!duplicated(df$sample), c("sample", "color_group")]
  rownames(sample_info) <- sample_info$sample
  sample_info <- sample_info[samples, , drop = FALSE]

  ## Create annotation data frame
  annot <- data.frame(row.names = samples)
  annot[[legend_title]] <- sample_info$color_group

  ## Order samples by group for better visualization
  sample_order <- order(sample_info$color_group)
  expr_matrix <- expr_matrix[, sample_order, drop = FALSE]
  annot <- annot[sample_order, , drop = FALSE]

  ## Use playbase function to create heatmap
  ## scale = "none" since data may already be transformed (log2 or z-score)
  plt <- playbase::pgx.splitHeatmapFromMatrix(
    X = expr_matrix,
    annot = annot,
    splitx = annot[[legend_title]],
    scale = "none",
    row_clust = (nrow(expr_matrix) > 1),
    show_legend = TRUE,
    rowcex = 1,
    colcex = 0  ## hide column labels (too many samples)
  )

 ## Convert iheatmapr to plotly widget
  plt <- plt %>%
    iheatmapr::to_plotly_list() %>%
    plotly::as_widget()

  plt <- plt %>%
    plotly::config(
      modeBarButtonsToRemove = setdiff(all.plotly.buttons, "toImage"),
      displaylogo = FALSE
    )

  return(plt)
}

