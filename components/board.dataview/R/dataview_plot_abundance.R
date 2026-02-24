##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

dataview_plot_abundance_ui <- function(
  id,
  label = "",
  height,
  width,
  title,
  info.text,
  caption
) {
  ns <- shiny::NS(id)

  menu_grouped <- "<code>grouped</code>"

  options <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("show_overall_prop"), "Show overall proportions", FALSE
      ),
      "Show proportions of major feature (gene, protein) types in the whole dataset.",
      placement = "top"
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    caption = caption,
    plotlib = "plotly",
    info.text = info.text,
    options = options,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "grouped_barplot"
  )
}

dataview_plot_abundance_server <- function(id,
                                           getCountsTable,
                                           r.samples = reactive(""),
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      res <- getCountsTable()
      samples <- r.samples()
      shiny::req(res)
      if (!input$show_overall_prop) {
        return(list(prop.counts = res$prop.counts))
      } else {
        return(list(prop.counts = res$prop.counts, gset.genes = res$gset.genes))
      }
    })

    output$rank_list <- shiny::renderUI({
      res <- plot_data()
      shiny::req(res)
      sortable::bucket_list(
        header = NULL,
        class = "default-sortable custom-sortable",
        sortable::add_rank_list(
          input_id = session$ns("rank_list_basic"),
          text = NULL,
          labels = colnames(res$prop.counts)
        )
      )
    })

    output$custom_palette_ui <- shiny::renderUI({
      shiny::req(input$palette == "custom")
      res <- plot_data()
      shiny::req(res)
      genes <- rownames(head(res$prop.counts, 5))
      default_clrs <- omics_pal_d(palette = "expanded")(length(genes))
      pickers <- lapply(seq_along(genes), function(i) {
        colourpicker::colourInput(
          session$ns(paste0("custom_color_", i)),
          label = genes[i],
          value = default_clrs[i]
        )
      })
      shiny::tagList(pickers)
    })

    plotly.RENDER <- function(return_csv = FALSE) {
      res <- plot_data()
      shiny::req(res)

      if (!input$show_overall_prop) {
        long.data <- reshape2::melt(head(res$prop.counts, 5))
        colnames(long.data) <- c("gene", "sample", "value")
        long.data$sample <- as.character(long.data$sample)
        if (return_csv) {
          return(long.data)
        }

        ## Apply sample ordering
        samples <- colnames(res$prop.counts)
        bars_order <- input$bars_order
        if (!is.null(bars_order)) {
          if (bars_order == "ascending") {
            top_vals <- head(res$prop.counts, 1)[1, samples]
            samples <- names(sort(top_vals))
          } else if (bars_order == "descending") {
            top_vals <- head(res$prop.counts, 1)[1, samples]
            samples <- names(sort(top_vals, decreasing = TRUE))
          } else if (bars_order == "custom" && !is.null(input$rank_list_basic) &&
            all(input$rank_list_basic %in% colnames(res$prop.counts))) {
            samples <- input$rank_list_basic
          }
        }
        long.data$sample <- factor(long.data$sample, levels = samples)

        ## Resolve palette
        n_genes <- length(unique(long.data$gene))
        palette <- if (!is.null(input$palette)) input$palette else "default"
        if (palette %in% c("default", "original")) {
          colors_vec <- omics_pal_d(palette = "expanded")(n_genes)
        } else if (palette == "custom") {
          colors_vec <- sapply(seq_len(n_genes), function(j) {
            val <- input[[paste0("custom_color_", j)]]
            if (is.null(val)) omics_pal_d(palette = "expanded")(8)[(j - 1) %% 8 + 1] else val
          })
        } else {
          colors_vec <- omics_pal_d(palette = palette)(n_genes)
        }

        ## stacked barchart
        fig <-
          plotly::plot_ly(
            data = long.data,
            x = ~sample,
            y = ~value,
            type = "bar",
            color = ~gene,
            colors = colors_vec,
            hovertemplate = ~ paste0(
              "Sample: <b>", sample, "</b><br>",
              "Gene: <b>", gene, "</b><br>",
              "Cum. proportion: <b>", sprintf("%2.1f", value), "%</b>",
              "<extra></extra>"
            )
          ) %>%
          plotly::layout(
            barmode = "stack",
            xaxis = list(title = FALSE),
            yaxis = list(title = "Cumulative proportion", ticksuffix = "%"),
            font = list(family = "Lato"),
            margin = list(l = 10, r = 10, b = 10, t = 10)
          ) %>%
          plotly_default()
        fig
      } else {
        avg.prop <- head(rowMeans(res$prop.counts, na.rm = TRUE), 15)
        genes <- head(res$gset.genes, 15)
        family <- paste0(names(avg.prop), "  ")
        family <- factor(family, levels = family)
        df <- data.frame(family = family, prop = avg.prop, genes = genes)
        if (return_csv) {
          return(df)
        }

        ## stacked barchart
        fig <-
          plotly::plot_ly(
            data = df,
            x = ~prop,
            y = ~family,
            type = "bar",
            marker = list(
              color = omics_colors("brand_blue")
            ),
            hovertemplate = ~ paste0(
              "Gene family: <b>", family, "</b><br>",
              "Genes: <b>", genes, "</b>",
              "<extra></extra>"
            )
          ) %>%
          plotly::layout(
            yaxis = list(title = FALSE),
            xaxis = list(title = "Proportion", ticksuffix = "%"),
            font = list(family = "Lato"),
            margin = list(l = 10, r = 10, b = 10, t = 10),
            showlegend = FALSE
          ) %>%
          plotly_default()
        fig
      }
    }

    modal_plotly.RENDER <- function() {
      fig <- plotly.RENDER() %>%
        plotly_modal_default() %>%
        plotly::layout(showlegend = FALSE)
      fig
    }

    plot_data_csv <- function() {
      df <- plotly.RENDER(return_csv = TRUE)
      return(df)
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data_csv,
      res = c(90, 170),
      pdf.width = 6,
      pdf.height = 6,
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
