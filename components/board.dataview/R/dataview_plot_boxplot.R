## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.

dataview_plot_boxplot_ui <- function(
  id,
  label = "",
  height,
  title,
  caption,
  info.text
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::radioButtons(
      inputId = ns("group_by_feature_class"),
      label = "Group by feature class (available when {Group by} is used)",
      choices = "a"
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    plotlib = "plotly",
    label = label,
    caption = caption,
    info.text = info.text,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    options = options,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "expression_barplot"
  )
}

dataview_plot_boxplot_server <- function(id,
                                         parent.input,
                                         getCountsTable,
                                         r.samples = reactive(""),
                                         r.annot = reactive(""),
                                         r.data_type,
                                         watermark = FALSE) {

  moduleServer(id, function(input, output, session) {
    
    shiny::observe({
      annot <- r.annot()
      display <- if (is.null(annot)) "none" else ""
      shinyjs::runjs(sprintf(
        'var modal = document.getElementById("%s");
         var card = modal ? modal.closest(".card") : null;
         var btn = card ? card.querySelector("[class*=\'fa-bars\']") : null;
         if (btn) btn.closest(".btn").style.display = "%s";',
        session$ns("pltmod-plotPopup"), display
      ))
      if (!is.null(annot)) {
        shiny::updateRadioButtons(
          session,
          "group_by_feature_class",
          choices = c("<ungrouped>", colnames(annot)) 
        )
      }
    })
    
    plot_data <- shiny::reactive({
      res <- getCountsTable()
      samples <- r.samples()
      annot <- r.annot()
      shiny::req(res)
      list(counts = res$log2counts, sample = colnames(res$log2counts), annot = annot)
    })
    
    output$rank_list <- shiny::renderUI({
      res <- plot_data()
      shiny::req(res)
      rank_list_ui(as.character(res$sample), session$ns)
    })

    plotly.RENDER <- function() {

      res <- plot_data()
      shiny::req(res)
      data_type <- r.data_type()

      if (data_type == "counts") {
        ylab <- "Counts (log2CPM)"
      } else if (data_type == "abundance") {
        ylab <- "Abundance"
      } else {
        ylab <- "Abundance (log2)"
      }

      df <- res$counts[, , drop = FALSE]

      if (nrow(df) > 1000) {
        sel <- sample(nrow(df), 1000)
        df <- df[sel, , drop = FALSE]
      }
      long.df <- reshape2::melt(df)
      colnames(long.df) <- c("gene", "sample", "value")
      long.df$sample <- as.character(long.df$sample)

      annot <- res$annot
      groupedByClass <- FALSE
      if (!is.null(annot)) {
        kk <- input$group_by_feature_class
        if (kk != "<ungrouped>") {
          jj <- which(long.df$gene %in% rownames(annot))
          if (kk %in% colnames(annot) & length(jj) > 0) {
            long.df <- long.df[jj, , drop = FALSE]
            jj <- match(long.df$gene, rownames(annot))
            long.df$class <- annot[jj, kk]
            groupedByClass <- TRUE
          }
        }
      }
      
      bar_color <- get_editor_color(input, "scatter_color", "secondary")
      fill_color <- adjustcolor(bar_color, alpha.f = 0.35)
      bars_order <- input$bars_order
      samples <- res$sample

      ## Apply sample ordering
      if (!is.null(bars_order)) {
        if (bars_order == "ascending") {
          medians <- tapply(long.df$value, long.df$sample, median, na.rm = TRUE)
          samples <- names(sort(medians))
        } else if (bars_order == "descending") {
          medians <- tapply(long.df$value, long.df$sample, median, na.rm = TRUE)
          samples <- names(sort(medians, decreasing = TRUE))
        } else if (bars_order == "custom" && !is.null(input$rank_list_basic) &&
          all(input$rank_list_basic %in% res$sample)) {
          samples <- input$rank_list_basic
        }
      }
      long.df$sample <- factor(long.df$sample, levels = samples)

      split <- if (groupedByClass) "class" else NULL
      fig <- playbase::pgx.boxplot.PLOTLY(
        data = long.df,
        x = "sample",
        y = "value",
        split = split,
        yaxistitle = ylab,
        color = bar_color,
        fillcolor = fill_color,
        linecolor = bar_color
      ) %>%
        plotly_default()

      fig

    }

    modal_plotly.RENDER <- function() {
      plotly.RENDER() %>%
        plotly_modal_default()
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 170), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
