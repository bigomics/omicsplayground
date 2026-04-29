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
      choices = "<ungrouped>"
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
    plot_type = "expression_boxplot"
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

    output$custom_palette_ui <- shiny::renderUI({
      shiny::req(input$palette == "custom")
      res <- plot_data()
      shiny::req(res)
      annot <- res$annot
      shiny::req(annot)
      class_col <- input$group_by_feature_class
      shiny::req(class_col, class_col != "<ungrouped>", class_col %in% colnames(annot))
      classes <- sort(unique(stats::na.omit(as.character(annot[, class_col]))))
      shiny::req(length(classes) > 0)
      theme_palette <- shiny::isolate(get_color_theme()$palette)
      if (is.null(theme_palette) || theme_palette == "") theme_palette <- "default"
      default_clrs <- tryCatch(
        omics_pal_d(palette = theme_palette)(length(classes)),
        error = function(e) omics_pal_d("default")(length(classes))
      )
      custom_palette_pickers(classes, session$ns, default_clrs)
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
        df <- df[sample(nrow(df), 1000), , drop = FALSE]
      }
      long.df <- reshape2::melt(df)
      colnames(long.df) <- c("gene", "sample", "value")
      long.df$sample <- as.character(long.df$sample)

      split <- NULL
      annot <- res$annot
      if (!is.null(annot)) {
        class_col <- input$group_by_feature_class
        if (class_col != "<ungrouped>" && class_col %in% colnames(annot)) {
          idx <- which(long.df$gene %in% rownames(annot))
          if (length(idx) > 0) {
            long.df <- long.df[idx, , drop = FALSE]
            long.df$class <- annot[match(long.df$gene, rownames(annot)), class_col]
            split <- "class"
          }
        }
      }

      palette_colors <- NULL
      if (!is.null(split)) {
        classes <- sort(unique(stats::na.omit(as.character(long.df$class))))
        n_classes <- length(classes)
        if (n_classes > 0) {
          theme_palette <- shiny::isolate(get_color_theme()$palette)
          if (is.null(theme_palette) || theme_palette == "") theme_palette <- "default"
          fallback_pal <- tryCatch(
            omics_pal_d(palette = theme_palette)(n_classes),
            error = function(e) omics_pal_d("default")(n_classes)
          )
          palette_colors <- resolve_palette_colors(
            input, n_classes,
            fallback_colors = fallback_pal
          )
          if (!is.null(palette_colors)) {
            palette_colors <- stats::setNames(palette_colors, classes)
          }
        }
      }

      bar_color <- if (!is.null(split)) "black" else get_editor_color(input, "scatter_color", "secondary")
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

      gp <- extract_ggprism_params(input)

      if (gp$use_ggprism) {
        ## --- ggplot2 + ggprism path ---
        p <- playbase::pgx.boxplot.GGPLOT(
          data = long.df, x = "sample", y = "value",
          split = split,
          yaxistitle = ylab,
          color = bar_color,
          fillcolor = fill_color,
          linecolor = bar_color,
          colors = palette_colors
        )
        p <- apply_ggprism_theme(p, gp, x_angle = 90)
        p <- apply_editor_theme(p, input)
        fig <- ggplot_as_plotly_image(p)
      } else {
        ## --- existing plotly path ---
        fig <- playbase::pgx.boxplot.PLOTLY(
          data = long.df,
          x = "sample",
          y = "value",
          split = split,
          yaxistitle = ylab,
          color = bar_color,
          fillcolor = fill_color,
          linecolor = bar_color,
          colors = palette_colors
        ) %>%
          plotly_default()
        fig <- apply_prism_plotly(fig, gp)
      }

      if (!gp$use_ggprism) fig <- apply_plotly_editor_theme(fig, input)
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
