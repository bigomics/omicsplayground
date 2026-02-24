##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

compare_plot_expression_ui <- function(
  id,
  label = "",
  height = c(600, 800),
  title,
  info.text
) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    plotlib = "plotly",
    label = "a",
    info.text = info.text,
    height = height,
    width = c("auto", "100%"),
    download.fmt = c("png", "pdf", "svg"),
    editor = TRUE,
    ns_parent = ns,
    plot_type = "grouped_barplot",
    palette_default = "original"
  )
}

compare_plot_expression_server <- function(id,
                                           pgx,
                                           dataset2,
                                           contrast1,
                                           contrast2,
                                           hilightgenes,
                                           getMatrices,
                                           getScoreTable,
                                           selected,
                                           ## compute,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Override bars_order to "custom" (only drag-and-drop makes sense for contrasts)
    shiny::observeEvent(TRUE, once = TRUE, {
      shiny::updateSelectInput(
        session, "bars_order",
        selected = "custom"
      )
    })

    ## Shared reactive: compute group names and contrast names for editor UI
    gene_data <- shiny::reactive({
      dt <- tryCatch(getScoreTable(), error = function(w) FALSE)
      shiny::req(dt)
      shiny::req(contrast1(), contrast2(), hilightgenes(), selected())

      pgx1 <- pgx
      pgx2 <- dataset2()
      mat <- getMatrices()
      X1 <- mat$X1
      X2 <- mat$X2
      xgenes <- intersect(rownames(X1), rownames(X2))
      sel.genes <- head(intersect(selected(), xgenes), 8)

      e1 <- playbase::pgx.getContrastMatrix(pgx1)[, contrast1(), drop = FALSE]
      e2 <- playbase::pgx.getContrastMatrix(pgx2)[, contrast2(), drop = FALSE]
      ## get group names and contrast names from the first gene
      rn <- NULL
      contrast_names <- NULL
      if (length(sel.genes) > 0) {
        x1 <- X1[sel.genes[1], ]
        x2 <- X2[sel.genes[1], ]
        m1 <- apply(e1, 2, function(y) tapply(x1, y, mean))
        m2 <- apply(e2, 2, function(y) tapply(x2, y, mean))
        mm <- cbind(m1, m2)
        rn <- rownames(mm)
        contrast_names <- colnames(mm)
      }

      list(sel.genes = sel.genes, rn = rn, contrast_names = contrast_names)
    })

    ## Editor: dynamic color pickers for custom palette
    output$custom_palette_ui <- shiny::renderUI({
      shiny::req(input$palette == "custom")
      gd <- gene_data()
      shiny::req(gd$rn)
      rn <- gd$rn
      default_pal <- omics_pal_d()(length(rn))
      pickers <- lapply(seq_along(rn), function(i) {
        colourpicker::colourInput(
          ns(paste0("custom_color_", i)),
          label = rn[i],
          value = default_pal[i]
        )
      })
      shiny::tagList(pickers)
    })

    ## Editor: rank list for custom drag-and-drop contrast ordering
    output$rank_list <- shiny::renderUI({
      gd <- gene_data()
      shiny::req(gd$contrast_names)
      labels <- gsub("_vs_", " vs ", gd$contrast_names)
      names(labels) <- gd$contrast_names
      sortable::bucket_list(
        header = NULL,
        class = "default-sortable custom-sortable",
        sortable::add_rank_list(
          input_id = ns("rank_list_order"),
          text = NULL,
          labels = labels
        )
      )
    })

    plotly_multibarplot.RENDER <- shiny::reactive({
      dt <- tryCatch(
        {
          getScoreTable()
        },
        error = function(w) {
          FALSE
        }
      )
      shiny::validate(shiny::need(dt, "Please select contrasts and run 'Compute'"))
      shiny::req(contrast1())
      shiny::req(contrast2())
      shiny::req(getScoreTable())
      shiny::req(hilightgenes())
      shiny::req(selected())
      pgx1 <- pgx
      pgx2 <- dataset2()

      org1 <- playbase::pgx.getOrganism(pgx1)
      org2 <- playbase::pgx.getOrganism(pgx2)

      ct1 <- head(names(pgx1$gx.meta$meta), 3)
      ct2 <- head(names(pgx2$gx.meta$meta), 3)
      ct1 <- contrast1()
      ct2 <- contrast2()

      mat <- getMatrices()
      X1 <- mat$X1
      X2 <- mat$X2

      df <- getScoreTable()
      sel.genes <- selected()
      xgenes <- intersect(rownames(X1), rownames(X2))
      sel.genes <- head(intersect(sel.genes, xgenes), 8)

      e1 <- playbase::pgx.getContrastMatrix(pgx1)[, ct1, drop = FALSE]
      e2 <- playbase::pgx.getContrastMatrix(pgx2)[, ct2, drop = FALSE]

      ## Editor: resolve group colors from first gene
      palette <- input$palette

      ## Compute rn (group names) from the first gene
      x1_first <- X1[sel.genes[1], ]
      m1_first <- apply(e1, 2, function(y) tapply(x1_first, y, mean))
      rn <- rownames(m1_first)

      ## Resolve group colors based on palette selection
      group_colors <- NULL ## NULL means use plotly defaults (original)
      if (!is.null(palette) && palette == "custom") {
        default_pal <- omics_pal_d()(length(rn))
        group_colors <- sapply(seq_along(rn), function(i) {
          val <- input[[paste0("custom_color_", i)]]
          if (is.null(val)) default_pal[i] else val
        })
        names(group_colors) <- rn
      } else if (!is.null(palette) && !palette %in% c("original", "")) {
        group_colors <- omics_pal_d(palette = palette)(length(rn))
        names(group_colors) <- rn
      }

      # Build plots
      sub_plots <- vector("list", length(sel.genes))
      names(sub_plots) <- sel.genes

      for (gene_i in sel.genes) {
        # Get data
        x1 <- X1[gene_i, ]
        x2 <- X2[gene_i, ]
        m1 <- apply(e1, 2, function(y) tapply(x1, y, mean))
        m2 <- apply(e2, 2, function(y) tapply(x2, y, mean))
        mm <- cbind(m1, m2)

        # Assemble barplots
        title_y <- 1.1 * max(mm, na.rm = TRUE)
        rn <- rownames(mm)
        plt <- plotly::plot_ly()
        if (!is.null(group_colors)) {
          plt <- plotly::plot_ly(colors = group_colors)
        }
        for (i in seq_len(NCOL(mm))) {
          col_i <- mm[, i, drop = FALSE]
          name_i <- gsub(pattern = "_vs_", replacement = " vs\n", x = colnames(col_i))
          show_legend1 <- (i == 1 && gene_i == sel.genes[1])

          # Add bars
          plt <- plotly::add_trace(
            plt,
            x = name_i,
            y = col_i[, 1] + 0.00001, # We add 0.00001 so that columns don't get ignored
            name = rn,
            color = rn,
            type = "bar",
            width = 0.25,
            showlegend = show_legend1
          ) %>%
            plotly::layout(
              xaxis = list(
                titlefont = list(size = 5),
                tickangle = 45
              ),
              legend = list(orientation = "h", bgcolor = "transparent", y = 1.2),
              barmode = "group"
            )
        }

        plt <- plotly::add_annotations(plt,
          text = paste("<b>", playbase::probe2symbol(gene_i, pgx$genes, "gene_name", fill_na = TRUE), "</b>"),
          font = list(size = 9),
          showarrow = FALSE,
          xanchor = "left",
          yanchor = "bottom",
          x = 0.5,
          y = title_y
        )
        sub_plots[[gene_i]] <- plt
      }

      # Put plots together
      suppressWarnings(
        all_plt <- plotly::subplot(sub_plots,
          nrows = 2, margin = 0.025,
          titleX = TRUE, titleY = TRUE, shareX = TRUE
        )
      ) %>%
        plotly::layout(margin = list(l = 0, b = 10, r = 0))

      ## Editor: apply contrast ordering via categoryorder on all subplot axes
      if (!is.null(input$rank_list_order)) {
        ordered_cats <- gsub("_vs_", " vs\n", input$rank_list_order)
        axis_layout <- list(categoryorder = "array", categoryarray = ordered_cats)
        layout_args <- list(all_plt)
        for (k in seq_len(length(sel.genes))) {
          ax_name <- if (k == 1) "xaxis" else paste0("xaxis", k)
          layout_args[[ax_name]] <- axis_layout
        }
        all_plt <- do.call(plotly::layout, layout_args)
      }

      return(all_plt)
    })


    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = plotly_multibarplot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(95, 130),
      add.watermark = watermark,
      parent_session = session
    )
  })
}
