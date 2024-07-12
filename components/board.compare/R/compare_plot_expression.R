##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

compare_plot_expression_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)

  info_text <- "<b>Multi barplots.</b> Barplots of expression values for multiple comparisons in the two datasets (blue and green). "

  PlotModuleUI(
    ns("plot"),
    title = "Expression",
    plotlib = "plotly",
    label = "a",
    info.text = info_text,
    height = height,
    width = c("auto", "100%"),
    download.fmt = c("png", "pdf")
  )
}

compare_plot_expression_server <- function(id,
                                           pgx,
                                           dataset2,
                                           input.contrast1,
                                           input.contrast2,
                                           hilightgenes,
                                           getOmicsScoreTable,
                                           score_table,
                                           compute,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    contrast1 <- shiny::reactiveVal(FALSE)
    contrast2 <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(compute(), {
      contrast1(input.contrast1())
      contrast2(input.contrast2())
    })

    plotly_multibarplot.RENDER <- shiny::reactive({
      dt <- tryCatch({
        getOmicsScoreTable()
      }, error = function(w) {
        FALSE
      })
      shiny::validate(shiny::need(dt, "Please select contrasts and run 'Compute'"))
      shiny::req(contrast1())
      shiny::req(contrast2())
      shiny::req(getOmicsScoreTable())
      shiny::req(hilightgenes())
      shiny::req(score_table())
      pgx1 <- pgx
      pgx2 <- dataset2()

      org1 <- playbase::pgx.getOrganism(pgx1)
      org2 <- playbase::pgx.getOrganism(pgx2)

      ct1 <- head(names(pgx1$gx.meta$meta), 3)
      ct2 <- head(names(pgx2$gx.meta$meta), 3)
      ct1 <- contrast1()
      ct2 <- contrast2()
      if (is.null(pgx1$version) && is.null(pgx2$version)) {
        target_col1 <- target_col2 <- "gene_name"
      } else if (org1 == org2) {
        target_col1 <- target_col2 <- "symbol"
      } else if (org1 != org2) {
        target_col1 <- target_col2 <- "human_ortholog"
        if (!target_col1 %in% colnames(pgx1$genes)) target_col1 <- "gene_name"
        if (!target_col2 %in% colnames(pgx2$genes)) target_col2 <- "gene_name"
      }

      X1 <- playbase::rename_by(pgx1$X, pgx1$genes, target_col1)
      X2 <- playbase::rename_by(pgx2$X, pgx2$genes, target_col2)

      genes <- rownames(X1)
      genes <- hilightgenes()

      df <- getOmicsScoreTable()

      sel <- score_table() ## from module
      genes <- rownames(df)[sel]

      xgenes <- intersect(rownames(X1), rownames(X2))
      genes <- head(intersect(genes, xgenes), 8)
      e1 <- playbase::contrastAsLabels(pgx1$model.parameters$exp.matrix[, ct1, drop = FALSE])
      e2 <- playbase::contrastAsLabels(pgx2$model.parameters$exp.matrix[, ct2, drop = FALSE])

      # Build plots
      sub_plots <- vector("list", length(genes))
      names(sub_plots) <- genes
      for (gene_i in genes) {
        # Get data
        x1 <- pgx1$X[gene_i, ]
        x2 <- pgx2$X[gene_i, ]
        m1 <- lapply(e1, function(y) tapply(x1, y, mean))
        m2 <- lapply(e2, function(y) tapply(x2, y, mean))
        mm <- cbind(do.call(cbind, m1), do.call(cbind, m2))

        # Assemble barplots
        title_y <- max(mm) + max(mm) * .1
        rn <- rownames(mm)
        # Add bars
        plt <- plotly::plot_ly()
        for (i in seq_len(ncol(mm))) {
          col_i <- mm[, i, drop = FALSE]
          name_i <- gsub(pattern = "_vs_", replacement = " vs\n", x = colnames(col_i))
          show_legend1 <- (i == 1 && gene_i == genes[1])
          plt <- plotly::add_trace(plt,
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
          text = paste("<b>", gene_i, "</b>"),
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

      return(all_plt)
    })


    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = plotly_multibarplot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(95, 130),
      add.watermark = watermark
    )
  })
}
