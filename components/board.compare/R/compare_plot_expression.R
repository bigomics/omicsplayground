##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

compare_plot_expression_ui <- function(
    id,
    label = "",
    height = c(600, 800),
    title,
    info.text) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    plotlib = "plotly",
    label = "a",
    info.text = info.text,
    height = height,
    width = c("auto", "100%"),
    download.fmt = c("png", "pdf")
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
