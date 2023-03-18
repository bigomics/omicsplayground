##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

dataview_plot_tsne_ui <- function(id, label = "", height = c(350, 600)) {
  ns <- shiny::NS(id)

  info_text <- paste0("<b>T-SNE clustering</b> of samples (or cells) colored by an expression of the gene selected in the <code>search_gene</code> dropdown menu. The red color represents an over-expression of the selected gene across samples (or cells).")

  PlotModuleUI(
    ns("pltmod"),
    plotlib = "plotly",
    info.text = info_text,
    download.fmt = c("png", "pdf", "csv"),
    width = c("auto", "100%"),
    height = height,
    label = label,
    title = "t-SNE clustering"
  )
}

dataview_plot_tsne_server <- function(id,
                                      pgx,
                                      r.gene = reactive(""),
                                      r.samples = reactive(""),
                                      r.data_type = reactive("counts"),
                                      r.groupby = reactive(""),
                                      watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_dl <- reactiveValues()

    plot_data <- shiny::reactive({
      shiny::req(pgx$X, pgx$Y, pgx$genes, pgx$counts, pgx$samples, pgx$tsne2d)

      ##            gene <- parent.input$search_gene
      ##            samplefilter <- parent.input$data_samplefilter
      ##            data_type <- parent.input$data_type
      ##            groupby <- parent.input$groupby

      ## dereference reactives
      gene <- r.gene()
      samples <- r.samples()
      data_type <- r.data_type()
      groupby <- r.groupby()
      shiny::req(gene, data_type)

      if (samples[1] == "") samples <- colnames(pgx$X)

      ## precompute
      pp <- rownames(pgx$genes)[1]
      sel <- match(gene, pgx$genes$gene_name)
      pp <- rownames(pgx$genes)[ifelse(is.na(sel), 1, sel)]

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

      pos <- pgx$tsne2d[samples, ]

      fc1 <- tanh(0.99 * scale(gx)[, 1])
      fc1 <- tanh(0.99 * scale(gx, center = FALSE)[, 1])
      ## fc1 <- tanh(0.99 * gx/sd(gx))
      fc2 <- (fc1 - min(fc1))

      data <- data.frame(
        pos_x = pos[, 1],
        pos_y = pos[, 2],
        expression = gx,
        fc2 = fc2,
        name = rownames(pos)
      )

      grp <- NULL
      filt.groupby <- groupby
      if (!is.null(filt.groupby) && filt.groupby %in% colnames(pgx$samples)) {
        grp <- factor(pgx$samples[samples, filt.groupby])
        data$group <- grp
      }

      return(list(data = data, gene = gene))
    })

    plot.RENDER <- function() {
      data <- plot_data()[[1]]
      shiny::req(data)

      fig_base <-
        ggplot(data, aes(pos_x, pos_y)) +
        labs(x = "tSNE1", y = "tSNE2") +
        scale_color_viridis_c(
          option = "rocket",
          direction = -1,
          begin = .05, end = .97,
          limits = c(0, 1),
          labels = function(x) sprintf("%1.2f", x),
          name = "Expression"
        ) +
        guide_continuous(aes = "color", type = "steps", width = .4) +
        theme_omics(base_size = 12, axis_num = "xy", legendnum = TRUE)

      plot_dl$base <- fig_base

      if (!is.null(plot_data()$group)) {
        fig <- fig_base +
          ggforce::geom_mark_hull(
            aes(
              fill = stage(group, after_scale = colorspace::desaturate(fill, 1)),
              label = group
            ),
            color = "grey33",
            size = .4,
            alpha = .33 / length(unique(plot_data()$group)),
            expand = unit(2.7, "mm"),
            con.cap = unit(.01, "mm"),
            con.colour = "grey33",
            label.buffer = unit(2, "mm"),
            label.fontsize = 12.5,
            label.fontface = "plain"
          ) +
          geom_point(
            aes(
              color = stage(fc2, after_scale = colorspace::darken(color, .35)),
              fill = after_scale(color)
            ),
            size = 1.8,
            shape = 21,
            stroke = .5
          ) +
          scale_x_continuous(expand = c(.4, .4)) +
          scale_y_continuous(expand = c(.4, .4)) +
          scale_fill_discrete(guide = "none")
      } else {
        fig <- fig_base +
          geom_point(
            aes(
              color = stage(fc2, after_scale = colorspace::darken(color, .35)),
              fill = after_scale(color)
            ),
            size = 2.3,
            shape = 21,
            stroke = .5
          )
      }

      plot_dl$plot <- fig
      ## fig <- plotly::ggplotly(fig)
      fig
    }

    modal_plot.RENDER <- function() {
      fig <- plot.RENDER() +
        guide_continuous(aes = "color", type = "steps", width = .7) +
        theme_omics(base_size = 20, axis_num = "xy", legendnum = TRUE)
      ## plotly::ggplotly(fig)
      fig
    }

    plotly.RENDER0 <- function() {
      data <- plot_data()
      shiny::req(data)

      df <- data[[1]]
      gene <- data[[2]]
      symbols <- c("circle", "square", "cross", "diamond", "triangle-down", "star", "x", "trianlge-up",
        "star-diamond", "square-cross", "diamond-wide")

      if (!is.null(df$group)) {
        fig <-
          plotly::plot_ly(
            data = df,
            type = "scatter",
            mode = "markers",
            x = ~pos_x,
            y = ~pos_y,
            symbol = ~group,
            symbols = symbols[1:length(unique(df$group))],
            color = ~expression,
            colors = omics_pal_c(palette = "brand_blue")(100),
            marker = list(
              size = 10,
              line = list(
                color = omics_colors("super_dark_grey"),
                width = 1.2
              )
            ),
            hovertemplate = ~ paste(
              "Gene:<b>", gene,
              "</b><br>Sample:<b>", name, "</b><br>",
              "Expression:<b>", sprintf("%1.3f", expression), "</b><br>",
              "Group:<b>", sprintf("%s", as.character(df$group)), "</b>",
              "<extra></extra>"
            )
          )
      } else {
        fig <-
          plotly::plot_ly(
            data = df,
            type = "scatter",
            mode = "markers",
            x = ~pos_x,
            y = ~pos_y,
            color = ~expression,
            colors = omics_pal_c(palette = "brand_blue")(100),
            marker = list(
              size = 10,
              line = list(
                color = omics_colors("super_dark_grey"),
                width = 1.2
              )
            ),
            hovertemplate = ~ paste(
              "Gene:<b>", gene,
              "</b><br>Sample:<b>", name,
              "</b><br>Expression:<b>", sprintf("%1.3f", expression),
              "</b><extra></extra>"
            )
          )
      }
      fig %>%
        plotly::layout(
          xaxis = list(title = "tSNE-x"),
          yaxis = list(title = "tSNE-y"),
          margin = list(l = 10, r = 10, b = 10, t = 10)
        ) %>%
        plotly::colorbar(
          title = "<b>Expression:</b>",
          width = .001,
          ticklen = 6
        ) %>%
        plotly_default1() ## %>% toWebGL()
    }

    plotly.RENDER <- function() {
      fig <- plotly.RENDER0() # %>%
      #    plotly::hide_colorbar()
      fig
    }

    modal_plotly.RENDER <- function() {
      fig <- plotly.RENDER0() %>%
        plotly_modal_default() %>%
        plotly::style(
          marker.size = 20
        )
      fig
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      ## func = plot.RENDER,
      ## func2 = modal_plot.RENDER,
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(100, 300) * 1, ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      ## label = label, title = "t-SNE clustering",
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
