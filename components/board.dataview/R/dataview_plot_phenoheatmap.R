##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

dataview_plot_phenoheatmap_ui <- function(
  id,
  label = "",
  height,
  width,
  title,
  info.text,
  info.methods,
  caption
) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("clustsamples"), "cluster samples", TRUE),
      "Cluster samples.",
      placement = "top"
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    info.text = info.text,
    info.methods = info.methods,
    caption = caption,
    options = opts,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

dataview_plot_phenoheatmap_server <- function(id, pgx, r.samples, watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## extract data from pgx object
    plot_data <- shiny::reactive({
      shiny::req(pgx$samples, r.samples())

      samples <- r.samples()
      do.clust <- input$clustsamples

      annot <- pgx$samples
      annot <- annot[samples, , drop = FALSE]

      list(
        annot = annot,
        do.clust = do.clust
      )
    })

    plot.RENDER <- function() {
      res <- plot_data()
      shiny::req(res$annot)

      annot.ht <- ifelse(ncol(res$annot) > 10, 5, 6)
      annot.ht <- ifelse(ncol(res$annot) > 20, 4, annot.ht)
      annot.ht <- ifelse(ncol(res$annot) > 30, 3, annot.ht)

      check_diversity_in_colums <- function(df) {
        sum(
          unlist(
            apply(
              df, 2, function(x) {
                length(unique(x)) > 1
              }
            )
          )
        ) > 1
      }

      if (check_diversity_in_colums(res$annot) && is.data.frame(res$annot)) {
        ## TODO: Color palettes should be unique, not the same for condition and time
        ## NOTE: the package doesnt allow to change the typeface, the position of the legend, the label placement, ...
        ## TODO: reimplement in plotly (not me as code is complex and not intuitive at all)
        plt <- playbase::pgx.plotPhenotypeMatrix0(
          annot = res$annot,
          annot.ht = annot.ht,
          cluster.samples = res$do.clust
        )
        #
        plt
      } else {
        shiny::validate(shiny::need(nrow(res) > 0, "The filters have no diference across samples,please choose another filter."))
        return(NULL)
      }
    }

    modal_plot.RENDER <- function() {
      plot.RENDER()
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "base",
      plotlib2 = "base",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot,
      res = c(60, 120), ## resolution of plots
      pdf.width = 8, pdf.height = 5,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
