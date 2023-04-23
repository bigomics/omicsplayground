##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

dataview_plot_phenoassociation_ui <- function(
  id,
  label = "",
  height,
  width,
  title,
  info.text,
  caption) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("phenoclustsamples"), "cluster samples", TRUE),
      "Cluster samples.",
      placement = "top"
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    options = opts,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
  )
}

dataview_plot_phenoassociation_server <- function(id, pgx, r.samples, watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X, pgx$Y)
      samples <- r.samples()
      annot <- pgx$samples
      annot <- annot[samples, , drop = FALSE]
      list(annot = annot)
    })

    plot.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)
      
      check_diversity_in_colums <- function(df){
        sum(
          unlist(
            apply(
              df, 2, function(x){
                length(unique(x))>1
                }
                )
              )
            ) >1
      }
      

      if (check_diversity_in_colums(res$annot) && is.data.frame(res$annot)) {
        ## NOTE: the package doesnt allow to change the typeface, the spacing of the legend, sizes + formatting of labels, ...
        ## TODO: reimplement in plotly (not me as code is complex and not intuitive at all)
        pq <- playbase::pgx.testPhenoCorrelation(
          df = res$annot, 
          plot = TRUE,
          cex = 1
          ) 
        return(pq)  
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
      csvFunc = plot_data,   ##  *** downloadable data as CSV
      renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot,
      res = c(100, 170) * 0.85, ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
