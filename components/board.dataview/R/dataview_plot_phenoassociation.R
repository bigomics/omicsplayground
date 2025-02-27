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
    info.methods,
    info.extra_link,
    caption) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    withTooltip(
      shiny::selectInput(ns("vars"), "show variables:", choices = NULL, multiple = TRUE),
      "Select phenotype variables to show.",
      placement = "top"
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    caption = caption,
    options = opts,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

dataview_plot_phenoassociation_server <- function(id, pgx, r.samples, watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    observeEvent(pgx$samples, {
      vars <- colnames(pgx$samples)
      vars <- unique(c(grep("^[.]|cluster", vars, invert = TRUE, value = TRUE), vars))
      sel.vars <- head(vars, 5)
      shiny::updateSelectInput(session, "vars", choices = vars, selected = sel.vars)
    })

    plot_data <- shiny::reactive({
      shiny::req(pgx$samples, input$vars)
      vars <- input$vars
      samples <- r.samples()
      annot <- pgx$samples
      if (!all(vars %in% colnames(annot))) {
        return(NULL)
      }
      if (!all(samples %in% rownames(annot))) {
        return(NULL)
      }
      annot <- annot[samples, vars, drop = FALSE]
      list(annot = annot)
    })

    plot.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)

      check_diversity_in_colums <- function(df) {
        sum(unlist(apply(df, 2, function(x) length(unique(x)) > 1))) > 1
      }

      if (check_diversity_in_colums(res$annot) && is.data.frame(res$annot)) {
        ## NOTE: the package doesnt allow to change the typeface, the spacing of the legend, sizes + formatting of labels, ...
        ## TODO: reimplement in plotly (not me as code is complex and not intuitive at all)
        ## TODO: use na.omit to prevent plot error. This removes all rows where any value is NA. shall we use imputation?
        clean_annot <- na.omit(res$annot)
        pq <- playbase::pgx.testPhenoCorrelation(
          df = clean_annot,
          plot = TRUE,
          cex = 0.8
        )
        return(NULL)
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
      res = c(100, 170) * 0.85, ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
