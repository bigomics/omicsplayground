##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_moduleheatmap_ui <- function(
  id,
  title = "",
  info.text = "",
  caption = "",
  label = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("split"),
      label = "Split heatmap by data types",
      value = TRUE
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

mofa_plot_moduleheatmap_server <- function(id,
                                           mofa,
                                           pgx,
                                           input_factor = reactive(1),
                                           show_types = reactive(NULL),
                                           ntop = c(40, 40),
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    render_plot <- function(legend, maxchar, n, show_colnames) {
      res <- mofa()
      k <- input_factor()
      factors <- colnames(res$F)
      if (!is.null(k)) shiny::req(k %in% factors)

      show_types <- show_types()
      dtypes <- names(res$ww)
      if (is.null(show_types)) show_types <- dtypes
      shiny::validate(need(
        length(show_types) > 0,
        "Please select at least one datatype"
      ))

      playbase::mofa.plot_heatmap(
        res,
        k = k, ## main=k,
        gene_table = pgx$genes,
        ntop = n,
        split = input$split,
        type = "splitmap",
        annot = "pheno",
        maxchar = maxchar,
        show_types = show_types,
        show_legend = legend,
        show_colnames = show_colnames,
        mar = c(3, 0, 0, 0),
        annot.ht = 0.9,
        cexRow = 0.9
      )
    }

    plot.RENDER <- function() {
      render_plot(
        legend = FALSE, maxchar = 30, n = ntop[1],
        show_colnames = FALSE
      )
    }

    plot.RENDER2 <- function() {
      render_plot(
        legend = TRUE, maxchar = 80, n = ntop[2],
        show_colnames = TRUE
      )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,
      pdf.width = 8,
      pdf.height = 12,
      res = c(80, 100),
      add.watermark = watermark
    )
  })
}
