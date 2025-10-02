##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_pathwayheatmap_ui <- function(
    id,
    title = "",
    info.text = "",
    info.methods = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
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
    info.methods = info.methods,
    options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

mofa_plot_pathwayheatmap_server <- function(id,
                                            mofa,
                                            pgx,
                                            input_factor = reactive(1),
                                            selected = reactive(1),
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function(show_legend = FALSE) {
      mofa <- mofa()
      k <- input_factor()
      factors <- colnames(mofa$F)
      if (!is.null(k)) shiny::req(k %in% factors)

      pw <- selected()
      shiny::validate(need(length(pw) == 1, "Please select a pathway"))

      features <- NULL
      if (pw %in% colnames(pgx$GMT)) {
        genes <- names(which(pgx$GMT[, pw] != 0))
        genes <- playbase::mofa.strip_prefix(genes)
        features <- playbase::map_probes(pgx$genes, genes)
      }

      playbase::mofa.plot_heatmap(
        mofa,
        gene_table = pgx$genes,
        k = k, ## main=k,
        features = features,
        ntop = 40,
        # split = input$split,
        split = FALSE,
        type = "splitmap",
        annot = "pheno",
        maxchar = 40,
        show_legend = show_legend,
        show_types = NULL,
        mar = c(3, 0, 0, 0),
        annot.ht = 0.9,
        cexRow = 0.9
      )
    }

    plot.RENDER2 <- function() {
      plot.RENDER(show_legend = TRUE)
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,
      pdf.width = 8, pdf.height = 12,
      res = c(80, 100),
      add.watermark = watermark
    )
  })
}
