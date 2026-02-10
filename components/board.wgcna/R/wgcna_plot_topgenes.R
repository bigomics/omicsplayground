##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_topgenes_ui <- function(
  id,
  label = "",
  title = "",
  info.text = "",
  caption = "",
  height,
  width
) {
  ns <- shiny::NS(id)


  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    # options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_plot_topgenes_server <- function(id,
                                       enrichTable,
                                       watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    get_data <- function() {
      df <- enrichTable$data()
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      ii <- enrichTable$rows_all()
      shiny::req(ii)
      df <- df[ii, , drop = FALSE]
      df <- head(df, 20)
      genes <- unlist(strsplit(df$genes, split = "\\|"))
      gtable <- sort(table(genes), decreasing = TRUE)
      gtable
    }

    plot.RENDER <- function() {
      gtable <- get_data()

      par(mar = c(8, 4, 2, 0))
      barplot(
        head(gtable, 35),
        las = 3,
        cex.names = 0.9,
        xlab = "",
        ylab = "Frequency in top"
      )
    }


    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      csvFunc = get_data,
      pdf.width = 10, pdf.height = 5,
      res = c(80, 110),
      add.watermark = watermark
    )
  })
}
