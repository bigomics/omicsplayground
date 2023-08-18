##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

upload_plot_contraststats_ui <- function(id,
                                         label = "",
                                         height,
                                         width,
                                         title,
                                         caption,
                                         info.text) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    plotlib = "base",
    info.text = info.text,
    caption = caption,
    options = NULL,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
  )
}

upload_plot_contraststats_server <- function(id, checkTables, uploaded, watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## extract data from pgx object
    plot_data <- shiny::reactive({
      ct <- uploaded$contrasts.csv
      has.contrasts <- !is.null(ct) && NCOL(ct) > 0
      check <- checkTables()

      status.ok <- check["contrasts.csv", "status"]
      status.ds <- tolower(check["contrasts.csv", "description"])
      error.msg <- paste(
        toupper(status.ok), "\nPlease upload 'contrasts.csv' (Optional):",
        status.ds
      )

      shiny::validate(
        shiny::need(
          status.ok == "OK" && has.contrasts,
          error.msg
        )
      )

      contrasts <- uploaded$contrasts.csv
      return(contrasts)
    })

    plot.RENDER <- function() {
      df <- plot_data()
      px <- head(colnames(df), 20) ## maximum to show??
      df <- data.frame(df[, px, drop = FALSE], check.names = FALSE)
      tt2 <- paste(nrow(df), "samples x", ncol(df), "comparisons")

      p1 <- df %>%
        inspectdf::inspect_cat() %>%
        inspectdf::show_plot()

      p1 <- p1 + ggplot2::ggtitle("COMPARISONS", subtitle = tt2) +
        ggplot2::theme(
          #
          axis.text.y = ggplot2::element_text(
            size = 12,
            margin = ggplot2::margin(0, 0, 0, 25),
            hjust = 1
          )
        )

      return(p1)
    }

    modal_plot.RENDER <- function() {
      plot.RENDER()
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "base",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 90), ## resolution of plots
      pdf.width = 4, pdf.height = 4,
      filename = "upload_plot_contraststats",
      add.watermark = watermark
    )

  }) ## end of moduleServer
}
