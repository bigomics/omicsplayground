##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

TimeSeriesBoard.parcoord_plot_ui <- function(
    id,
    label = "label",
    title = "title",
    info.text = "info.text",
    info.methods = "info.methods",
    info.references = list(),
    info.extra_link = "extra.link",
    caption = "caption",
    height = c("calc(100vh - 310px)", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")
    ) {
  ns <- shiny::NS(id)

  parcoord_opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(ns("average"), tspan("Average by gene module"), FALSE),
      "Average gene by gene module"
    ),
    withTooltip(
      shiny::checkboxInput(ns("scale"), "Scale values", TRUE),
      "Scale expression values to mean=0 and SD=1."
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    options = parcoord_opts,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

TimeSeriesBoard.parcoord_table_ui <- function(
    id,
    label = "label",
    title = "title",
    info.text = "info.text",
    caption = "caption",
    height = c("40%", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")
    ) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("table"),
    info.text = info.text,
    height = height,
    caption = caption,
    width = width,
    title = title,
    label = "b"
  )
}

TimeSeriesBoard.parcoord_server <- function(id,
                                            pgx,
                                            data,
                                            select_module,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- function() {
      res <- data()
      selmod <- select_module()
      shiny::req(selmod)

      dbg("[TimeSeriesBoard.parcoord_server:plot_data] selmod=", selmod)
      
      ii <- which(res$colors %in% selmod)
      timeX <- res$X[ii,,drop=FALSE]
      colors <- res$colors[ii]

      res <- list(
        timeX = timeX,
        colors = colors
      )
    }

    plot.RENDER <- function() {
      res <- plot_data()
      timeX <- res$timeX
      
      dbg("[TimeSeriesBoard.parcoord_server:plot.RENDER] dim.timeX=", dim(timeX))
      dbg("[TimeSeriesBoard.parcoord_server:plot.RENDER] head.colors=", head(res$colors))      
      dimensions <- list()
      for(i in 1:ncol(timeX)) {
        d <- list(
          range = range(timeX),
          label = colnames(timeX)[i],
          values = timeX[,i]
        )
        dimensions[[i]] <- d
      }
      
      df <- data.frame(timeX, check.names=FALSE)
      int.colors <- as.integer(factor(res$colors))

      plt <- plotly::plot_ly(
        df, 
        source = "pcoords"
      ) %>%
        plotly::add_trace(
          type = "parcoords",
          line = list(
            color = int.colors,
            colorscale = "Jet",
            showscale = FALSE,
            width = 40
          ),
          dimensions = dimensions
        )
      plt

    }

    plot.RENDER_MODAL <- function() {
      plot.RENDER() %>%
        plotly_modal_default()
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plot.RENDER,
      func2 = plot.RENDER_MODAL,
      csvFunc = plot_data,
      res = c(90, 170),
      pdf.width = 8,
      pdf.height = 8,
      add.watermark = watermark
    )

    ##-----------------------------------------------------
    ##------------------- gene table ----------------------
    ##-----------------------------------------------------
    table.RENDER <- function() {

      res <- plot_data()
      timeX <- res$timeX

      dbg("[TimeSeriesBoard.parcoord_server:table.RENDER] 0:")
      
      df <- data.frame(
        module = res$colors,
        feature = rownames(timeX),
        timeX,
        check.names = FALSE
      )
      dbg("[TimeSeriesBoard.parcoord_server:table.RENDER] 1:")

      symbol <- pgx$genes[rownames(timeX),"symbol"]
      if(mean(symbol == rownames(timeX), na.rm=TRUE) < 0.2) {
        df$symbol <- symbol
      }

      dbg("[TimeSeriesBoard.parcoord_server:table.RENDER] 2:")
      
      cols <- c("module","feature","symbol",colnames(df))
      cols <- intersect(cols, colnames(df))
      df <- df[,cols]
      
      numeric.cols <- colnames(timeX)
      DT::datatable(
        df,
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = "23vh",
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    table.RENDER_modal <- function() {
      dt <- parcoord_table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    }

    TableModuleServer(
      "table",
      func = table.RENDER,
      func2 = table.RENDER_modal,
      selector = "none"
    )
  })
}
