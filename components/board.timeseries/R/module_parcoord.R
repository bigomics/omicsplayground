##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

TimeSeriesBoard.parcoord_plot_ui <- function(
    id,
    label = "label",
    title = "title",
    info.text = "Parallel line plot displaying the average expression per time point of features mapped within the selected gene module. ",
    info.methods = "The normalized and log2-transformed expression data are scaled and centered. Per each feature, the average expression across samples is calculated per each time point. The plot displays the average expression per time point of features mapped within the selected gene module.",
    # info.references = list(),
    info.extra_link = "extra.link",
    caption = "caption",
    height = c("calc(100vh - 310px)", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")) {
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
    # info.references = info.references,
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
    info.text = "Table reporting the features mapped within the selected time series clustering module. Table includes the average feature expression (log2-scale) across samples per each variable, and standard deviation of these average values.",
    caption = "Table reporting the features mapped within the selected time series clustering module. Table includes the average feature expression (log2-scale) across samples per each variable, and standard deviation of these average values.",
    height = c("40%", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")) {
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

      ii <- which(res$modules %in% selmod)
      timeZ <- res$Z[ii, , drop = FALSE]
      timeX <- res$X[ii, , drop = FALSE]
      modules <- res$modules[ii]

      res <- list(
        timeX = timeX,
        timeZ = timeZ,
        modules = modules
      )
    }

    plot.RENDER <- function() {
      res <- plot_data()
      timeZ <- res$timeZ

      dimensions <- list()
      for (i in 1:ncol(timeZ)) {
        d <- list(
          range = range(timeZ),
          label = colnames(timeZ)[i],
          values = timeZ[, i]
        )
        dimensions[[i]] <- d
      }

      df <- data.frame(timeZ, check.names = FALSE)
      int.modules <- as.integer(factor(res$modules))

      plt <- plotly::plot_ly(
        df,
        source = "pcoords"
      ) %>%
        plotly::add_trace(
          type = "parcoords",
          line = list(
            color = int.modules,
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

    ## -----------------------------------------------------
    ## ------------------- gene table ----------------------
    ## -----------------------------------------------------
    table.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)

      timeX <- res$timeX
      feature1 <- gsub(";.*", ";...", rownames(timeX)) ## shorten
      sdx <- matrixStats::rowSds(timeX)

      df <- data.frame(
        module = res$modules,
        feature = feature1,
        SD = sdx,
        timeX,
        check.names = FALSE
      )

      symbol <- pgx$genes[rownames(timeX), "symbol"]
      if (mean(symbol == rownames(timeX), na.rm = TRUE) < 0.2) {
        df$symbol <- symbol
      }
      df <- df[order(-df$SD), ]

      cols <- c("module", "feature", "symbol", colnames(df))
      cols <- intersect(cols, colnames(df))
      df <- df[, cols]

      numeric.cols <- c("SD", colnames(timeX))
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
      dt <- table.RENDER()
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
