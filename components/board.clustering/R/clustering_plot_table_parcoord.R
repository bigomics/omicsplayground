##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##



## Annotate clusters ############

clustering_plot_parcoord_ui <- function(
  id,
  label = "",
  title,
  info.text,
  caption,
  height,
  width) {
  ns <- shiny::NS(id)

  parcoord_opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("hm_pcscale"), "Scale values", TRUE),
      "Scale expression values to mean=0 and SD=1.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
      ns("pltmod"),
      title = title,
      label = label,
      plotlib = "plotly",
      info.text = info.text,
      caption = caption,
      options = parcoord_opts,
      download.fmt = c("png", "pdf", "csv"),
      width = width,
      height = height
  )
  
}

clustering_table_parcoord_ui <- function(
  id,
  label = "",
  title,
  info.text,
  caption,
  height,
  width) {
    ns <- shiny::NS(id)
  
    TableModuleUI(
        ns("datasets"),
        info.text = info.text,
        height = height,
        caption = caption,
        width = width,
        title = title,
        label = "b"
    )
  
}

clustering_plot_table_parcoord_server <- function(id,
                                                     parcoord.matrix,
                                                     watermark = FALSE,
                                                     getTopMatrix) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    shiny::observeEvent(plotly::event_data("plotly_restyle", source = "pcoords"), {
      ## From: https://rdrr.io/cran/plotly/src/inst/examples/shiny/event_data_parcoords/app.R
      ##
      d <- plotly::event_data("plotly_restyle", source = "pcoords")
      ## what is the relevant dimension (i.e. variable)?
      dimension <- as.numeric(stringr::str_extract(names(d[[1]]), "[0-9]+"))
      ## If the restyle isn't related to a dimension, exit early.
      if (!length(dimension)) {
        return()
      }
      if (is.na(dimension)) {
        return()
      }

      pc <- parcoord.matrix()
      shiny::req(pc)
      ## careful of the indexing in JS (0) versus R (1)!
      dimension_name <- colnames(pc$mat)[[dimension + 1]]
      ## a given dimension can have multiple selected ranges
      ## these will come in as 3D arrays, but a list of vectors
      ## is nicer to work with
      info <- d[[1]][[1]]
      if (length(dim(info)) == 3) {
        parcoord.ranges[[dimension_name]] <- lapply(seq_len(dim(info)[2]), function(i) info[, i, ])
      } else {
        parcoord.ranges[[dimension_name]] <- list(as.numeric(info))
      }
    })


    parcoord.ranges <- shiny::reactiveValues()

    parcoord.matrix <- shiny::reactive({
      filt <- getTopMatrix()
      shiny::req(filt)
      zx <- filt$mat[, ]
      if (input$hm_pcscale) {
        zx <- t(scale(t(zx)))
      }
      rr <- shiny::isolate(shiny::reactiveValuesToList(parcoord.ranges))
      nrange <- length(rr)
      for (i in names(rr)) parcoord.ranges[[i]] <- NULL
      zx <- round(zx, digits = 4)
      list(mat = zx, clust = filt$idx)
    })

    parcoord.selected <- shiny::reactive({
      mat <- parcoord.matrix()$mat
      clust <- parcoord.matrix()$clust
      shiny::req(mat)
      keep <- TRUE
      for (i in names(parcoord.ranges)) {
        range_ <- parcoord.ranges[[i]]
        range_ <- range_[sapply(range_, length) > 0]
        if (length(range_) > 0) {
          keep_var <- FALSE
          for (j in seq_along(range_)) {
            rng <- range_[[j]]
            keep_var <- keep_var | dplyr::between(mat[, i], min(rng), max(rng))
          }
          keep <- keep & keep_var
        }
      }
      list(mat = mat[keep, , drop = FALSE], clust = clust[keep])
    })

    plot_data <- function() {
      parcoord.matrix()
    }

    parcoord.RENDER <- function() {
      pc <- plot_data()
      zx <- pc$mat
      ## build dimensions
      dimensions <- list()
      for (i in 1:ncol(zx)) {
        dimensions[[i]] <- list(
          range = c(min(zx[, i]), max(zx[, i])),
          visible = TRUE,
          label = colnames(zx)[i],
          values = zx[, i]
        )
      }

      clust.id <- as.integer(factor(pc$clust))
      table(clust.id)

      df <- data.frame(clust.id = clust.id, zx)
      klrpal <- rep(RColorBrewer::brewer.pal(8, "Set2"), 99)
      klrpal <- klrpal[1:max(clust.id)]
      ## klrpal <- setNames(klrpal, sort(unique(clust.id)))
      klrpal2 <- lapply(1:length(klrpal), function(i) c((i - 1) / (length(klrpal) - 1), klrpal[i]))

      plt <- plotly::plot_ly(df, source = "pcoords") %>%
        plotly::add_trace(
          type = "parcoords",
          line = list(
            color = ~clust.id,
            colorscale = klrpal2,
            cmin = min(clust.id), cmax = max(clust.id),
            showscale = FALSE
            ## reversescale = TRUE
          ),
          dimensions = dimensions
        )
      plt <- plt %>%
        plotly::layout(margin = list(l = 60, r = 60, t = 0, b = 30)) %>%
        ## config(displayModeBar = FALSE) %>%
        ## config(modeBarButtonsToRemove = setdiff(all.plotly.buttons,"toImage") ) %>%
        plotly::config(toImageButtonOptions = list(format = "svg", width = 900, height = 350, scale = 1.2)) %>%
        plotly::config(displaylogo = FALSE) %>%
        plotly::event_register("plotly_restyle")

      plt
    }

    parcoord.RENDER_MODAL <- function() {
      parcoord.RENDER() %>%
        plotly_modal_default()
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = parcoord.RENDER,
      func2 = parcoord.RENDER_MODAL,
      csvFunc = plot_data,
      res = c(90, 170),
      pdf.width = 8,
      pdf.height = 8,
      add.watermark = watermark
    )

    ## Table ------------------------------------------------------------------
    parcoord_table.RENDER <- function() {
      parcoord <- parcoord.selected()

      mat <- parcoord$mat
      clust <- parcoord$clust
      df <- data.frame(cluster = clust, mat, check.names = FALSE)
      numeric.cols <- 2:ncol(df)
      DT::datatable(
        df,
        rownames = TRUE, ## escape = c(-1,-2),
        extensions = c("Buttons", "Scroller"),
        plugins = 'scrollResize',        
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip", ## buttons = c('copy','csv','pdf'),
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE, ## scrollY = TRUE,
          ## scrollY = 170,
          scrollY = "23vh",
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    parcoord_table.RENDER_modal <- function() {
      dt <- parcoord_table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    }

    TableModuleServer(
      "datasets",
      func = parcoord_table.RENDER,
      func2 = parcoord_table.RENDER_modal,
      selector = "none"
    )
  })
}
