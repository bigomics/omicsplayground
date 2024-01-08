##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


clustering_plot_phenoplot_ui <- function(
    id,
    title,
    info.text,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  phenoplot.opts <- shiny::tagList(
    shiny::checkboxInput(ns("showlabels"), "Show group labels", TRUE)
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    ##    plotlib = "plotly",
    plotlib = "generic",
    outputFunc = function(...) uiOutput(..., fill = TRUE),
    info.text = info.text,
    caption = caption,
    options = phenoplot.opts,
    ##    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
  )
}

clustering_plot_phenoplot_server <- function(id,
                                             pgx,
                                             selected_phenotypes,
                                             clustmethod,
                                             selected_samples,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- reactive({
      pgx <- pgx
      shiny::req(pgx$Y)

      ## get t-SNE positions
      clustmethod1 <- paste0(clustmethod(), "2d")
      pos <- pgx$cluster$pos[[clustmethod1]]
      colnames(pos) <- c("x", "y")
      jj <- selected_samples()
      kk <- selected_phenotypes()
      pos <- pos[jj, ]
      Y <- pgx$Y[jj, kk, drop = FALSE]
      ## complete dataframe for downloading
      df <- data.frame(pos, Y)
      return(df)
    })

    create_plots <- function(cex = 1) {
      pd <- plot_data()
      showlabels <- input$showlabels
      pheno <- selected_phenotypes()
      Y <- pd[, pheno, drop = FALSE]
      pos <- pd[, c("x", "y")]

      ## points size depending on how many points we have
      ncex <- cut(nrow(pos), breaks = c(-1, 40, 200, 1000, 1e10))
      cex1 <- 0.8 * cex * c(1.8, 1.3, 0.8, 0.5)[ncex]
      cex1 <- cex1 * ifelse(length(pheno) > 6, 0.8, 1)
      cex1 <- cex1 * ifelse(length(pheno) > 12, 0.8, 1)

      plt <- list()
      for (i in 1:min(20, length(pheno))) {
        ## ------- set colors
        colvar <- factor(Y[, pheno[i]])
        colvar[which(colvar %in% c(NA, "", " ", "NA", "na"))] <- NA
        colvar <- factor(as.character(colvar))
        klrpal <- playdata::COLORS
        klr1 <- klrpal[colvar]
        klr1 <- paste0(gplots::col2hex(klr1), "99")
        jj <- which(is.na(klr1))
        if (length(jj)) klr1[jj] <- "#AAAAAA22"
        tt <- tolower(pheno[i])

        ## ------- start plot
        p <- playbase::pgx.scatterPlotXY.PLOTLY(
          pos,
          var = colvar,
          col = klrpal,
          cex = cex1,
          xlab = "",
          ylab = "",
          title = tt,
          cex.title = cex * 1.2,
          cex.clust = cex * 1.1,
          label.clusters = showlabels
        ) %>%
          plotly_default() %>%
          plotly::layout(
            plot_bgcolor = "#f8f8f8",
            margin = list(l = 0, r = 0, b = 0, t = 20) # lrbt
          )

        plt[[i]] <- p
      }
      return(plt)
    }

    plotly.RENDER <- function() {
      plt <- create_plots(cex = 0.85)
      nc <- floor(sqrt(length(plt)))
      cw <- 12 / nc
      page <- bslib::layout_columns(col_widths = cw, !!!plt)
      return(page)
    }

    plotly_modal.RENDER <- function() {
      plt <- create_plots(cex = 1.3)
      nc <- ceiling(sqrt(length(plt)))
      cw <- 12 / nc
      page <- bslib::layout_columns(col_widths = cw, !!!plt)
      return(page)
    }

    PlotModuleServer(
      "pltmod",
      #      plotlib = "plotly",
      plotlib = "generic",
      renderFunc = shiny::renderUI,
      func = plotly.RENDER,
      func2 = plotly_modal.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(85), ## resolution of plots
      pdf.width = 6,
      pdf.height = 9,
      add.watermark = watermark
    )
  })
}
