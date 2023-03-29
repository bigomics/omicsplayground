##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

featuremap_plot_table_geneset_map_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)

  info_text <- "<b>Geneset UMAP.</b> UMAP clustering of genesets colored by standard-deviation of log-expression(sd.X), or standard-deviation of the fold-change (sd.FC). The distance metric is covariance of the geneset expression. Genesets that are clustered nearby have high covariance."

  info_text_table <- "<b>Geneset table.</b> The contents of this table can be subsetted by selecting (by click&drag) on the <b>Geneset map</b> plot."

  plot.opts <- shiny::tagList(
    shiny::selectInput(ns("gsmap_nlabel"), "nr labels:",
      choices = c(0, 10, 20, 50, 100, 1000), selected = 20
    ),
    shiny::sliderInput(ns("gsmap_gamma"), "color gamma:",
      min = 0.1, max = 1.2, value = 0.4, step = 0.1
    ),
    shiny::radioButtons(ns("gsmap_colorby"), "color by:",
      ##choices = c("sd.X", "sd.FC", "mean.FC"),
      choices = c("sd.X", "sd.FC"),      
      selected = "sd.X", inline = TRUE
    )
  )

  div(
    PlotModuleUI(
      ns("gset_map"),
      title = "Geneset UMAP",
      label = "a",
      plotlib = "plotly",
      plotlib2 = "plotly",
      info.text = info_text,
      options = plot.opts,
      height = c(600, 700),
      width = c("auto", "100%"),
      download.fmt = c("png", "pdf")
    ),
    TableModuleUI(
      ns("gset_table"),
      info.text = info_text_table,
      height = c(280, TABLE_HEIGHT_MODAL),
      width = c("auto", "90%"),
      title = "Geneset table",
      label = "c"
    )
  )
}

featuremap_plot_table_geneset_map_server <- function(id,
                                                     pgx,
                                                     getGsetUMAP,
                                                     plotUMAP,
                                                     filter_gsets,
                                                     sigvar,
                                                     watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    selGsets <- shiny::reactive({
      shiny::req(pgx)
      db <- filter_gsets()
      gsets <- rownames(pgx$gsetX)
      gsets <- grep(paste0("^", db, ":"), gsets, value = TRUE)
      gsets
    })

    gsetUMAP.RENDER.SAVE <- shiny::reactive({

      pos <- getGsetUMAP()
      hilight <- NULL
      colgamma <- as.numeric(input$gsmap_gamma)

      F <- pgx.getMetaMatrix(pgx, level = "geneset")$fc
      F <- scale(F, center = FALSE)
      colorby <- input$gsmap_colorby
      if (colorby == "sd.FC") {
        fc <- (rowMeans(F**2))**0.2
      } else if (colorby == "mean.FC") {
        fc <- rowMeans(F)
      } else {
        cX <- pgx$gsetX - rowMeans(pgx$gsetX, na.rm = TRUE)
        fc <- sqrt(rowMeans(cX**2))
      }
      fc <- sign(fc) * abs(fc / max(abs(fc)))**colgamma

      ## filter on table
      hilight <- selGsets()
      nlabel <- as.integer(input$gsmap_nlabel)

      par(mfrow = c(1, 1))
      p <- plotUMAP(pos, fc, hilight,
        nlabel = nlabel, title = colorby,
        cex = 0.9, source = "", plotlib = "base"
      )
      p
    })

    plot_data <- shiny::reactive({
      pos <- getGsetUMAP()
      colnames(pos) <- c("x","y")

      hilight <- selGsets()
      colgamma <- as.numeric(input$gsmap_gamma)
      nlabel <- as.integer(input$gsmap_nlabel)
      colorby <- input$gsmap_colorby
      
      F <- pgx.getMetaMatrix(pgx, level = "geneset")$fc
      F <- scale(F, center = FALSE)
      if (colorby == "var.FC") {
        fc <- (rowMeans(F**2))**0.5
      } else if (colorby == "mean.FC") {
        fc <- rowMeans(F)
      } else {
        cX <- pgx$gsetX - rowMeans(pgx$gsetX, na.rm = TRUE)
        fc <- sqrt(rowMeans(cX**2))
      }
      fc <- sign(fc) * abs(fc / max(abs(fc)))**colgamma

      pd <- list(
        df = data.frame(pos, fc=fc),
        hilight = hilight,
        colgamma = colgamma,
        nlabel = nlabel,
        colorby = colorby
      )

    })
    
    render_gsetUMAP <- function() {

      pd  <- plot_data()
      pos <- pd$df[,c("x","y")]
      fc  <- setNames(pd$df[,"fc"], rownames(pd$df))
      hilight <- pd$hilight
      nlabel  <- pd$nlabel
      colorby <- pd$colorby
      
      ## filter on table
      p <- plotUMAP(
        pos,
        fc,
        hilight,
        nlabel = nlabel,
        title = colorby,
        cex = 1.2,
        source =  ns("geneset_filter"),
        plotlib = "plotly"
      ) %>%
        plotly::layout(dragmode = "select")
      p
    }

    gsetUMAP.RENDER <- function() {
      p <- render_gsetUMAP() %>%
        plotly::config(
          modeBarButtons = list(list("toImage", "zoom2d", "select2d", "resetScale2d"))
        ) %>%
        plotly_default()
      p
    }

    gsetUMAP.RENDER2 <- function() {
      p <- render_gsetUMAP() %>%
        plotly::config(
          modeBarButtons = list(list("toImage", "zoom2d", "select2d", "resetScale2d"))
        ) %>%
        plotly_modal_default()
      p
    }
    
    PlotModuleServer(
      "gset_map",
      plotlib = "plotly",
      plotlib2 = "plotly",
      func = gsetUMAP.RENDER,
      func2 = gsetUMAP.RENDER2,
      csvFunc = plot_data,
      pdf.width = 5, pdf.height = 5,
      add.watermark = watermark
    )

    # Table

    gsetTable.RENDER <- shiny::reactive({
      shiny::req(pgx)
      if (is.null(pgx$drugs)) {
        return(NULL)
      }

      pos <- getGsetUMAP()

      ## detect brush
      sel.gsets <- NULL
      b <- plotly::event_data("plotly_selected", source = ns("geneset_filter"))

      if (!is.null(b) & length(b)) {
        sel <- b$key
        sel.gsets <- rownames(pos)[rownames(pos) %in% sel]
      }

      pheno <- "tissue"
      pheno <- sigvar()
      is.fc <- FALSE
      if (pheno %in% colnames(pgx$samples)) {
        X <- pgx$gsetX - rowMeans(pgx$gsetX)
        y <- pgx$samples[, pheno]
        F <- do.call(cbind, tapply(1:ncol(X), y, function(i) {
          rowMeans(X[, i, drop = FALSE])
        }))
        is.fc <- FALSE
      } else {
        F <- pgx.getMetaMatrix(pgx, level = "geneset")$fc
        is.fc <- TRUE
      }

      if (!is.null(sel.gsets)) {
        sel.gsets <- intersect(sel.gsets, rownames(F))
        F <- F[sel.gsets, ]
      }
      F <- F[order(-rowMeans(F**2)), ]

      F <- cbind(sd.X = sqrt(rowMeans(F**2)), F)
      if (is.fc) colnames(F)[1] <- "sd.FC"
      F <- round(F, digits = 3)
      gs <- substring(rownames(F), 1, 100)
      df <- data.frame(geneset = gs, F, check.names = FALSE)

      DT::datatable(df,
        rownames = FALSE,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE, ## scrollY = TRUE,
          scrollY = 240,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    gsetTable.RENDER_modal <- shiny::reactive({
      dt <- gsetTable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "gset_table",
      func = gsetTable.RENDER,
      func2 = gsetTable.RENDER_modal,
      selector = "none"
    )
  })
}
