##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

featuremap_plot_geneset_map_ui <- function(
  id,
  label = "",
  title,
  info.text,
  caption,
  height,
  width) {
  ns <- shiny::NS(id)

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

  PlotModuleUI(
      ns("gset_map"),
      title = title,
      label = "a",
      plotlib = "plotly",
      plotlib2 = "plotly",
      info.text = info.text,
      caption = caption,
      options = plot.opts,
      height = height,
      width = width,
      download.fmt = c("png", "pdf")
  )
}

featuremap_table_geneset_map_ui <- function(
  id,
  label = "",
  title,
  info.text,
  caption,
  height,
  width) {
  ns <- shiny::NS(id)

  TableModuleUI(
      ns("gset_table"),
      info.text = info.text,
      height = height,
      caption = caption,
      width = width,
      title = title,
  
      label = "c"
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

    plot_data <- shiny::reactive({
      pos <- getGsetUMAP()
      colnames(pos) <- c("x","y")

      hilight <- selGsets()
      colgamma <- as.numeric(input$gsmap_gamma)
      nlabel <- as.integer(input$gsmap_nlabel)
      colorby <- input$gsmap_colorby
      
      F <- playbase::pgx.getMetaMatrix(pgx, level = "geneset")$fc
      F <- scale(F, center = FALSE)
      if (colorby == "sd.FC") {
        fc <- (rowMeans(F**2))**0.5
      } else {
        cX <- pgx$gsetX - rowMeans(pgx$gsetX, na.rm = TRUE)
        fc <- sqrt(rowMeans(cX**2))
      }
      fc <- sign(fc) * abs(fc / max(abs(fc)))**colgamma

      ## conform
      gg <- intersect(rownames(pos), names(fc))
      pos <- pos[gg,]
      fc <- fc[gg]
      
      pd <- list(
        df = data.frame(pos, fc=fc),
        fc = fc,  
        hilight = hilight,
        colgamma = colgamma,
        nlabel = nlabel,
        colorby = colorby
      )

    })
    
    render_gsetUMAP <- function(cex=1, cex.label=1) {

      pd  <- plot_data()
      pos <- pd$df[,c("x","y")]
      fc <- pd$fc
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
        cex = cex,
        cex.label = cex.label,
        source =  ns("geneset_filter"),
        plotlib = "plotly"
      ) %>%
        plotly::layout(
          dragmode = "select",
          margin = list(l = 5, r = 5, b = 5, t = 20)                    
        )
      p
    }

    gsetUMAP.RENDER <- function() {
      p <- render_gsetUMAP(cex=1, cex.label=0.9) %>%
        plotly::config(
          modeBarButtons = list(list("toImage", "zoom2d", "select2d", "resetScale2d"))
        ) %>%
        plotly_default()
      p
    }

    gsetUMAP.RENDER2 <- function() {
      p <- render_gsetUMAP(cex=1.2, cex.label=1.3) %>%
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
        F <- playbase::pgx.getMetaMatrix(pgx, level = "geneset")$fc
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
