##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

featuremap_plot_geneset_map_ui <- function(
    id,
    label = "",
    title,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
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
      choices = c("sd.X", "rms.FC"),
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
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
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
                                                     plotUMAP,
                                                     filter_gsets,
                                                     sigvar,
                                                     r_fulltable,
                                                     watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## In this setup, PlotModule and TableModule servers are combined
    ## into one. Logically this might be preferred but is also
    ## convenient if PlotModule and TableModule share functions or
    ## share data. Notice the UI for table and plot is still seperate
    ## to allow flexibility in the placement of the widget. Both use
    ## the same module server id.
    ##
    ns <- session$ns
    sel.row <- 1
    getUMAP <- function() {
      pos <- pgx$cluster.gsets$pos[["umap2d"]]
      colnames(pos) <- c("x", "y")
      pos
    }

    filteredGsets <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::validate(shiny::need(
        filter_gsets(),
        tspan("Please input at least one value in Annotate genesets!", js = FALSE)
      ))
      db <- filter_gsets()
      gsets <- rownames(pgx$gsetX)
      if (!"<all>" %in% db) {
        filt_genesets <- unlist(lapply(db, function(geneset) {
          grep(paste0("^", geneset, ":"), gsets, value = TRUE)
        }))
        gsets <- filt_genesets
      }
      gsets
    })

    plot_data <- shiny::reactive({
      pos <- getUMAP()
      colnames(pos) <- c("x", "y")

      hilight <- filteredGsets()
      nlabel <- as.integer(input$gsmap_nlabel)

      F <- playbase::pgx.getMetaMatrix(pgx, level = "geneset")$fc
      F <- scale(F, center = FALSE)
      fc <- sqrt(rowMeans(F**2, na.rm = TRUE))

      ## conform
      gg <- intersect(rownames(pos), names(fc))
      pos <- pos[gg, ]
      fc <- fc[gg]

      pd <- list(
        df = data.frame(pos, fc = fc),
        fc = fc,
        pos = pos,
        hilight = hilight,
        nlabel = nlabel
      )
    })

    render_gsetUMAP <- function(cex = 1, cex.label = 1) {
      pd <- plot_data()
      pos <- pd$pos
      fc <- pd$fc
      hilight <- pd$hilight
      nlabel <- pd$nlabel
      colgamma <- as.numeric(input$gsmap_gamma)
      fc <- sign(fc) * abs(fc / max(abs(fc), na.rm = TRUE))**colgamma
      if (length(setdiff(names(fc), hilight))) {
        fc[!names(fc) %in% hilight] <- NA
      }

      ## filter on table
      p <- plotUMAP(
        pos,
        fc,
        hilight,
        nlabel = nlabel,
        title = "rms.FC",
        cex = cex,
        cex.label = cex.label,
        source = ns("geneset_umap"),
        plotlib = "plotly"
      ) %>%
        plotly::layout(
          dragmode = "select",
          margin = list(l = 5, r = 5, b = 5, t = 20)
        )
      p
    }

    gsetUMAP.RENDER <- function() {
      p <- render_gsetUMAP(cex = 1, cex.label = 0.9) %>%
        plotly::config(
          modeBarButtons = list(list("toImage", "zoom2d", "select2d", "resetScale2d"))
        ) %>%
        plotly_default()
      p
    }

    gsetUMAP.RENDER2 <- function() {
      p <- render_gsetUMAP(cex = 1.2, cex.label = 1.3) %>%
        plotly::config(
          modeBarButtons = list(list("toImage", "zoom2d", "select2d", "resetScale2d"))
        ) %>%
        plotly_modal_default()
      p
    }

    plotmodule <- PlotModuleServer(
      "gset_map",
      plotlib = "plotly",
      plotlib2 = "plotly",
      func = gsetUMAP.RENDER,
      func2 = gsetUMAP.RENDER2,
      csvFunc = plot_data,
      pdf.width = 5, pdf.height = 5,
      add.watermark = watermark
    )

    table_data <- shiny::reactive({
      shiny::req(pgx$X)
      pos <- getUMAP()

      ## detect brush
      sel.gsets <- NULL
      b <- plotly::event_data("plotly_selected", source = ns("geneset_umap"))

      if (!is.null(b) & length(b)) {
        sel <- b$key
        sel.gsets <- rownames(pos)[rownames(pos) %in% sel]
      } else {
        sel.gsets <- rownames(pos)
      }

      if (!is.null(sel.gsets)) {
        filt.gsets <- filteredGsets()
        sel.gsets <- intersect(sel.gsets, filt.gsets)
      }

      pheno <- sigvar()
      is.fc <- FALSE
      if (any(pheno %in% colnames(pgx$samples))) {
        gg <- intersect(sel.gsets, rownames(pgx$gsetX))
        X <- pgx$gsetX[gg, , drop = FALSE]
        X <- X - rowMeans(X, na.rm = TRUE)
        y <- pgx$samples[, pheno]
        if (nrow(X) == 1) {
          F <- tapply(1:ncol(X), y, function(i) {
            rowMeans(X[, c(i, i), drop = FALSE])
          })
          F <- data.frame(t(F))
          rownames(F) <- gg
        } else {
          F <- do.call(cbind, tapply(1:ncol(X), y, function(i) {
            rowMeans(X[, c(i, i), drop = FALSE])
          }))
        }
        is.fc <- FALSE
      } else {
        F <- playbase::pgx.getMetaMatrix(pgx, level = "geneset")$fc
        gg <- intersect(sel.gsets, rownames(F))
        F <- F[gg, , drop = FALSE]
        is.fc <- TRUE
      }

      F <- F[order(-rowMeans(F**2, na.rm = TRUE)), , drop = FALSE]
      F <- cbind(sd.X = sqrt(rowMeans(F**2)), F)
      if (is.fc) colnames(F)[1] <- "rms.FC"
      F <- round(F, digits = 3)
      gs.db <- sub(":.*", "", rownames(F))
      gs <- sub(".*[:]", "", rownames(F))
      gs <- substring(gs, 1, 100)
      df <- data.frame(DB = gs.db, geneset = gs, F, check.names = FALSE)
      rownames(df) <- rownames(F)
      return(df)
    })

    # Table
    gsetTable.RENDER <- shiny::reactive({
      df <- table_data()
      gset_link <- playbase::wrapHyperLink(
        rep_len("<i class='fa-solid fa-arrow-up-right-from-square'></i>", nrow(df)),
        rownames(df)
      ) |> HandleNoLinkFound(
        NoLinkString = "<i class='fa-solid fa-arrow-up-right-from-square'></i>",
        SubstituteString = "<i class='fa-solid fa-arrow-up-right-from-square blank_icon'></i>"
      )
      df$geneset <- paste(df$geneset, "&nbsp;", gset_link)


      DT::datatable(df,
        rownames = FALSE,
        escape = c(-1, -2),
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = NULL),
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE, #
          scrollY = 240,
          scrollResize = TRUE,
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

    table_data_csv <- function() {
      df <- table_data()
      return(df)
    }

    tablemodule <- TableModuleServer(
      "gset_table",
      func = gsetTable.RENDER,
      func2 = gsetTable.RENDER_modal,
      csvFunc = table_data_csv,
      selector = "none"
    )


    ## combined module return value for any downstream connections
    list(
      plotmodule = plotmodule,
      tablemodule = tablemodule
    )
  })
}
