##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

featuremap_plot_gene_map_ui <- function(
    id,
    title,
    info.text,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  plot.opts <- shiny::tagList(
    shiny::selectInput(ns("umap_nlabel"), "nr labels:",
      c(0, 10, 20, 50, 100, 1000),
      selected = 50
    ),
    shiny::sliderInput(ns("umap_gamma"), "tweak colors:",
      min = 0.1, max = 1.2, value = 0.4, step = 0.1
    ),
    shiny::radioButtons(ns("umap_colorby"), "color by:",
      #
      choices = c("sd.X", "rms.FC"),
      selected = "sd.X", inline = TRUE
    )
  )

  PlotModuleUI(
    ns("gene_map"),
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

featuremap_table_gene_map_ui <- function(
    id,
    label = "",
    title,
    caption,
    info.text,
    height,
    width) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("gene_table"),
    info.text = info.text,
    height = height,
    width = width,
    caption = caption,
    title = title,
    label = "c"
  )
}

featuremap_plot_gene_map_server <- function(id,
                                            pgx,
                                            plotUMAP,
                                            sigvar,
                                            filter_genes,
                                            r_fulltable,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    getUMAP <- function() {
      pos <- pgx$cluster.genes$pos[["umap2d"]]
      colnames(pos) <- c("x", "y")
      pos
    }

    filteredGenes <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::validate(need(filter_genes(), "Please input at least one value in Annotate genes!"))
      sel <- filter_genes()
      filtgenes <- c()
      filtgenes <- unlist(lapply(sel, function(genes) playdata::FAMILIES[[genes]]))
      filtgenes
    })

    plot_data <- shiny::reactive({
      pos <- getUMAP()
      hilight <- filteredGenes()
      colorby <- input$umap_colorby
      nlabel <- as.integer(input$umap_nlabel)

      ## select on table filter
      F <- playbase::pgx.getMetaMatrix(pgx)$fc
      F <- scale(F, center = FALSE)
      if (colorby == "rms.FC") {
        fc <- sqrt(rowMeans(F**2,na.rm=TRUE))
      } else {
        cX <- pgx$X - rowMeans(pgx$X, na.rm = TRUE)
        fc <- sqrt(rowMeans(cX**2))
      }

      ## conform
      gg  <- intersect(rownames(pos), names(fc))
      pos <- pos[gg, ]
      fc  <- fc[gg]

      pd <- list(
        df = data.frame(pos, fc = fc),
        fc = fc,
        hilight = hilight,
        nlabel = nlabel,
        colorby = colorby
      )
    })

    render_geneUMAP <- function(cex = 1, cex.label = 1) {
      pd <- plot_data()
      pos <- getUMAP()
      fc <- pd$fc

      hilight  <- pd$hilight
      nlabel   <- pd$nlabel
      colorby  <- pd$colorby
      colgamma <- as.numeric(input$umap_gamma)

      ## dim non hilighted genes
      fc <- sign(fc) * abs(fc / max(abs(fc),na.rm=TRUE))**colgamma
      if(length(setdiff(names(fc),hilight))) {      
        fc[!names(fc) %in% hilight] <- NA
      }
      
      p <- plotUMAP(
        pos,
        fc,
        hilight,
        nlabel = nlabel,
        title = colorby,
        cex = cex,
        cex.label = cex.label,
        plotlib = "plotly",
        source = ns("gene_umap")
      ) %>%
        plotly::layout(
          dragmode = "select",
          margin = list(l = 5, r = 5, b = 5, t = 20)
        )
      p
    }

    geneUMAP.RENDER <- function() {
      p <- render_geneUMAP(cex = 1, cex.label = 1) %>%
        plotly::config(
          modeBarButtons = list(list("toImage", "zoom2d", "select2d", "resetScale2d"))
        ) %>%
        plotly_default()
      p
    }

    geneUMAP.RENDER2 <- function() {
      p <- render_geneUMAP(cex = 1.2, cex.label = 1.5) %>%
        plotly::config(
          modeBarButtons = list(list("toImage", "zoom2d", "select2d", "resetScale2d"))
        ) %>%
        plotly_modal_default()
      p
    }

    PlotModuleServer(
      "gene_map",
      plotlib = "plotly",
      plotlib2 = "plotly",
      func = geneUMAP.RENDER,
      func2 = geneUMAP.RENDER2,
      csvFunc = plot_data,
      pdf.width = 5, pdf.height = 5,
      add.watermark = watermark
    )

    #================================================================================
    #============================= Table server module ==============================
    #================================================================================    
    geneTable.RENDER <- shiny::reactive({
      shiny::req(pgx$X)

      ## detect brush
      pos <- getUMAP()
      sel.genes <- NULL
      b <- plotly::event_data("plotly_selected", source = ns("gene_umap"))
      if (!is.null(b) & length(b) > 0) {
        sel <- b$key
        sel.genes <- rownames(pos)[rownames(pos) %in% sel]
      } else {
        sel.genes <- rownames(pos)
      }

      if (!r_fulltable()) {
        if (!is.null(sel.genes)) {
          filt.genes <- filteredGenes()
          sel.genes_aux <- match(filt.genes |> stringr::str_to_upper(), sel.genes |> stringr::str_to_upper())
          sel.genes_aux <- sel.genes_aux[!is.na(sel.genes_aux)]
          sel.genes <- sel.genes[sel.genes_aux]
        } else {
          sel.genes <- filteredGenes()
        }
      }

      pheno <- "tissue"
      pheno <- sigvar()
      is.fc <- FALSE
      if (any(pheno %in% colnames(pgx$samples))) {
        gg <- intersect(sel.genes, rownames(pgx$X))
        X <- pgx$X[gg, , drop = FALSE]
        X <- X - rowMeans(X)
        y <- pgx$samples[, pheno]
        F <- do.call(cbind, tapply(1:ncol(X), y, function(i) {
          rowMeans(X[, i, drop = FALSE])
        }))
        is.fc <- FALSE
      } else {
        F <- playbase::pgx.getMetaMatrix(pgx, level = "gene")$fc
        gg <- intersect(sel.genes, rownames(F))
        F <- F[gg, , drop = FALSE]
        is.fc <- TRUE
      }
      F <- F[order(-rowMeans(F**2)), , drop = FALSE]

      tt <- playbase::shortstring(pgx$genes[rownames(F), "gene_title"], 60)
      tt <- as.character(tt)
      F <- cbind(sd.X = sqrt(rowMeans(F**2)), F)
      if (is.fc) colnames(F)[1] <- "rms.FC"
      F <- round(F, digits = 3)

      df <- data.frame(
        gene = rownames(F), title = tt, F,
        check.names = FALSE
      )

      DT::datatable(df,
        rownames = FALSE,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = NULL),
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = 240,
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    geneTable.RENDER_modal <- shiny::reactive({
      dt <- geneTable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "gene_table",
      func = geneTable.RENDER,
      func2 = geneTable.RENDER_modal,
      selector = "none"
    )
  }) ## moduleServer
}
