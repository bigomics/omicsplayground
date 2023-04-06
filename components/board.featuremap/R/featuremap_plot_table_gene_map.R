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
    shiny::sliderInput(ns("umap_gamma"), "color gamma:",
      min = 0.1, max = 1.2, value = 0.4, step = 0.1
    ),
    shiny::radioButtons(ns("umap_colorby"), "color by:",
      ##choices = c("sd.X", "var.FC", "mean.FC"),
      choices = c("sd.X", "sd.FC"),
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
                                            getGeneUMAP,
                                            plotUMAP,
                                            sigvar,
                                            filter_genes,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    selGenes <- shiny::reactive({
      shiny::req(pgx)
      sel <- filter_genes()
      selgenes <- FAMILIES[[sel]]
      selgenes
    })

    plot_data <- shiny::reactive({

      ## pos <- pgx$cluster.genes$pos[['umap2d']]
      pos <- getGeneUMAP()
      colnames(pos) <- c("x","y")
      hilight <- selGenes()
      colgamma <- as.numeric(input$umap_gamma)
      colorby <- input$umap_colorby
      nlabel <- as.integer(input$umap_nlabel)

      ## select on table filter
      F <- playbase::pgx.getMetaMatrix(pgx)$fc
      F <- scale(F, center = FALSE)
      if (colorby == "sd.FC") {
        fc <- (rowMeans(F**2))**0.5
      } else {
        cX <- pgx$X - rowMeans(pgx$X, na.rm = TRUE)
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
    
    render_geneUMAP <- function(cex.label=1) {

      pd  <- plot_data()
      pos <- pd$df[,c("x","y")]
      fc  <- pd$fc
      hilight <- pd$hilight
      nlabel  <- pd$nlabel
      colorby <- pd$colorby

      p <- plotUMAP(
        pos,
        fc,
        hilight,
        nlabel = nlabel,
        title = colorby,
        cex = 1.2,
        cex.label = cex.label,
        plotlib = "plotly",
        source = ns("gene_filter")
      ) %>%
        plotly::layout(dragmode = "select")
      p
    }

    geneUMAP.RENDER <- function() {
      p <- render_geneUMAP(cex.label=1) %>%
        plotly::config(
          modeBarButtons = list(list("toImage", "zoom2d", "select2d", "resetScale2d"))
        ) %>%
        plotly_default()
      p
    }

    geneUMAP.RENDER2 <- function() {
      p <- render_geneUMAP(cex.label=1.5) %>%
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

    # Table

    geneTable.RENDER <- shiny::reactive({

      shiny::req(pgx)
      if (is.null(pgx$drugs)) {
        return(NULL)
      }

      pos <- getGeneUMAP()

      ## detect brush
      sel.genes <- NULL
      b <- plotly::event_data("plotly_selected", source = ns("gene_filter"))
      if (!is.null(b) & length(b) > 0) {
        sel <- b$key
        sel.genes <- rownames(pos)[rownames(pos) %in% sel]
      }

      pheno <- "tissue"
      pheno <- sigvar()
      is.fc <- FALSE
      if (pheno %in% colnames(pgx$samples)) {
        X <- pgx$X - rowMeans(pgx$X)
        y <- pgx$samples[, pheno]
        F <- do.call(cbind, tapply(1:ncol(X), y, function(i) {
          rowMeans(X[, i, drop = FALSE])
        }))
        is.fc <- FALSE
      } else {
        F <- playbase::pgx.getMetaMatrix(pgx, level = "gene")$fc
        is.fc <- TRUE
      }

      if (!is.null(sel.genes)) {
        sel.genes <- intersect(sel.genes, rownames(F))
        F <- F[sel.genes, , drop = FALSE]
      }
      F <- F[order(-rowMeans(F**2)), ]

      tt <- playbase::shortstring(pgx$genes[rownames(F), "gene_title"], 60)
      tt <- as.character(tt)
      F <- cbind(sd.X = sqrt(rowMeans(F**2)), F)
      if (is.fc) colnames(F)[1] <- "sd.FC"
      F <- round(F, digits = 3)
      df <- data.frame(
        gene = rownames(F), title = tt, F,
        check.names = FALSE
      )

      DT::datatable(df,
        rownames = FALSE,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = 240,
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
      
  })  ## moduleServer
}
