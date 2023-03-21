##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

featuremap_plot_gene_map_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)

  info_text <- "<b>Gene map.</b> UMAP clustering of genes colored by standard-deviation (sd.X), variance (var.FC) or mean of fold-change (mean.FC). The distance metric is covariance. Genes that are clustered nearby exihibit high covariance and may have similar biological function."

  info_text_table <- "<b>Gene table.</b> The contents of this table can be subsetted by selecting (by click&drag) on the <b>Gene map</b> plot."

  plot.opts <- shiny::tagList(
    shiny::selectInput(ns("umap_nlabel"), "nr labels:",
      c(0, 10, 20, 50, 100, 1000),
      selected = 50
    ),
    shiny::sliderInput(ns("umap_gamma"), "color gamma:",
      min = 0.1, max = 1.2, value = 0.4, step = 0.1
    ),
    shiny::radioButtons(ns("umap_colorby"), "color by:",
      choices = c("sd.X", "var.FC", "mean.FC"),
      selected = "sd.X", inline = TRUE
    )
  )

  div(
    PlotModuleUI(
      ns("gene_map"),
      title = "Gene UMAP",
      label = "a",
      # outputFunc = function(x, width, height) {
      #   plotOutput(x,
      #     brush = ns("geneUMAP_brush"), width = width,
      #     height = height
      #   )
      # },
      plotlib = "plotly",
      plotlib2 = "plotly",
      info.text = info_text,
      options = plot.opts,
      height = c(600, 700),
      width = c("auto", "100%"),
      download.fmt = c("png", "pdf")
    ),
    TableModuleUI(
      ns("datasets"),
      info.text = info_text_table,
      height = c(280, TABLE_HEIGHT_MODAL),
      width = c("auto", "90%"),
      title = "Gene table",
      label = "c"
    )
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
    selGenes <- shiny::reactive({
      shiny::req(pgx)
      sel <- filter_genes()
      selgenes <- FAMILIES[[sel]]
      selgenes
    })

    geneUMAP.RENDER <- shiny::reactive({
      shiny::req(pgx)

      pos <- getGeneUMAP()
      hilight <- NULL
      colgamma <- as.numeric(input$umap_gamma)

      ## select on table filter
      F <- pgx.getMetaMatrix(pgx)$fc
      F <- scale(F, center = FALSE)
      colorby <- input$umap_colorby
      if (colorby == "sd.FC") {
        fc <- (rowMeans(F**2))**0.2
      } else if (colorby == "mean.FC") {
        fc <- rowMeans(F)
      } else {
        ## sdX
        cX <- pgx$X - rowMeans(pgx$X, na.rm = TRUE)
        fc <- sqrt(rowMeans(cX**2))
      }
      fc <- sign(fc) * abs(fc / max(abs(fc)))**colgamma
      hilight <- names(sort(-fc))

      ## filter on table
      hilight <- selGenes()
      nlabel <- as.integer(input$umap_nlabel)

      p <- plotUMAP(pos, fc, hilight,
        nlabel = nlabel, title = colorby,
        cex = 0.9, source = "", plotlib = "base"
      )
      p
    })

    geneUMAP.RENDER2 <- shiny::reactive({
      shiny::req(pgx)

      ## pos <- pgx$cluster.genes$pos[['umap2d']]
      pos <- getGeneUMAP()
      hilight <- NULL
      colgamma <- as.numeric(input$umap_gamma)

      ## select on table filter
      F <- pgx.getMetaMatrix(pgx)$fc
      F <- scale(F, center = FALSE)
      colorby <- input$umap_colorby
      if (colorby == "var.FC") {
        fc <- (rowMeans(F**2))**0.2
      } else if (colorby == "mean.FC") {
        fc <- rowMeans(F)
      } else {
        cX <- pgx$X - rowMeans(pgx$X, na.rm = TRUE)
        fc <- sqrt(rowMeans(cX**2))
      }
      fc <- sign(fc) * abs(fc / max(abs(fc)))**colgamma

      ## filter on table
      hilight <- selGenes()
      nlabel <- as.integer(input$umap_nlabel)

      p <- plotUMAP(pos, fc, hilight,
        nlabel = nlabel, title = colorby,
        cex = 1.2, source = "", plotlib = "plotly"
      )
      p
    })

    PlotModuleServer(
      "gene_map",
      plotlib = "plotly",
      plotlib2 = "plotly",
      func = geneUMAP.RENDER2,
      func2 = geneUMAP.RENDER2,
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
      ## b <- input$ftmap-geneUMAP_brush  ## ugly??
      b <- NULL
      b <- input[["geneUMAP_brush"]] ## ugly??
      if (!is.null(b) & length(b) > 0) {
        sel <- which(pos[, 1] > b$xmin & pos[, 1] < b$xmax &
          pos[, 2] > b$ymin & pos[, 2] < b$ymax)
        sel.genes <- rownames(pos)[sel]
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
        F <- pgx.getMetaMatrix(pgx, level = "gene")$fc
        is.fc <- TRUE
      }

      if (!is.null(sel.genes)) {
        sel.genes <- intersect(sel.genes, rownames(F))
        F <- F[sel.genes, , drop = FALSE]
      }
      F <- F[order(-rowMeans(F**2)), ]

      tt <- shortstring(pgx$genes[rownames(F), "gene_title"], 60)
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
          scrollX = TRUE, ## scrollY = TRUE,
          scrollY = "70vh",
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
      "datasets",
      func = geneTable.RENDER,
      func2 = geneTable.RENDER_modal,
      selector = "none"
    )
  })
}
