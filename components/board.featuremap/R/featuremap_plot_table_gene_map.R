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
      #
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
      if (is.null(pgx$version) | pgx$organism == "Human") {
        filtgenes <- unlist(lapply(sel, function(genes) playdata::FAMILIES[[genes]]))
      } else {
        filtgenes <- unlist(lapply(sel, function(genes) {
          if (genes == "<all>") {
            x <- pgx$genes$symbol
          } else {
            x <- playdata::FAMILIES[[genes]]
            x <- pgx$genes$symbol[match(x, pgx$genes$human_ortholog, nomatch = 0)]
          }
          return(x)
        }))
      }
      filtgenes
    })

    plot_data <- shiny::reactive({
      pos <- getUMAP()
      colnames(pos) <- c("x", "y")

      hilight <- filteredGenes()
      colorby <- input$umap_colorby
      nlabel <- as.integer(input$umap_nlabel)

      ## select on table filter
      FC <- playbase::pgx.getMetaMatrix(pgx)$fc
      FC <- scale(FC, center = FALSE)
      if (colorby == "sd.FC") {
        fc <- (rowMeans(FC**2))**0.5
      } else {
        cX <- pgx$X - rowMeans(pgx$X, na.rm = TRUE)
        fc <- sqrt(rowMeans(cX**2))
      }

      ## conform
      pos <- playbase::rename_by(pos, pgx$genes, "symbol")
      names(fc) <- pgx$genes$symbol[match(names(fc), rownames(pgx$genes), nomatch = 0)]
      fc <- fc[!duplicated(names(fc))]
      pos <- pos[!duplicated(rownames(pos)), , drop = FALSE]
      gg <- intersect(rownames(pos), names(fc))
      pos <- pos[gg, ]
      fc <- fc[gg]

      pd <- list(
        df = data.frame(pos, fc = fc),
        pos = pos,
        fc = fc,
        hilight = hilight,
        nlabel = nlabel,
        colorby = colorby
      )
      return(pd)
    })

    render_geneUMAP <- function(cex = 1, cex.label = 1) {
      pd <- plot_data()
      pos <- pd$pos
      fc <- pd$fc

      hilight <- pd$hilight
      nlabel <- pd$nlabel
      colorby <- pd$colorby
      colgamma <- as.numeric(input$umap_gamma)

      ## dim non hilighted genes
      fc <- sign(fc) * abs(fc / max(abs(fc), na.rm = TRUE))**colgamma
      if (length(setdiff(names(fc), hilight))) {
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

    # ================================================================================
    # ============================= Table server module ==============================
    # ================================================================================
    geneTable.RENDER <- shiny::reactive({
      shiny::req(pgx$X)
      X <- pgx$X
      pos <- getUMAP()
      pos <- playbase::rename_by(pos, pgx$genes, "symbol")
      sel.genes <- NULL
      ## detect brush
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

      pheno <- sigvar()
      is.fc <- FALSE
      X <- playbase::rename_by(X, pgx$genes, "symbol")
      if (any(pheno %in% colnames(pgx$samples))) {
        gg <- intersect(sel.genes, rownames(X))
        X <- X[gg, , drop = FALSE]
        X <- X - rowMeans(X)
        y <- pgx$samples[, pheno]
        FC <- do.call(cbind, tapply(1:ncol(X), y, function(i) {
          rowMeans(X[, i, drop = FALSE])
        }))
        is.fc <- FALSE
      } else {
        FC <- playbase::pgx.getMetaMatrix(pgx, level = "gene")$fc
        gg <- intersect(sel.genes, rownames(FC))
        if (length(gg) == 0) {
          FC <- playbase::rename_by(FC, pgx$genes, "symbol")
          gg <- intersect(sel.genes, rownames(FC))
        }
        FC <- FC[gg, , drop = FALSE]
        is.fc <- TRUE
      }
      FC <- FC[order(-rowMeans(FC**2)), , drop = FALSE]
      gene_table <- pgx$genes
      if (all(gene_table$human_ortholog == rownames(gene_table)) | all(is.na(gene_table$human_ortholog))) {
        gene_table_cols <- c("feature", "symbol", "gene_title")
      } else {
        gene_table_cols <- c("feature", "symbol", "human_ortholog", "gene_title")
      }

      # Retreive gene table with rownames (symbols)
      rowids <- pgx$genes$gene_name[match(rownames(FC), pgx$genes$symbol, nomatch = 0)]
      tt <- pgx$genes[rowids, gene_table_cols, drop = FALSE]
      tt <- apply(tt, MARGIN = 2, playbase::shortstring, n = 60)

      FC <- cbind(sd.X = sqrt(rowMeans(FC**2)), FC)
      if (is.fc) colnames(FC)[1] <- "sd.FC"
      FC <- round(FC, digits = 3)

      df <- data.frame(tt, FC,
        check.names = FALSE
      )

      # Remove feature column if it is the same as symbol column
      if (sum(df$feature %in% df$symbol) > nrow(df) * .8) {
        df$feature <- NULL
      }

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
