##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

loading_tsne_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("pltmod"),
    outputFunc = plotly::plotlyOutput,
    outputFunc2 = plotly::plotlyOutput,
    info.text = info.text,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    caption = caption,
    height = height,
    label = label,
    title = title
  )
}

loading_tsne_server <- function(id, pgx.dirRT, info.table,
                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot_data <- shiny::reactive({

      pgx.dir <- pgx.dirRT()
      info.table <- info.table()
      
      dbg("[loading_tsne_server] reacted! id =",id)
      
      tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")
      ## pgx.files <- sub("[.]pgx$", "", dir(pgx.dir, pattern = ".pgx$"))
      pgx.files <- info.table$dataset
      
      pos <- NULL
      if (file.exists(tsne.file)) {
        pos <- read.csv(tsne.file, row.names = 1)
        dim(pos)
        pos.files <- unique(gsub("^\\[|\\].*", "", rownames(pos)))
        if (!all(pgx.files %in% pos.files)) {
          ## missing pgx files.. need recompute
          pos <- NULL
        } else {
          ## OK
        }
      }

      ## if no t-SNE file exists, we need to calculate it
      if (is.null(pos)) {
        dbg("[loading_tsneplot.R] recalculating tSNE positions...")
        F <- data.table::fread(file.path(pgx.dir, "datasets-allFC.csv"))
        F <- as.matrix(F[, -1], rownames = F[[1]])
        dim(F)
        ## F[is.na(F)] <- 0
        F <- apply(F, 2, rank, na.last = FALSE)
        corF <- cor(F, use = "pairwise") ## slow...
        dim(corF)
        px <- max(min(30, floor(ncol(corF) / 4)), 1)
        dbg("[loading_tsne_server] plot_data : dim(corF) = ", dim(corF))
        dbg("[loading_tsne_server] plot_data : px = ", px)        
        
        pos <- try(Rtsne::Rtsne(1 - abs(corF),
          perplexity = px,
          check_duplicates = FALSE, is_distance = TRUE
        )$Y)

        if("try-error" %in% class(pos)) {
          pos <- svd(corF)$u[,1:2]
          rownames(pos) <- colnames(F)
        }
        
        ## pos <- umap::umap(1-corF)$layout
        pos <- round(pos, digits = 4)
        rownames(pos) <- colnames(F)
        colnames(pos) <- c("x", "y")        
        write.csv(pos, file = tsne.file)
      }

      ## filter out non-existing entries
      pos.pgx <- gsub("^\\[|\\].*", "", rownames(pos))
      table(pos.pgx %in% pgx.files)
      setdiff(pos.pgx, pgx.files)
      pos <- pos[which(pos.pgx %in% pgx.files), , drop = FALSE]
      dim(pos)

      ##
      dset <- gsub("^\\[|\\].*", "", rownames(pos))
      comparison <- gsub("^.*\\]", "", rownames(pos))
      df <- data.frame(pos, dataset = dset, comparison = comparison)

      ## compute medioid of datasets
      dpos <- apply(pos, 2, function(x) tapply(x, dset, median, na.rm = TRUE))
      if (length(unique(dset)) == 1) {
        dpos <- matrix(dpos, nrow = 1, ncol = 2)
        rownames(dpos) <- dset[1]
      }
      colnames(dpos) <- c("x", "y")

      pdata <- list(
        df = df,
        pos = dpos
      )
      
      return(pdata)
    })

    plot.RENDER <- function() {

      pdata <- plot_data()
      shiny::req(pdata)      
      pos <- pdata[['pos']]
      df <- pdata[['df']]
      
      marker_size <- ifelse( nrow(df) > 50, 5, 10)
      
      fig <- plotly::plot_ly(
        data = df,
        x = ~x,
        y = ~y,
        text = ~ paste("Dataset:", dataset, "<br>Comparison:", comparison),
        color = ~dataset,
        ## colors = omics_pal_c(palette = "brand_blue")(100),
        marker = list(
          size = marker_size,
          line = list(
            color = omics_colors("super_dark_grey"),
            width = 1.0
          )
        )
      )

      dy <- diff(range(pos[,"y"]))
      dbg("[loading_tsneplot.R] range.y=",dy)
      
      fig <- fig %>%
        plotly::add_annotations(
          x = pos[,"x"],
          y = pos[,"y"],
          text = rownames(pos),
          xref = "x",
          yref = "y",          
          ## textposition = 'top',
          xanchor = "middle",
          yanchor = "bottom",
          yshift = 0.02*dy,
#          ax = 0,
#          ay = -0.05 * dy,
          showarrow = FALSE
        )

      fig <- fig %>%
        plotly::layout(
          showlegend = FALSE,
          xaxis = list(
            title = "tsne-x",
            zeroline = FALSE,
            showticklabels = FALSE
          ),
          yaxis = list(
            title = "tsne-y",
            zeroline = FALSE,            
            showticklabels = FALSE
          )
        )

      fig
    }

    modal_plot.RENDER <- function() {
      pdata <- plot_data()
      shiny::req(pdata)
      df <- pdata[['df']]
      marker_size <- ifelse( nrow(df) > 50, 6, 12)
      p <- plot.RENDER() %>%
        plotly::layout(
          showlegend = TRUE,
          font = list(
            size = 16
          )
        )
      p <- plotly::style(p, marker.size = marker_size)
      p
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      plotlib2 = "plotly",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      renderFunc = plotly::renderPlotly,
      renderFunc2 = plotly::renderPlotly,
      ## res = c(100,300)*1,              ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      ## label = label, title = "t-SNE clustering",
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
