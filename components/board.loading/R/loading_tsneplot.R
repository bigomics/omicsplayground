##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

loading_tsne_ui <- function(id,
                            label = "",
                            height,
                            width,
                            card_height = "100%") {
  ns <- shiny::NS(id)

  info_text <- paste0("<b>Similarity clustering</b> of fold-change signatures colored by data sets using t-SNE. Each dot corresponds to a specific comparison. Signatures/datasets that are clustered closer together, are more similar.")

  PlotModuleUI(
    ns("pltmod"),
    outputFunc = plotly::plotlyOutput,
    outputFunc2 = plotly::plotlyOutput,
    info.text = info_text,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height,
    label = label,
    title = "Dataset explorer",
    card_height = card_height
  )
}

loading_tsne_server <- function(id, pgx.dirRT,
                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot_data <- shiny::reactive({

      pgx.dir <- pgx.dirRT()
      tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")
      pgx.files <- sub("[.]pgx$", "", dir(pgx.dir, pattern = ".pgx$"))

      pos <- NULL
      if (file.exists(tsne.file)) {
        pos <- read.csv(tsne.file, row.names = 1)
        dim(pos)
        pos.files <- unique(gsub("^\\[|\\].*", "", rownames(pos)))
        if (!all(pgx.files %in% pos.files)) {
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
        pos <- Rtsne::Rtsne(1 - abs(corF),
          perplexity = px,
          check_duplicates = FALSE, is_distance = TRUE
        )$Y
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

      return(list(df, dpos))
    })

    plot.RENDER <- function() {
      df <- plot_data()
      shiny::req(df)
      
      dataset_pos <- df[[2]]
      marker_size <- ifelse( nrow(df[[1]]) > 50, 5, 10)
      
      fig <- plotly::plot_ly(
        data = df[[1]],
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

      dy <- diff(range(dataset_pos[,"y"]))
      dbg("[loading_tsneplot.R] range.y=",dy)
      
      fig <- fig %>%
        plotly::add_annotations(
          x = dataset_pos[,"x"],
          y = dataset_pos[,"y"],
          text = rownames(dataset_pos),
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
            showticklabels = FALSE
          ),
          yaxis = list(
            title = "tsne-y",
            showticklabels = FALSE
          )
        )

      fig
    }

    modal_plot.RENDER <- function() {
      df <- plot_data()
      shiny::req(df)
      marker_size <- ifelse( nrow(df[[1]]) > 50, 6, 12)
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
