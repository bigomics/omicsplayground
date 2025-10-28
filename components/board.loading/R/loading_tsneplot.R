##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

loading_tsne_ui <- function(
  id,
  title,
  info.text,
  info.references = NULL,
  info.methods = NULL,
  info.extra_link = NULL,
  caption,
  label = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("pltmod"),
    outputFunc = plotly::plotlyOutput,
    outputFunc2 = plotly::plotlyOutput,
    info.text = info.text,
    info.references = info.references,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    caption = caption,
    height = height,
    label = label,
    title = title
  )
}

loading_tsne_server <- function(id, pgx.dirRT, info.table, r_selected,
                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      pgx.dir <- pgx.dirRT()
      info.table <- info.table()
      validate(need(nrow(info.table) > 0, "Need at least one dataset!"))

      tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")
      #
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

      allfc.file <- file.path(pgx.dir, "datasets-allFC.csv")
      if (!file.exists(allfc.file)) {
        return(NULL)
      }

      ## if no t-SNE file exists, we need to calculate it
      if (is.null(pos) && file.exists(allfc.file)) {
        shiny::withProgress(
          message = "Calculating signature t-SNE...",
          value = 0.33,
          {
            F <- data.table::fread(allfc.file)
            F <- as.matrix(F[, -1], rownames = F[[1]])
            fnames <- colnames(F)
            ## 1: Make this fast as possible!!! (IK)
            ## 2: Should we calculate few layouts/methods?
            F[is.na(F)] <- 0 ## really??

            ## get top 2000
            sel <- head(order(-rowMeans(F**2, na.rm = TRUE)), 2000)
            F <- F[sel, ]

            if (NCOL(F) == 1) {
              pos <- matrix(0, 1, 2)
              rownames(pos) <- colnames(F)
              colnames(pos) <- c("x", "y")
            } else {
              rmsF <- (sqrt(colSums(F**2, na.rm = TRUE)) + 1e-8)
              F <- F %*% Matrix::Diagonal(x = 1 / rmsF) ## fast scale
              F <- as.matrix(F)
              F[is.na(F)] <- 0 ## really??
              colnames(F) <- fnames ## might be lost...
              ppx <- max(min(30, floor(ncol(F) / 4)), 1)
              pos <- try(Rtsne::Rtsne(t(abs(F)),
                perplexity = ppx,
                check_duplicates = FALSE,
                is_distance = FALSE
              )$Y)
              ## safe...
              if ("try-error" %in% class(pos)) {
                pos <- svd(F)$v[, 1:2]
              }
            }
            colnames(pos) <- c("x", "y")
            rownames(pos) <- colnames(F)
          }
        )

        #
        pos <- round(pos, digits = 4)
        colnames(pos) <- c("x", "y")
        write.csv(pos, file = tsne.file)
      }

      ## filter out non-existing entries
      pos.pgx <- gsub("^\\[|\\].*", "", rownames(pos))
      pos <- pos[which(pos.pgx %in% pgx.files), , drop = FALSE]

      dataset <- gsub("^\\[|\\].*", "", rownames(pos))
      comparison <- gsub("^.*\\]", "", rownames(pos))
      colnames(pos) <- c("x", "y")
      df <- data.frame(pos, dataset = dataset, comparison = comparison)

      ## compute medioid of datasets
      dataset_pos <- apply(pos, 2, function(x) tapply(x, dataset, median, na.rm = TRUE))
      if (length(unique(dataset)) == 1) {
        dataset_pos <- matrix(dataset_pos, nrow = 1, ncol = 2)
        rownames(dataset_pos) <- dataset[1]
      }
      colnames(dataset_pos) <- c("x", "y")

      pdata <- list(
        df = df,
        dataset_pos = dataset_pos
      )

      return(pdata)
    })

    plot.RENDER <- function() {
      pdata <- plot_data()
      shiny::req(pdata)
      dataset_pos <- pdata[["dataset_pos"]]
      df <- pdata[["df"]]

      ## filter with datatable active rows
      active_datasets <- info.table()$dataset[r_selected()]
      df <- df[which(df$dataset %in% active_datasets), , drop = FALSE]
      sel <- which(rownames(dataset_pos) %in% active_datasets)
      dataset_pos <- dataset_pos[sel, , drop = FALSE]

      marker_size <- ifelse(nrow(df) > 60, 8, 11)
      marker_size <- ifelse(nrow(df) > 120, 5, marker_size)
      font_size <- marker_size**0.55 * 5

      fig <- plotly::plot_ly(
        type = "scatter",
        mode = "markers",
        data = df,
        x = ~x,
        y = ~y,
        text = ~ I(paste(
          ifelse(nrow(df), "<b>Dataset:</b>", "Whoops!"), dataset,
          ifelse(nrow(df), "<br><b>Comparison:</b>", ""), comparison
        )),
        color = ~dataset,
        marker = list(
          size = marker_size,
          line = list(
            color = omics_colors("super_dark_grey"),
            width = 1.0
          )
        ),
        hovertemplate = "%{text}<br><b>tsne-x:</b> %{x:.2f}<br><b>tsne-y:</b> %{y:.2f}<extra></extra>"
      )

      fig <- fig %>%
        plotly::add_annotations(
          x = dataset_pos[, "x"],
          y = dataset_pos[, "y"],
          text = rownames(dataset_pos),
          font = list(size = font_size),
          xref = "x",
          yref = "y",
          xanchor = "middle",
          yanchor = "bottom",
          yshift = 0,
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
      fig <- plotly_default(fig)
      fig
    }

    modal_plot.RENDER <- function() {
      pdata <- plot_data()
      shiny::req(pdata)
      df <- pdata[["df"]]
      marker_size <- ifelse(nrow(df) > 60, 9, 13)
      marker_size <- ifelse(nrow(df) > 120, 6, marker_size)
      p <- plot.RENDER() %>%
        plotly::layout(
          showlegend = TRUE,
          font = list(
            size = 16
          )
        )
      p <- plotly::style(p, marker.size = marker_size)
      p <- plotly_modal_default(p)
      p
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      plotlib2 = "plotly",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      ## res = c(100,300)*1,              ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      ## label = label, title = "t-SNE clustering",
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
