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
      validate(need(nrow(info.table)>0, 'Need at least one dataset!'))      
            
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

      allfc.file <- file.path(pgx.dir, "datasets-allFC.csv")
      if(!file.exists(allfc.file)) {
        return(NULL)
      }

      ## if no t-SNE file exists, we need to calculate it
      if (is.null(pos) && file.exists(allfc.file)) {

        shiny::withProgress(
          message = "Calculating signature similarities...", value = 0.33, {

            F <- data.table::fread(allfc.file)
            F <- as.matrix(F[, -1], rownames = F[[1]])
            fnames <- colnames(F)
            ## 1: Make this fast as possible!!! (IK)
            ## 2: Should we calculate few layouts/methods?
            F[is.na(F)] <- 0 ## really??

            ## get top 2000
            sel <- head(order(-rowMeans(F**2)),2000)
            F <- F[sel,]
            
            if(NCOL(F)==1) {
              pos <- matrix(0, 1, 2)
              rownames(pos) <- colnames(F)
              colnames(pos) <- c("x","y")
            } else {            
              ##F <- apply(F, 2, rank, na.last = "keep")
              rmsF <- (sqrt(colSums(F**2,na.rm=TRUE)) + 1e-8)
              F <- F %*% Matrix::Diagonal(x=1/rmsF)  ## fast scale
              F <- as.matrix(F)
              F[is.na(F)] <- 0 ## really??
              colnames(F) <- fnames  ## might be lost...
              
              if(0) {
                system.time( corF <- Rfast::Crossprod(F,F)) ## fast innerprod
                corF <- abs(corF) ## negative corr is also good...
                ##corF <- corF / max(diag(corF),na.rm=TRUE)  ## normalize diag=1??
                corF[is.na(corF)] <- 0
                distF <- pmax(-log(corF + 1e-8),0)  ## transform to range [0,inf]                        
                ppx <- max(min(30, floor(ncol(corF) / 4)), 1)
                pos <- try( Rtsne::Rtsne( distF,
                                         perplexity = ppx,
                                         check_duplicates = FALSE,
                                         is_distance = TRUE
                                         )$Y )
              } else {
                ppx <- max(min(30, floor(ncol(F) / 4)), 1)                
                pos <- try( Rtsne::Rtsne( t(abs(F)),
                                         perplexity = ppx,
                                         check_duplicates = FALSE,
                                         is_distance = FALSE
                                         )$Y )
              }
              ## safe...
              if("try-error" %in% class(pos)) {
                pos <- svd(F)$v[,1:2]
              }
            } 
            colnames(pos) <- c("x","y")
            rownames(pos) <- colnames(F)
        })
        
        ##plot(pos)
        pos <- round(pos, digits = 4)
        colnames(pos) <- c("x", "y")        
        write.csv(pos, file = tsne.file)
      }

      ## filter out non-existing entries
      pos.pgx <- gsub("^\\[|\\].*", "", rownames(pos))
      pos <- pos[which(pos.pgx %in% pgx.files), , drop = FALSE]

      dset <- gsub("^\\[|\\].*", "", rownames(pos))
      comparison <- gsub("^.*\\]", "", rownames(pos))
      colnames(pos) <- c("x", "y")
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
      
      marker_size <- ifelse( nrow(df) > 50, 8, 11)
      marker_size <- ifelse( nrow(df) > 100, 5, marker_size)
      font_size <- marker_size**0.55 * 5
      
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
          font = list( size=font_size ),
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
      marker_size <- ifelse( nrow(df) > 50, 9, 13)
      marker_size <- ifelse( nrow(df) > 100, 6, marker_size)
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
      ## res = c(100,300)*1,              ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      ## label = label, title = "t-SNE clustering",
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
