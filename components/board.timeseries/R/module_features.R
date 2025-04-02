##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
TimeSeriesBoard.features_plot <- function(
    id,
    label = "label",
    title = "title",
    caption = "caption",
    info.text = "info.text",
    info.methods = "info.methods",
    info.references = list(),
    info.extra_link = "extra.link",
    height = c("calc(100vh - 310px)", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")
    ) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::checkboxInput(ns("show_others"), "Show other groups",FALSE)
  )
  
  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    ##plotlib = "plotly",
    options = options,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

TimeSeriesBoard.features_table <- function(
    id,
    label = "label",
    title = "title",
    info.text = "info.text",
    caption = "caption",
    height = c("40%", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")
    ) {
  ns <- shiny::NS(id)


  options <- tagList(
    withTooltip(
      shiny::checkboxInput(ns("show_statdetails"), "Show detailed statistical methods"),
      title = "Show detailed statistical methods."
    )
  )
  
  TableModuleUI(
    ns("table"),
    info.text = info.text,
    options = options,
    height = height,
    caption = caption,
    width = width,
    title = title,
    label = "b"
  )
}


#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#' @param comp
#' @param pgx
#' @param res
#' @param ii
#' @param watermark
#'
#'
#'
#' @export
TimeSeriesBoard.features_server <- function(id,
                                            pgx,
                                            # data,
                                            timevar,
                                            contrast,
                                            watermark = FALSE) {

  moduleServer(id, function(input, output, session) {

    plot_data <- shiny::reactive({
      
      sel.timevar <- timevar()      
      genes <- rownames(pgx$X)
      genes <- table_module$rownames_all()
      genes <- head(genes, 16)

      expr  <- pgx$X[genes,,drop=FALSE]
      time <- pgx$samples[,sel.timevar]

      ct <- contrast()
      group <- pgx$contrasts[,ct]
      group[is.na(group)] <- "others"

      if(!input$show_others) {
        kk <- which(!is.na(pgx$contrasts[,ct]))      
        expr  <- expr[,kk,drop=FALSE]
        time <- time[kk]
        group <- group[kk]
      }

      ngenes <- length(genes)
      xgenes <- as.vector(mapply(rep, genes, ncol(expr)))
      xexpr  <- as.vector(t(expr))
      xtime  <- rep(time, ngenes)
      xgroup <- rep(group, ngenes)
      
      ## long format
      df <- data.frame(gene=xgenes, expr = xexpr, time=xtime, group=xgroup)      
      df
    })

    stats_data <- shiny::reactive({
      k <- contrast()
      shiny::req(k)
      kstats <- as.matrix(pgx$gx.meta$meta[[k]]) #[,1:5]
      cols <- c("meta.fx","meta.p","meta.q","avg.0","avg.1")
      stats <- kstats[,cols]
      colnames(stats) <- c("log2FC","p.value","q.value","avg.0","avg.1")
      if(input$show_statdetails) {
        i=1;
        pq.tables <- kstats[, grep("^p[.]|^q[.]",colnames(kstats)), drop = FALSE]
        stats <- cbind(stats, pq.tables)
      }
      stats <- as.data.frame(stats, check.names=FALSE)

      ##stats <- stats[order(-abs(stats$log2FC)),]
      stats <- stats[order(stats$p.value),]

      return(stats)
    })
    
    ##-----------------------------------------------------
    ##----------------------- Plot -----------------------
    ##-----------------------------------------------------
    
    render_plot <- function() {
      library(ggplot2)
      library(plotly)
      df <- plot_data()
      shiny::req(df)

      timevar <- timevar()
      
      par(mfrow=c(3,3), mar=c(5,4,2,1))
      genes <- head(unique(df$gene),9)
      for(g in genes) {
        ii <- which(df$gene == g)
        tt <- df$time[ii]
        xx <- df$expr[ii]
        gr <- df$group[ii]
        playbase::plotTimeSeries.groups(
          time=tt, y=xx, group=gr, main=g, lwd=3,
          xlab=timevar, time.factor=TRUE)
      }
      
    }

    PlotModuleServer(
      "plot",
      func = render_plot,
      plotlib = "base",
      ##csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 110), ## resolution of plots
      pdf.width = 14,
      pdf.height = 3.5,
      add.watermark = watermark
    )

    ##-----------------------------------------------------
    ##----------------------- Table -----------------------
    ##-----------------------------------------------------

    render_table <- function() {

      df <- stats_data()
      shiny::req(df)
      ft <- gsub("[;].*",";...",rownames(df))
      df <- as.data.frame(df, check.names=FALSE)

      ## do not show symbol column if symbol==feature
      symbol <- pgx$genes[rownames(df),"symbol"]
      if(mean(symbol == ft, na.rm=TRUE) > 0.9) {
        df1 <- cbind(feature=ft, df)
      } else {
        df1 <- cbind(feature=ft, symbol=symbol, df)        
      }
      numeric.cols <- colnames(df)
      DT::datatable(
        df1,
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = "23vh",
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    table_module <- TableModuleServer(
      "table",
      func = render_table,
      selector = "none"
    )
    
  }) ## end of moduleServer
}
