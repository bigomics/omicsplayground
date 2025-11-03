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
  caption = "Boxplot of feature expression across time points. Within each boxplot, each dot is a sample. The colors distinguish the groups in the selected contrast. A spline interpolation (smooth curve) of expression across time points is performed and added within each group.",
  info.text = "Boxplot of feature expression across time points. Within each boxplot, each dot is a sample. The colors distinguish the groups in the selected contrast. A spline interpolation (smooth curve) of expression across time points is performed and added within each group.",
  info.methods = "A spline interpolation (smooth curve) of median expression data across time points is performed using the stats::spline R function with default parameters.",
  info.references = list(
    list("splinefun:", "https://stat.ethz.ch/R-manual/R-devel/library/stats/html/splinefun.html")
  ),
  info.extra_link = "extra.link",
  height = c("calc(100vh - 310px)", TABLE_HEIGHT_MODAL),
  width = c("auto", "100%")
) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::checkboxInput(
      ns("show_others"),
      "Show other groups",
      FALSE
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    ## plotlib = "plotly",
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
  info.text = "Table reporting results of differential gene expression testing for the main effect and, if available, the interaction with time. P-value and q-value columns show the meta p and meta q value, respectively, corresponding to max p and max q among selected methods. Avg0 and Avg1 columns report the average feature expression in each group. Interaction with time is tested using a design formula containing natural cubic spline of the 'time' variable detected in the metadata (i.e. ~ phenotype * spline(time)).",
  caption = "Table reporting results of differential gene expression testing for the main effect and, if available, the interaction with time. P-value and q-value correspond to max p and max q among selected methods. Avg0 and Avg1 columns report the average feature expression in each group. Interaction with time is tested using a design formula containing natural cubic spline of the 'time' variable detected in the metadata (i.e. ~ phenotype * spline(time)).",
  height = c("40%", TABLE_HEIGHT_MODAL),
  width = c("auto", "100%")
) {
  ns <- shiny::NS(id)

  options <- tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("show_statdetails"),
        "Show detailed statistical methods"
      ),
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
                                            data,
                                            timevar,
                                            contrast,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      sel.timevar <- timevar()
      genes <- rownames(pgx$X)
      genes <- table_module$rownames_all()
      genes <- head(genes, 16)

      expr <- pgx$X
      if (any(is.na(expr))) expr <- playbase::imputeMissing(expr, method = "SVD2")
      expr <- expr[genes, , drop = FALSE]
      time <- pgx$samples[, sel.timevar]

      ct <- contrast()
      group <- pgx$contrasts[, ct]
      group[is.na(group)] <- "others"

      if (!input$show_others) {
        kk <- which(!is.na(pgx$contrasts[, ct]))
        expr <- expr[, kk, drop = FALSE]
        time <- time[kk]
        group <- group[kk]
      }

      ngenes <- length(genes)
      xgenes <- as.vector(mapply(rep, genes, ncol(expr)))
      xexpr <- as.vector(t(expr))
      xtime <- rep(time, ngenes)
      xgroup <- rep(group, ngenes)

      ## long format
      df <- data.frame(gene = xgenes, expr = xexpr, time = xtime, group = xgroup)
      df
    })

    stats_data <- shiny::reactive({
      k <- contrast()
      shiny::req(k)

      kstats.full <- as.matrix(pgx$gx.meta$meta[[k]])
      sel <- grep("^p[.]*", colnames(kstats.full))
      kstats.full[, "meta.p"] <- apply(kstats.full[, sel, drop = FALSE], 1, function(x) x[which.max(x)])
      sel <- grep("^q[.]*", colnames(kstats.full))
      kstats.full[, "meta.q"] <- apply(kstats.full[, sel, drop = FALSE], 1, function(x) x[which.max(x)])
      cols <- c("meta.fx", "meta.p", "meta.q", "avg.0", "avg.1")
      cols <- intersect(cols, colnames(kstats.full))
      kstats <- kstats.full[, cols]
      col_mapping <- c("meta.fx" = "log2FC", "meta.p" = "p.value", "meta.q" = "q.value", "avg.0" = "avg.0", "avg.1" = "avg.1")
      new_names <- col_mapping[cols]
      colnames(kstats) <- new_names

      ik <- paste0("IA:", k)
      if (ik %in% names(pgx$gx.meta$meta)) {
        ikstats.full <- as.matrix(pgx$gx.meta$meta[[ik]])
        sel <- grep("^p[.]*", colnames(ikstats.full))
        ikstats.full[, "meta.p"] <- apply(ikstats.full[, sel, drop = FALSE], 1, function(x) x[which.max(x)])
        sel <- grep("^q[.]*", colnames(ikstats.full))
        ikstats.full[, "meta.q"] <- apply(ikstats.full[, sel, drop = FALSE], 1, function(x) x[which.max(x)])
        ikstats <- ikstats.full[, c("meta.p", "meta.q")]
        colnames(ikstats) <- c("p.interaction", "q.interaction")
        cm <- intersect(rownames(kstats), rownames(ikstats))
        kstats <- cbind(kstats[cm, ], ikstats[cm, ])
        cols <- c(
          "log2FC", "p.value", "q.value", "p.interaction",
          "q.interaction", "avg.0", "avg.1"
        )
        cols <- intersect(cols, colnames(kstats))
        if (length(cols) > 0) kstats <- kstats[, cols, drop = FALSE]
      }

      if (input$show_statdetails) {
        sel.p <- grep("^p[.]", colnames(kstats.full))
        sel.q <- grep("^q[.]", colnames(kstats.full))
        pq.tables <- kstats.full[, c(sel.p, sel.q), drop = FALSE]
        if (ik %in% names(pgx$gx.meta$meta)) {
          sel.p <- grep("^p[.]", colnames(ikstats.full))
          sel.q <- grep("^q[.]", colnames(ikstats.full))
          ik.pq.tables <- ikstats.full[, c(sel.p, sel.q), drop = FALSE]
          colnames(ik.pq.tables) <- paste0(colnames(ik.pq.tables), ".interaction")
          pq.tables <- cbind(pq.tables, ik.pq.tables)
        }
        kstats <- cbind(kstats, pq.tables)
      }

      kstats <- as.data.frame(kstats, check.names = FALSE)
      kstats <- kstats[, unique(colnames(kstats))]
      kstats <- kstats[order(-kstats$log2FC), ]

      return(kstats)
    })

    ## -----------------------------------------------------
    ## ----------------------- Plot -----------------------
    ## -----------------------------------------------------

    render_plot <- function() {
      library(ggplot2)
      library(plotly)
      df <- plot_data()
      shiny::req(df)

      timevar <- timevar()

      par(mfrow = c(3, 3), mar = c(5, 4, 2, 1))
      genes <- head(unique(df$gene), 9)
      for (g in genes) {
        ii <- which(df$gene == g)
        tt <- df$time[ii]
        xx <- df$expr[ii]
        gr <- df$group[ii]
        g <- playbase::probe2symbol(g, pgx$genes, "gene_name", fill_na = TRUE)
        playbase::plotTimeSeries.groups(
          time = tt, y = xx, group = gr, main = g, lwd = 3,
          xlab = timevar, time.factor = TRUE
        )
      }
    }

    PlotModuleServer(
      "plot",
      func = render_plot,
      plotlib = "base",
      ## csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 110), ## resolution of plots
      pdf.width = 14,
      pdf.height = 3.5,
      add.watermark = watermark
    )

    ## -----------------------------------------------------
    ## ----------------------- Table -----------------------
    ## -----------------------------------------------------

    render_table <- function(full = FALSE) {
      df <- stats_data()
      shiny::req(df)
      ft <- gsub("[;].*", ";...", rownames(df))
      df <- as.data.frame(df, check.names = FALSE)

      ## add module information
      module <- data()$modules
      module <- module[match(rownames(df), names(module))]
      module[is.na(module)] <- "-"

      ## do not show symbol column if symbol==feature
      symbol <- pgx$genes[rownames(df), "symbol"]
      na.symbol <- is.na(symbol) | symbol == ""
      if (!full && mean(!na.symbol) > 0.66) {
        df1 <- cbind(module = module, symbol = symbol, df)
      } else if (!full) {
        df1 <- cbind(module = module, feature = ft, df)
      } else {
        df1 <- cbind(module = module, feature = ft, symbol = symbol, df)
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

    render_table2 <- function() {
      render_table(full = TRUE)
    }

    table_module <- TableModuleServer(
      "table",
      func = render_table,
      func2 = render_table2,
      selector = "none"
    )
  }) ## end of moduleServer
}
