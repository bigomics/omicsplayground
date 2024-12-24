##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

dataview_plot_correlation_ui <- function(
    id,
    label = "",
    title,
    height,
    width,
    caption,
    info.text,
    info.methods,
    info.extra_link) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("pltsrv"),
    title = title,
    label = label,
    plotlib = "plotly",
    caption = caption,
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    options = NULL,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
  )
}

dataview_plot_correlation_server <- function(id,
                                             pgx,
                                             r.gene = reactive(""),
                                             r.samples = reactive(NULL),
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    getTopCorrelatedGenes <- function(pgx, gene, n, samples) {
      ## precompute
      if (is.null(samples)) samples <- colnames(pgx$X)
      if (!all(samples %in% colnames(pgx$X))) {
        return(NULL)
      }
      if (!gene %in% rownames(pgx$X)) {
        return(NULL)
      }

      samples <- intersect(samples, colnames(pgx$X))
      probe <- gene

      ## corr always in log.scale and restricted to selected samples subset
      ## should match exactly the rawtable!!
      if (probe %in% rownames(pgx$X)) {
        rho <- cor(t(pgx$X[, samples]), pgx$X[probe, samples], use = "pairwise")[, 1]
      } else if (probe %in% rownames(pgx$counts)) {
        x0 <- playbase::logCPM(pgx$counts[, samples])
        x1 <- x0[probe, ]
        rho <- cor(t(x0), x1, use = "pairwise")[, 1]
      } else {
        rho <- rep(0, nrow(pgx$genes))
        names(rho) <- rownames(pgx$genes)
      }

      rho[is.na(rho)] <- 0
      jj <- head(order(-abs(rho)), n)
      jj <- c(head(order(-rho), n / 2), tail(order(-rho), n / 2))
      top.rho <- rho[jj]

      gx1 <- sqrt(rowSums(pgx$X[names(top.rho), samples, drop = FALSE]**2, na.rm = TRUE))
      gx1 <- (gx1 / max(gx1))
      klr1 <- omics_pal_c(palette = "brand_blue")(16)[1 + round(15 * gx1)]
      klr1[which(is.na(klr1))] <- unname(omics_colors("mid_grey"))

      names(top.rho) <- sub(".*:", "", names(top.rho))

      ## NOTE: currently some labels are pretty long; also the var names are cryptic
      ## TODO: check if names can be shortened and variable names can be formatted nicely
      getGeneAnnot <- function(ngs, genes) {
        ann <- ngs$genes[, grep("name|title|chr|map", colnames(ngs$genes), ignore.case = TRUE)]
        apply(ann[genes, ], 1, function(x) {
          paste(mapply(paste0, colnames(ann), ": <b>", x), collapse = "</b><br>")
        })
      }
      annot <- getGeneAnnot(pgx, names(top.rho))

      df <- data.frame(
        genes = names(top.rho),
        rho = top.rho,
        color = klr1,
        value = gx1,
        annot = paste0(annot, "<extra></extra>")
      )

      res <- list(df, gene)

      return(res)
    }

    plot_data <- shiny::reactive({
      shiny::req(pgx$X, pgx$Y)
      shiny::req(r.gene())
      samples <- r.samples()
      gene <- r.gene()
      pd <- getTopCorrelatedGenes(pgx, gene = gene, n = 40, samples = samples)
      pd
    })

    plotly_render <- function() {
      pd <- plot_data()
      shiny::req(pd)

      df <- pd[[1]]
      gg <- unique(df$genes)
      df <- df[match(gg, df$genes), , drop = FALSE]
      df$genes <- factor(df$genes, levels = df$genes)

      ay <- list(overlaying = "y", side = "right", title = "")

      ## plot as regular bar plot
      plotly::plot_ly(
        data = df,
        x = ~genes,
        y = ~rho,
        type = "bar",
        marker = list(
          color = ~color # ,
        ),
        hovertemplate = ~annot
      ) %>%
        plotly::add_trace(
          data = df,
          x = ~genes,
          y = ~rho,
          type = "bar",
          name = "",
          yaxis = "y2",
          mode = "lines+markers"
        ) %>%
        plotly::layout(
          xaxis = list(title = FALSE),
          yaxis2 = ay,
          showlegend = FALSE,
          bargap = .4,
          margin = list(l = 10, r = 30, b = 10, t = 10)
        )
    }

    plotly.RENDER <- function() {
      plotly_render() %>%
        plotly::layout(
          xaxis = list(tickfont = list(size = 10))
        ) %>%
        plotly_default()
    }

    modal_plotly.RENDER <- function() {
      plotly_render() %>%
        plotly::layout(
          xaxis = list(tickfont = list(size = 18))
        ) %>%
        plotly_modal_default()
    }

    plot_data_csv <- function() {
      pd <- plot_data()
      df <- pd[[1]]
      df$genes <- factor(df$genes, levels = df$genes)
      df$gene_title <- stringr::str_extract(df$annot, "(?<=gene_title: <b>)[^<]+")
      df$gene_name <- stringr::str_extract(df$annot, "(?<=gene_name: <b>)[^<]+")
      df$color <- NULL
      df$annot <- NULL
      return(df)
    }

    PlotModuleServer(
      "pltsrv",
      plotlib = "plotly",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data_csv, ##  *** downloadable data as CSV
      res = c(80, 170), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
