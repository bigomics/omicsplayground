##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


dataview_plot_tissue_ui <- function(
    id,
    label = "",
    title,
    height,
    width,
    caption,
    info.text,
    info.methods,
    info.references) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    caption = caption,
    options = NULL,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

dataview_plot_tissue_server <- function(id, pgx, r.gene, r.data_type, watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(r.gene(), r.data_type())

      ## dereference reactive
      gene <- r.gene()
      data_type <- r.data_type()

      # Find ortholog proportion
      n <- length(pgx$genes$human_ortholog)
      ortho <- sum(pgx$genes$human_ortholog != "")

      homologue_ratio <- ortho / n
      if (pgx$organism %in% c("Human", "human")) {
        hgnc.gene <- pgx$genes[gene, "symbol"]
      } else if (homologue_ratio > .9) {
        hgnc.gene <- pgx$genes[gene, "human_ortholog"]
      } else {
        shiny::validate(shiny::need(FALSE, "No tissue data available for this organism."))
      }

      tx <- tissue.klr <- grp <- NULL
      is.summary.available <- hgnc.gene %in% rownames(playdata::TISSUE)

      shiny::validate(shiny::need(is.summary.available, tspan("No tissue data available for this gene.", js = FALSE)))

      if (is.summary.available) {
        tx <- playdata::TISSUE[hgnc.gene, ]
        grp <- playdata::TISSUE_GRP[names(tx)]
        tissue.klr <- omics_pal_d()(8)[grp]
        ylab <- "expression (TPM)"
        if (data_type == "logCPM") {
          ylab <- "expression (log2TPM)"
          tx <- log(1 + tx)
        }
        jj <- 1:length(tx)
        sorting <- "no"
        if (sorting == "decr") jj <- order(-tx)
        if (sorting == "inc") jj <- order(tx)

        tx <- tx[jj]
        tissue.klr <- tissue.klr[jj]
      }

      df <- data.frame(
        tissue = names(tx),
        x = tx,
        group = grp,
        color = tissue.klr
      )
      df <- df[with(df, order(-x)), ]
      df <- head(df, 15) # select top 15 tissues

      return(
        list(
          df = df,
          gene = hgnc.gene,
          ylab = ylab
        )
      )
    })

    plot.RENDER <- function() {
      pdat <- plot_data()
      shiny::req(pdat)

      df <- pdat$df
      ylab <- stringr::str_to_sentence(pdat$ylab)
      if (pgx$organism %in% c("Human", "human")) {
        title <- FALSE
      } else {
        title <- "Expression in human tissue"
      }

      ## plot as regular bar plot
      df <- dplyr::mutate(
        df,
        tissue = forcats::fct_reorder(
          stringr::str_to_title(paste(tissue, " ")), x
        )
      )


      plotly::plot_ly(
        data = df,
        y = ~tissue,
        x = ~x,
        type = "bar",
        orientation = "h",
        color = ~color, ## TODO: use variable that encodes grouping
        colors = omics_pal_d()(length(unique(df$color))),
        hovertemplate = "%{y}: %{x}<extra></extra>"
      ) %>%
        plotly::layout(
          yaxis = list(title = title),
          xaxis = list(title = ylab),
          font = list(family = "Lato"),
          showlegend = FALSE,
          bargap = .4,
          margin = list(l = 10, r = 10, b = 10, t = 10)
        ) %>%
        plotly_default()
    }

    modal_plot.RENDER <- function() {
      plot.RENDER() %>%
        plotly_modal_default() %>%
        plotly::layout(
          showlegend = FALSE ## TODO: I guess a legend makes sense here?
        )
    }

    plot_data_csv <- function() {
      pdat <- plot_data()
      df <- pdat$df
      df$group <- NULL
      df$color <- NULL
      return(df)
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      csvFunc = plot_data_csv, ##  *** downloadable data as CSV
      res = c(90, 170), ## resolution of plots
      pdf.width = 8, pdf.height = 4,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
