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
  info.references
) {
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
    height = height,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "grouped_barplot",
    palette_default = "original"
  )
}

dataview_plot_tissue_server <- function(id, pgx, r.gene, r.data_type, watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Override bar ordering default to ascending
    shiny::observeEvent(TRUE, once = TRUE, {
      shiny::updateSelectInput(
        session, "bars_order",
        selected = "ascending"
      )
    })

    ## Original plot_data reactive — kept intact
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
      } else if (homologue_ratio > .5) {
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

    ## Editor: dynamic color pickers for custom palette
    output$custom_palette_ui <- shiny::renderUI({
      shiny::req(input$palette == "custom")
      pdat <- plot_data()
      shiny::req(pdat)
      grp_indices <- sort(unique(pdat$df$group))
      shiny::req(length(grp_indices) > 0)
      default_pal <- omics_pal_d()(8)
      pickers <- lapply(seq_along(grp_indices), function(i) {
        colourpicker::colourInput(
          ns(paste0("custom_color_", grp_indices[i])),
          label = paste("Group", grp_indices[i]),
          value = default_pal[grp_indices[i]]
        )
      })
      shiny::tagList(pickers)
    })

    ## Editor: rank list for custom drag-and-drop ordering
    output$rank_list <- shiny::renderUI({
      pdat <- plot_data()
      shiny::req(pdat)
      df <- pdat$df

      ## Apply current ordering so drag-and-drop starts from it
      bar_order <- if (!is.null(input$bars_order)) input$bars_order else "ascending"
      if (bar_order == "ascending" || bar_order == "custom") {
        df <- df[order(df$x), ]
      } else if (bar_order == "alphabetical") {
        df <- df[order(df$tissue), ]
      }
      ## "descending" keeps the default order

      labels <- as.character(df$tissue)
      sortable::bucket_list(
        header = NULL,
        class = "default-sortable custom-sortable",
        sortable::add_rank_list(
          input_id = ns("rank_list_order"),
          text = NULL,
          labels = labels
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

      ## Editor: bar ordering (default ascending)
      bar_order <- if (!is.null(input$bars_order)) input$bars_order else "ascending"
      if (bar_order == "ascending") {
        df <- df[order(df$x), ]
      } else if (bar_order == "alphabetical") {
        df <- df[order(df$tissue), ]
      } else if (bar_order == "custom" && !is.null(input$rank_list_order)) {
        custom_order <- input$rank_list_order
        custom_order <- intersect(custom_order, df$tissue)
        if (length(custom_order) > 0) {
          df <- df[match(custom_order, df$tissue), ]
        }
      }
      ## "descending" is already the default from plot_data

      ## Set tissue factor levels to preserve ordering
      df$tissue <- stringr::str_to_title(paste(df$tissue, " "))
      df$tissue <- factor(df$tissue, levels = df$tissue)

      ## Editor: palette override
      palette <- input$palette
      if (!is.null(palette) && palette == "custom") {
        ## Custom per-group colors
        grp_indices <- sort(unique(df$group))
        default_pal <- omics_pal_d()(8)
        bar_colors <- sapply(df$group, function(g) {
          val <- input[[paste0("custom_color_", g)]]
          if (is.null(val)) default_pal[g] else val
        })
        p <- plotly::plot_ly(
          data = df,
          y = ~tissue,
          x = ~x,
          type = "bar",
          orientation = "h",
          marker = list(color = bar_colors),
          hovertemplate = "%{y}: %{x}<extra></extra>"
        )
      } else if (!is.null(palette) && !palette %in% c("original", "")) {
        ## Named palette
        bar_colors <- omics_pal_d(palette = palette)(8)[df$group]
        p <- plotly::plot_ly(
          data = df,
          y = ~tissue,
          x = ~x,
          type = "bar",
          orientation = "h",
          marker = list(color = bar_colors),
          hovertemplate = "%{y}: %{x}<extra></extra>"
        )
      } else {
        ## Original behavior
        p <- plotly::plot_ly(
          data = df,
          y = ~tissue,
          x = ~x,
          type = "bar",
          orientation = "h",
          color = ~color,
          colors = omics_pal_d()(length(unique(df$color))),
          hovertemplate = "%{y}: %{x}<extra></extra>"
        )
      }

      p %>%
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
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
