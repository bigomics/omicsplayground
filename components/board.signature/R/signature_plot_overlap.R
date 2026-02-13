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
#'
#' @export
signature_plot_overlap_ui <- function(
  id,
  title,
  info.text,
  info.methods,
  info.extra_link,
  caption,
  height,
  width
) {
  ns <- shiny::NS(id)

  overlapScorePlot.opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("overlapScorePlot_ntop"),
        "Number of features", c(60, 120, 250),
        inline = TRUE
      ),
      "Specify the number to top features to show.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("overlapScorePlot_shownames"),
        "Show feature names", TRUE
      ),
      "Select to show/hide the feature names in the plot.",
      placement = "top", options = list(container = "body")
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    label = "a",
    plotlib = "plotly",
    options = overlapScorePlot.opts,
    caption = caption,
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    download.fmt = c("png", "pdf", "svg"),
    height = height,
    width = width,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "grouped_barplot",
    palette_default = "original"
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
signature_plot_overlap_server <- function(id,
                                          getOverlapTable,
                                          overlapTable,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Editor: override bars_order choices for this grouped plot
    shiny::observe({
      shiny::updateSelectInput(session, "bars_order",
        choices = c(
          "By group, score desc" = "group_desc",
          "By group, score asc" = "group_asc",
          "Value (descending)" = "descending",
          "Value (ascending)" = "ascending"
        ),
        selected = "group_desc"
      )
    })

    ## Shared reactive: filtered data after ntop (used by render + custom_palette_ui)
    plot_data <- shiny::reactive({
      df <- getOverlapTable()
      sel <- overlapTable$rows_all()
      shiny::req(df, sel)
      df1 <- df[sel, ]
      df1$geneset <- as.character(rownames(df1))
      ntop <- as.integer(input$overlapScorePlot_ntop)
      df1 <- df1[head(order(-df1$score), ntop), ]
      ## Re-factor db so only displayed groups remain
      df1$db <- factor(df1$db)
      df1
    })

    ## Editor: dynamic color pickers for custom palette
    output$custom_palette_ui <- shiny::renderUI({
      shiny::req(input$palette == "custom")
      df1 <- plot_data()
      shiny::req(df1)
      db_levels <- levels(df1$db)
      shiny::req(length(db_levels) > 0)
      default_clrs <- rep(RColorBrewer::brewer.pal(8, "Set2"), 10)[seq_along(db_levels)]
      pickers <- lapply(seq_along(db_levels), function(i) {
        colourpicker::colourInput(
          ns(paste0("custom_color_", i)),
          label = db_levels[i],
          value = default_clrs[i]
        )
      })
      shiny::tagList(pickers)
    })

    overlapScorePlot.RENDER <- shiny::reactive({
      df1 <- plot_data()
      shiny::req(df1)
      ntop <- as.integer(input$overlapScorePlot_ntop)

      ## Editor: bar ordering
      bars_order <- input$bars_order
      if (!is.null(bars_order)) {
        if (bars_order == "group_asc") {
          jj <- order(df1$db, df1$score)
        } else if (bars_order == "descending") {
          jj <- order(-df1$score)
        } else if (bars_order == "ascending") {
          jj <- order(df1$score)
        } else {
          ## default: group_desc
          jj <- order(df1$db, -df1$score)
        }
      } else {
        jj <- order(df1$db, -df1$score)
      }
      df1 <- df1[jj, ]

      df1$idx <- factor(1:nrow(df1), levels = 1:nrow(df1))
      df1$idx <- as.integer(df1$idx)

      ## Editor: palette override
      db_levels <- levels(df1$db)
      n_groups <- length(db_levels)
      palette <- input$palette
      if (!is.null(palette) && palette == "custom") {
        grp_clrs <- sapply(seq_len(n_groups), function(j) {
          val <- input[[paste0("custom_color_", j)]]
          if (is.null(val)) rep(RColorBrewer::brewer.pal(8, "Set2"), 10)[j] else val
        })
      } else if (!is.null(palette) && !palette %in% c("original", "")) {
        grp_clrs <- rep(omics_pal_d(palette = palette)(8), ceiling(n_groups / 8))[1:n_groups]
      } else {
        grp_clrs <- rep(RColorBrewer::brewer.pal(8, "Set2"), 10)[1:n_groups]
      }
      klr <- grp_clrs[as.integer(df1$db)]

      plt <- plotly::plot_ly(
        df1,
        x = ~idx, y = ~score,
        type = "bar",
        hovertemplate = paste0(df1$geneset, "<br>score: %{y}<extra></extra>"),
        marker = list(color = klr)
      ) %>%
        plotly::layout(
          showlegend = FALSE,
          dragmode = "select",
          yaxis = list(
            titlefont = list(size = 11),
            tickfont = list(size = 10),
            showgrid = TRUE,
            title = "overlap score"
          ),
          xaxis = list(
            title = "",
            showgrid = FALSE,
            showline = FALSE,
            showticklabels = FALSE,
            showgrid = FALSE,
            zeroline = FALSE
          )
        )

      if (min(nrow(df1), ntop) < 100 && input$overlapScorePlot_shownames) {
        ## labeling the y-axis inside bars
        plt <- plt %>%
          plotly::add_annotations(
            yref = "paper", xref = "x",
            x = ~idx, y = 0.005, yanchor = "bottom",
            text = substring(df1$geneset, 1, 35),
            textangle = -90,
            font = list(size = 10),
            showarrow = FALSE, align = "right"
          )
      }
      plt
    })

    PlotModuleServer(
      "plot",
      func = overlapScorePlot.RENDER,
      plotlib = "plotly",
      res = c(100, 100), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
