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
correlation_plot_barplot_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::selectInput(ns("order_opt"), "Order by:",
        choices = c(
          "both",
          "correlation",
          "partial Correlation"
        ),
        multiple = FALSE,
        selected = "both"
      ),
      "Sort order of groups based on correlation.",
      placement = "top"
    )
  )

  PlotModuleUI(
    id = ns("plot"),
    title = title,
    label = label,
    plotlib = "plotly",
    caption = caption,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
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
correlation_plot_barplot_server <- function(id,
                                            getPartialCorrelation,
                                            getGeneCorr,
                                            cor_table,
                                            watermark = FALSE,
                                            pgx,
                                            labeltype) {
  moduleServer(id, function(input, output, session) {
    # reactive function listeninng for changes in input
    plot_data <- shiny::reactive({
      df <- getPartialCorrelation()
      R <- getGeneCorr()

      sel <- cor_table$rownames_current()
      shiny::req(sel)

      ## cor_table!=R mismatch!!!
      if (length(sel) > nrow(R)) {
        return(NULL)
      }

      NTOP <- 40 ## how many genes to show in barplot
      sel <- intersect(sel, rownames(R))
      sel <- head(sel, NTOP)
      rho <- R[sel, "cor"]
      if (length(sel) == 1) names(rho) <- rownames(R)[sel]

      prho <- df$pcor
      names(prho) <- playbase::probe2symbol(rownames(df), pgx$genes, labeltype(), fill_na = TRUE)
      names(rho) <- playbase::probe2symbol(names(rho), pgx$genes, labeltype(), fill_na = TRUE)
      prho <- prho[match(names(rho), names(prho))]

      pd <- data.frame(
        "correlation" = rho,
        "partial correlation" = prho
      )

      return(pd)
    })

    render_barplot <- function() {
      pd <- plot_data()

      playbase::pgx.stackedBarplot(
        x = pd,
        ylab = "Correlation",
        xlab = "",
        showlegend = FALSE
      )
    }

    barplot.RENDER <- function() {
      render_barplot() %>%
        plotly_default()
    }

    barplot.RENDER2 <- function() {
      render_barplot() %>%
        plotly_modal_default()
    }

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = barplot.RENDER,
      func2 = barplot.RENDER2,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(63, 100), ## resolution of plots
      pdf.width = 6, pdf.height = 4,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
