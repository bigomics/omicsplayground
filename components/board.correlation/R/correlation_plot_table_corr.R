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
correlation_plot_corr_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::selectInput(ns("order_opt"), "Order by:",
        choices = c(
          "Both",
          "Correlation",
          "Partial Correlation"
        ),
        multiple = FALSE,
        selected = "Both"
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
      options = plot_opts,
      info.text = info.text,
      download.fmt = c("png", "pdf", "csv"),
      width = width,
      height = height
  )
}

#' @export
correlation_table_corr_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::selectInput(ns("order_opt"), "Order by:",
        choices = c(
          "Both",
          "Correlation",
          "Partial Correlation"
        ),
        multiple = FALSE,
        selected = "Both"
      ),
      "Sort order of groups based on correlation.",
      placement = "top"
    )
  )

  TableModuleUI(
      ns("datasets"),
      info.text = info.text,
      height = height,
      width = width,
      title = title,
      caption = caption,
      label = label
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
correlation_plot_table_corr_server <- function(id,
                                               getPartialCorrelation,
                                               getGeneCorr,
                                               cor_table,
                                               pgx,
                                               pcor_ntop,
                                               scrollY,
                                               watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    # reactive function listeninng for changes in input
    plot_data <- shiny::reactive({
      df <- getPartialCorrelation()
      R <- getGeneCorr()
      sel <- 1:pcor_ntop
      shiny::req(sel)

      rho <- R[sel, "cor"]
      if (length(sel) == 1) names(rho) <- rownames(R)[sel]

      prho <- df$pcor
      names(prho) <- rownames(df)
      prho <- prho[match(names(rho), names(prho))]
      names(prho) <- names(rho)

      return(list(
        rho, prho,
        order_opt = input$order_opt
      ))
    })

    cor_barplot.PLOTFUN <- function() {
      pd <- plot_data()

      pd_rho <- data.frame(genes = names(pd[[1]]), rho = pd[[1]])
      pd_prho <- data.frame(genes = names(pd[[2]]), rho = pd[[2]])

      pd_plot <- base::merge(pd_rho, pd_prho, by = "genes")

      pd_plot <- pd_plot[complete.cases(pd_plot), ]

      rownames(pd_plot) <- pd_plot$genes

      pd_plot$genes <- NULL

      colnames(pd_plot) <- c("Correlation", "Partial correlation")

      if (input$order_opt == "Correlation") {
        pd_plot <- pd_plot[order(pd_plot$Correlation, decreasing = TRUE), ]
      } else if (input$order_opt == "Partial Correlation") {
        pd_plot <- pd_plot[order(pd_plot$`Partial correlation`, decreasing = TRUE), ]
      } else if (input$order_opt == "Both") {
        total_sum_cor <- rowSums(pd_plot)

        pd_plot <- pd_plot[order(total_sum_cor, decreasing = TRUE), ]
      }

      playbase::pgx.stackedBarplot(
        x = pd_plot,
        ylab = "Correlation",
        showlegend = FALSE
      )
    }

    ### TABLE

    cor_table.RENDER <- shiny::reactive({
      shiny::req(pgx)

      R <- getGeneCorr()
      if (is.null(R)) {
        return(NULL)
      }

      P <- getPartialCorrelation()
      pcor <- P[match(rownames(R), rownames(P)), "pcor"]

      title <- pgx$genes[rownames(R), "gene_title"]
      title <- substring(title, 1, 80)
      df <- data.frame(gene = rownames(R), title = title, cor = R[, "cor"], pcor = pcor)

      numeric.cols <- colnames(df)[3:ncol(df)]
      ## selectmode <- ifelse(input$corGSEAtable_multiselect,'multiple','single')

      DT::datatable(
        df,
        rownames = FALSE, ## escape = c(-1),
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = c(1)),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = FALSE,
          scrollY = scrollY,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("cor",
          background = playbase::color_from_middle(
            df[, "cor"], "lightblue", "#f5aeae"
          ),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    cor_table.RENDER_modal <- shiny::reactive({
      dt <- cor_table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "datasets",
      func = cor_table.RENDER,
      func2 = cor_table.RENDER_modal,
      selector = "none"
    )

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = cor_barplot.PLOTFUN,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(63, 100), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
