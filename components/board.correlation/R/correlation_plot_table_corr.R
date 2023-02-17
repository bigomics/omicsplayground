##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics Sagl. All rights reserved.
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
correlation_plot_table_corr_ui <- function(id,
                                           label = "",
                                           height,
                                           width) {
  ns <- shiny::NS(id)
  info_text <- "<b>Top correlated genes.</b> Highest correlated genes in respect to the selected gene. The height of the bars correspond to the Pearson correlation value. The dark grey bars correspond to the 'partial correlation' which essentially corrects the correlation value for indirect effects and tries to estimate the amount of direct interaction."

  cor_table.info <- "<b>DGCA table.</b> Statistical results from the DGCA computation for differentially correlated gene pairs."

  div(
    PlotModuleUI(ns("plot"),
      title = "Top correlated genes",
      label = label,
      plotlib = "base",
      info.text = info_text,
      download.fmt = c("png", "pdf", "csv"),
      width = width,
      height = height
    ),
    TableModuleUI(
      ns("datasets"),
      info.text = cor_table.info,
      height = c(360, 700),
      width = c("auto", "90%"),
      title = "Correlation table",
      label = "b"
    )
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
                                               inputData,
                                               pcor_ntop,
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
        rho, prho
      ))
    })

    cor_barplot.PLOTFUN <- function() {
      df <- plot_data()
      rho <- df[[1]]
      prho <- df[[2]]

      ylim0 <- c(-1, 1) * max(abs(rho)) * 1.05

      par(mfrow = c(1, 1), mar = c(10, 4, 1, 0.5))
      barplot(rho,
        beside = FALSE, las = 3,
        ylim = ylim0,
        ylab = "correlation",
        cex.names = 0.85
      )
      barplot(prho,
        beside = FALSE, add = TRUE,
        col = "grey40", names.arg = ""
      )
      legend("topright",
        cex = 0.85, y.intersp = 0.85,
        inset = c(0.035, 0),
        c("correlation", "partial correlation"),
        fill = c("grey70", "grey40")
      )
    }

    PlotModuleServer(
      "plot",
      plotlib = "base",
      func = cor_barplot.PLOTFUN,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(63, 100), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )

    ### TABLE

    cor_table.RENDER <- shiny::reactive({
      ngs <- inputData()
      shiny::req(ngs)

      R <- getGeneCorr()
      if (is.null(R)) {
        return(NULL)
      }

      P <- getPartialCorrelation()
      pcor <- P[match(rownames(R), rownames(P)), "pcor"]

      title <- ngs$genes[rownames(R), "gene_title"]
      title <- substring(title, 1, 80)
      df <- data.frame(gene = rownames(R), title = title, cor = R[, "cor"], pcor = pcor)

      numeric.cols <- colnames(df)[3:ncol(df)]
      ## selectmode <- ifelse(input$corGSEAtable_multiselect,'multiple','single')

      DT::datatable(
        df,
        rownames = FALSE, ## escape = c(-1),
        extensions = c("Buttons", "Scroller"),
        ## selection=list(mode='multiple', target='row', selected=c(1)),
        selection = list(mode = "single", target = "row", selected = c(1)),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE, ## scrollY = TRUE,
          ## scrollY = 170,
          scrollY = "70vh",
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("cor",
          background = color_from_middle(
            df[, "cor"], "lightblue", "#f5aeae"
          ),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    TableModuleServer(
      "datasets",
      func = cor_table.RENDER,
      selector = "none"
    )
  }) ## end of moduleServer
}
