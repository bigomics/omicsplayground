##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' UI code for table code: expression board
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
pcsf_table_centrality_ui <- function(
    id,
    title,
    info.text,
    caption,
    width,
    height) {
  ns <- shiny::NS(id)

  table_opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("showq"), "show q-values", FALSE),
      "Show q-values next to FC values.",
      placement = "right", options = list(container = "body")
    )
  )

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    caption = caption,
    width = width,
    height = height,
    options = table_opts,
    title = title,
  )
}

#' Server side table code: expression board
#'
#' @param id
#' @param watermark
#'
#' @export
pcsf_table_centrality__server <- function(id,
                                          pgx,
                                          res, # filteredDiffExprTable
                                          metaFC,
                                          metaQ,
                                          height,
                                          scrollY,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    table.RENDER <- shiny::reactive({
      res <- res()
      if (is.null(res) || nrow(res) == 0) {
        return(NULL)
      }
      df <- data.frame()
      ## try(pcsf <- playbase::pgx.computePCSF(pgx, contrast = ct, ntop = 250))
      
      ##df <- playbase::pgx.getPCSFcentrality(
      ##  pgx, contrast=ct, pcsf=pcsf, n=10, plot=FALSE))

      
      dt <- DT::datatable(df,
        rownames = FALSE,
        extensions = c("Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = c(1)),
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE,
          scrollResize = TRUE,
          scrollY = scrollY,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatSignif(columns = fc.cols, digits = 4) %>%
        DT::formatStyle(
          "rms.FC",
          background = color_from_middle(fc.rms, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) %>%
        DT::formatStyle(
          fc.cols,
          background = color_from_middle(F, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )

      if (length(qv.cols) > 0) {
        dt <- dt %>%
          DT::formatSignif(columns = qv.cols, digits = 4)
      }
      return(dt)
    })

    table.RENDER_modal <- shiny::reactive({
      dt <- fctable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      return(dt)
    })

    TableModuleServer(
      "centrality",
      func = table.RENDER,
      func2 = table.RENDER_modal,
      selector = "none"
    )
  }) # end module server
} # end server
