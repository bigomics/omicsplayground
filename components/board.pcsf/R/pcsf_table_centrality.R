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

  table_opts <- shiny::tagList()

  TableModuleUI(
    ns("table"),
    info.text = info.text,
    caption = caption,
    width = width,
    height = height,
    ##    options = table_opts,
    title = title,
  )
}

#' Server side table code: expression board
#'
#' @param id
#' @param watermark
#'
#' @export
pcsf_table_centrality_server <- function(id,
                                         pgx,
                                         r_contrast,
                                         r_pcsf) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    table.RENDER <- function() {
      df <- playbase::pgx.getPCSFcentrality(
        pgx,
        contrast = r_contrast(),
        pcsf = r_pcsf(),
        n = 100,
        plot = FALSE
      )

      num.cols <- match(c("centrality", "logFC"), colnames(df))

      dt <- DT::datatable(df,
        rownames = FALSE,
        extensions = c("Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = c(1)),
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          pageLength = 9999,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE,
          scrollResize = TRUE,
          scrollY = 100,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatSignif(columns = num.cols, digits = 4) %>%
        DT::formatStyle(
          "logFC",
          background = color_from_middle(df$logFC, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) %>%
        DT::formatStyle(
          "centrality",
          background = color_from_middle(df$centrality, "white", "#fec34d"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )

      return(dt)
    }

    table.RENDER_modal <- function() {
      dt <- table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      return(dt)
    }

    TableModuleServer(
      "table",
      func = table.RENDER,
      func2 = table.RENDER_modal,
      selector = "none"
    )
  }) # end module server
} # end server
