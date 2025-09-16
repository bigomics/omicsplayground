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
clustering_table_clustannot_ui <- function(
    id,
    title,
    info.text,
    caption,
    width,
    height) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    width = width,
    height = height,
    info.text = info.text,
    title = title,
    caption = caption,
    label = "b"
  )
}

#' Server side table code: expression board
#'
#' @param id
#' @param watermark
#'
#' @export
clustering_table_clustannot_server <- function(
    id,
    getClustAnnotCorrelation,
    xann_level,
    watermark,
    scrollY) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    table_data <- reactive({
      rho <- getClustAnnotCorrelation()
      rho.name <- playbase::shortstring(sub(".*:", "", rownames(rho)), 60)
      df <- data.frame(feature = rho.name, round(as.matrix(rho), digits = 3))
      rownames(df) <- make.unique(rownames(rho))
      return(df)
    })

    clustannot_table.RENDER <- shiny::reactive({
      df <- table_data()
      rho <- getClustAnnotCorrelation()
      xann_level <- xann_level()
      if (is.null(rho)) {
        return(NULL)
      }

      if (xann_level == "geneset") {
        feature_link <- playbase::wrapHyperLink(
          rep_len("<i class='fa-solid fa-arrow-up-right-from-square'></i>", nrow(df)),
          rownames(df)
        ) |> HandleNoLinkFound(
          NoLinkString = "<i class='fa-solid fa-arrow-up-right-from-square'></i>",
          SubstituteString = "<i class='fa-solid fa-arrow-up-right-from-square blank_icon'></i>"
        )
      } else {
        feature_link <- ""
      }

      df$feature <- paste(df$feature, "&nbsp;", feature_link)

      DT::datatable(
        df,
        rownames = FALSE,
        escape = c(-1, -2),
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = c(1)),
        class = "compact hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip", buttons = c("copy", "csv", "pdf"),
          scrollX = TRUE,
          scrollY = scrollY,
          scroller = TRUE,
          scrollResize = TRUE,
          deferRender = TRUE
        ) ## end of options
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    clustannot_table.RENDER_modal <- shiny::reactive({
      dt <- clustannot_table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "datasets",
      func = clustannot_table.RENDER,
      func2 = clustannot_table.RENDER_modal,
      csvFunc = table_data,
      selector = "none"
    )
  }) # end module server
} # end server
