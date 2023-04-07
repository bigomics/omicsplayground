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
expression_table_genetable_ui <- function(
  id,
  title,
  info.text,
  caption,
  width,
  height) {
  ns <- shiny::NS(id)

  genetable_opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("gx_top10"), "top 10 up/down genes", FALSE),
      "Display only top 10 differentially (positively and negatively) expressed genes in the table.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(shiny::checkboxInput(ns("gx_showqvalues"), "show indivivual q-values", FALSE),
      "Show q-values of each indivivual statistical method in the table.",
      placement = "top", options = list(container = "body")
    )
  )

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    width = width,
    caption = caption,
    height = height,
    options = genetable_opts,
    title = title,
    label = "I"
  )
}

#' Server side table code: expression board
#'
#' @param id
#' @param watermark
#'
#' @export
expression_table_genetable_server <- function(id,
                                              res, # filteredDiffExprTable
                                              height,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    table.RENDER <- function() {
      res <- res()

      if (is.null(res) || nrow(res) == 0) {
        return(NULL)
      }

      fx.col <- grep("fc|fx|mean.diff|logfc|foldchange", tolower(colnames(res)))[1]
      fx.col
      fx <- res[, fx.col]

      if ("gene_title" %in% colnames(res)) res$gene_title <- playbase::shortstring(res$gene_title, 50)
      rownames(res) <- sub(".*:", "", rownames(res))

      if (!DEV) {
        kk <- grep("meta.fx|meta.fc|meta.p", colnames(res), invert = TRUE)
        res <- res[, kk, drop = FALSE]
      }
      if (!input$gx_showqvalues) {
        kk <- grep("^q[.]", colnames(res), invert = TRUE)
        res <- res[, kk, drop = FALSE]
      }

      numeric.cols <- which(sapply(res, is.numeric))
      numeric.cols <- colnames(res)[numeric.cols]

      DT::datatable(res,
        rownames = FALSE,
        ## class = 'compact cell-border stripe hover',
        class = "compact hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          # paging = TRUE,
          # pageLength = 16, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE,
          scrollY = 240,
          scroller = TRUE,
          deferRender = TRUE,
          search = list(
            regex = TRUE,
            caseInsensitive = TRUE
            ## , search = 'M[ae]'
          )
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(colnames(res)[fx.col],
          ## background = DT::styleColorBar(c(0,3), 'lightblue'),
          background = playbase::color_from_middle(fx, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }

    table.RENDER_modal <- function() {
      dt <- table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    }

    genetable <- TableModuleServer(
      "datasets",
      func = table.RENDER,
      func2 = table.RENDER_modal,
      selector = "single"
    )

    return(genetable)
  })
}
