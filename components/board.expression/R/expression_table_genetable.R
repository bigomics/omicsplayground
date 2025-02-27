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
    withTooltip(shiny::checkboxInput(ns("gx_top10"), tspan("top 10 up/down genes"), FALSE),
      "Display only top 10 differentially (positively and negatively) expressed genes in the table.",
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
                                              res,
                                              organism,
                                              show_pv,
                                              height,
                                              scrollY,
                                              cont,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    table_data <- function() {
      res <- res()
      req(res)
      
      if ("gene_title" %in% colnames(res)) {
        res$gene_title <- playbase::shortstring(res$gene_title, 50)
      }
      ##rownames(res) <- sub(".*:", "", rownames(res))

      if (show_pv()) {
        res <- res[, -grep(".q$", colnames(res)), drop = FALSE]
      } else {
        res <- res[, -grep(".p$", colnames(res)), drop = FALSE]
      }

      if (input$gx_top10) {
        res <- res[!is.na(res$logFC), ]
        res <- res[order(res$logFC, decreasing = TRUE), ]
        if (nrow(res) >= 20) {
          res <- rbind(res[1:10, ], res[(nrow(res) - 9):nrow(res), ])
        } else {
          pos <- any(res$logFC > 0)
          neg <- any(res$logFC < 0)
          if (pos && neg) {
            res <- rbind(res[res$logFC > 0, ], res[res$logFC < 0, ])
          }
        }
      }
      res
    }

    table.RENDER <- function(showdetails = FALSE) {
      df <- table_data()
      df$gene_name <- NULL

      if (organism %in% c("Human", "human")) {
        df$human_ortholog <- NULL
      }
      if (sum(df$feature %in% df$symbol) > nrow(df) * .8) {
        df$feature <- NULL
      }

      if (!showdetails) {
        hide.cols <- grep("^AveExpr|p$|q$", colnames(df))
        hide.cols <- setdiff(hide.cols, grep("^meta", colnames(df)))
        if (length(hide.cols)) df <- df[, -hide.cols]
      }

      numeric.cols <- which(sapply(df, is.numeric))
      numeric.cols <- colnames(df)[numeric.cols]
      fx.col <- grep("fc|fx|mean.diff|logfc|foldchange", tolower(colnames(df)))[1]
      fx <- df[, fx.col]

      DT::datatable(df,
        rownames = FALSE,
        class = "compact hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        plugins = "scrollResize",
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          scrollX = TRUE,
          scrollY = scrollY,
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE,
          search = list(
            regex = TRUE,
            caseInsensitive = TRUE
          )
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(colnames(df)[fx.col],
          background = color_from_middle(fx, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }

    table.RENDER_modal <- function() {
      dt <- table.RENDER(showdetails = TRUE)
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    }

    table_csv <- function() {
      dt <- table_data()
      dt$stars <- NULL
      return(dt)
    }

    genetable <- TableModuleServer(
      "datasets",
      func = table.RENDER,
      func2 = table.RENDER_modal,
      csvFunc = table_csv,
      selector = "single",
      download.contrast.name = cont
    )

    return(genetable)
  })
}
