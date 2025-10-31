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
expression_table_FDRtable_ui <- function(
  id,
  title,
  caption,
  info.text,
  width,
  height
) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    caption = caption,
    width = width,
    height = height,
    title = title
  )
}

#' Server side table code: expression board
#'
#' @param id
#' @param watermark
#'
#' @export
expression_table_FDRtable_server <- function(id,
                                             pgx,
                                             methods, # input$gx_statmethod
                                             height,
                                             scrollY,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    FDRtable.RENDER <- shiny::reactive({
      methods <- methods() # input$gx_statmethod

      if (is.null(methods)) {
        return(NULL)
      }

      #

      kk <- rownames(pgx$gx.meta$sig.counts[[1]][[1]])
      kk <- intersect(methods, rownames(pgx$gx.meta$sig.counts[[1]][[1]]))
      counts.up <- pgx$gx.meta$sig.counts$up
      counts.down <- pgx$gx.meta$sig.counts$down
      counts.up <- lapply(counts.up, function(x) x[kk, , drop = FALSE])
      counts.down <- lapply(counts.down, function(x) x[kk, , drop = FALSE])
      for (i in 1:length(counts.up)) {
        rownames(counts.up[[i]]) <- paste0(names(counts.up)[i], "::", rownames(counts.up[[i]]))
        rownames(counts.down[[i]]) <- paste0(names(counts.down)[i], "::", rownames(counts.down[[i]]))
      }
      sig.up <- do.call(rbind, counts.up)
      sig.down <- do.call(rbind, counts.down)

      sig.up <- sig.up[order(rownames(sig.up)), , drop = FALSE]
      sig.down <- sig.down[order(rownames(sig.down)), , drop = FALSE]
      colnames(sig.up)[1] <- paste("UP   FDR = ", colnames(sig.up)[1])
      colnames(sig.down)[1] <- paste("DOWN   FDR = ", colnames(sig.down)[1])
      colnames(sig.down) <- paste0("  ", colnames(sig.down))
      sigcount <- cbind(sig.down, sig.up[rownames(sig.down), , drop = FALSE])
      dim(sigcount)
      maxsig <- 0.99 * max(sigcount, na.rm = TRUE)

      contr <- sub("::.*", "", rownames(sigcount))
      #
      metd <- sub(".*::", "", rownames(sigcount))
      D <- data.frame(method = metd, contrast = contr, sigcount, check.names = FALSE)

      DT::datatable(D,
        rownames = FALSE,

        #
        fillContainer = TRUE,
        extensions = c("Scroller"),
        plugins = "scrollResize",
        options = list(
          dom = "lfrtip",
          pageLength = 999, #
          scrollX = TRUE,
          scrollY = scrollY,
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(colnames(sig.up),
          background = DT::styleColorBar(c(0, maxsig), "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) %>%
        DT::formatStyle(colnames(sig.down),
          background = DT::styleColorBar(c(0, maxsig), "lightblue"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    FDRtable.RENDER_modal <- shiny::reactive({
      dt <- FDRtable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "datasets",
      func = FDRtable.RENDER,
      func2 = FDRtable.RENDER_modal,
      selector = "none"
    )
  }) # end module server
} # end server
