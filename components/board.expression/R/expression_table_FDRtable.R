##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' UI code for table code: expression board
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
expression_table_FDRtable_ui <- function(id) {

  ns <- shiny::NS(id)

  tableWidget(ns("FDRtable"))

}

#' Server side table code: expression board
#'
#' @param id
#' @param watermark
#'
#' @export
expression_table_FDRtable_server <- function(id,
                                             ngs,
                                             GX.DEFAULTTEST,
                                             height, #c(tabH, 700)
                                             watermark=FALSE){
  moduleServer( id, function(input, output, session) {

    ns <- session$ns

    FDRtable.RENDER <- shiny::reactive({

      methods <- GX.DEFAULTTEST
      methods <- input$gx_statmethod
      ## methods = input$gx_statmethod
      if (is.null(methods)) {
        return(NULL)
      }

      ## comp <- input$gx_contrast
      ngs <- ngs()

      kk <- rownames(ngs$gx.meta$sig.counts[[1]][[1]])
      kk <- intersect(methods, rownames(ngs$gx.meta$sig.counts[[1]][[1]]))
      counts.up <- ngs$gx.meta$sig.counts$up
      counts.down <- ngs$gx.meta$sig.counts$down
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
      ## contr = rownames(sigcount)
      metd <- sub(".*::", "", rownames(sigcount))
      D <- data.frame(method = metd, contrast = contr, sigcount, check.names = FALSE)

      DT::datatable(D,
                    rownames = FALSE,
                    #                      class = 'compact cell-border stripe hover',
                    class = "compact hover",
                    fillContainer = TRUE,
                    extensions = c("Scroller"),
                    options = list(
                      dom = "lfrtip",
                      pageLength = 999, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                      scrollX = TRUE,
                      scrollY = tabV,
                      scroller = TRUE, deferRender = TRUE
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

  FDRtable_text <- "The <strong>FDR table</strong> tab reports the number of significant genes at different FDR thresholds for all contrasts within the dataset."

  FDRtable_caption <- "<b>Number of significant genes versus FDR.</b> This table reports the number of significant genes at different FDR thresholds for all contrasts and methods. This enables to quickly see which methods are more sensitive. The left part of the table (in blue) correspond to the number of significant down-regulated genes, the right part (in red) correspond to the number of significant overexpressed genes."

  shiny::callModule(
    tableModule,
    id = "FDRtable",
    func = FDRtable.RENDER,
    info.text = FDRtable_text,
    title = "Number of significant genes",
    caption = FDRtable_caption,
    height = height
  )
  })#end module server
}#end server