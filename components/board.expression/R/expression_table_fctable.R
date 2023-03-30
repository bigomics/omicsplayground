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
expression_table_fctable_ui <- function(id, width, height) {
  ns <- shiny::NS(id)

  fctable_text <- "The <strong>Foldchange (all)</strong> tab reports the gene fold changes for all contrasts in the selected dataset."

  fctable_opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("fctable_showq"), "show q-values", TRUE),
      "Show q-values next to FC values.",
      placement = "right", options = list(container = "body")
    )
  )

  TableModuleUI(
    ns("datasets"),
    info.text = fctable_text,
    width = width,
    height = height,
    options = fctable_opts,
    title = "Gene fold changes for all contrasts"
  )
}

#' Server side table code: expression board
#'
#' @param id
#' @param watermark
#'
#' @export
expression_table_fctable_server <- function(id,
                                            pgx,
                                            res, # filteredDiffExprTable
                                            metaFC,
                                            metaQ,
                                            height,
                                            tabV,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    fctable.RENDER <- shiny::reactive({
      res <- res()
      if (is.null(res) || nrow(res) == 0) {
        return(NULL)
      }

      ## F <- sapply(pgx$gx.meta$meta, function(x) unclass(x$fc)[,"trend.limma"])
      ## Q <- sapply(pgx$gx.meta$meta, function(x) x$meta.q)
      ## F <- sapply(pgx$gx.meta$meta, function(x) x$meta.fx)
      ## rownames(F)=rownames(Q)=rownames(pgx$gx.meta$meta[[1]])
      F <- metaFC()
      Q <- metaQ()

      fc.rms <- sqrt(F[, 1]**2)
      if (NCOL(F) > 1) {
        fc.rms <- round(sqrt(rowMeans(F**2)), digits = 4)
      }

      show.q <- TRUE
      show.q <- input$fctable_showq
      df <- NULL
      if (show.q) {
        F1 <- do.call(cbind, lapply(1:ncol(F), function(i) cbind(F[, i], Q[, i])))
        colnames(F1) <- as.vector(rbind(paste0("FC.", colnames(F)), paste0("q.", colnames(Q))))
        ## colnames(F1) <- sub("q.*","q",colnames(F1))
        df <- data.frame(gene = rownames(F), rms.FC = fc.rms, F1, check.names = FALSE)
      } else {
        F1 <- F
        colnames(F1) <- paste0("FC.", colnames(F))
        df <- data.frame(gene = rownames(F), rms.FC = fc.rms, F1, check.names = FALSE)
      }

      df <- df[intersect(rownames(df), rownames(res)), ] ## take intersection of current comparison
      df <- df[order(-df$rms.FC), ]
      colnames(df) <- gsub("_", " ", colnames(df)) ## so it allows wrap line
      colnames(F1) <- gsub("_", " ", colnames(F1)) ## so it allows wrap line
      qv.cols <- grep("^q", colnames(F1))
      fc.cols <- setdiff(which(colnames(df) %in% colnames(F1)), qv.cols)
      ## if(length(qv.cols)==0) qv = 0

      dt <- DT::datatable(df,
        rownames = FALSE,
        # class = 'compact cell-border stripe hover',
        class = "compact hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = c(1)),
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE,
          scrollY = "220",
          scroller = TRUE, deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatSignif(columns = fc.cols, digits = 3) %>%
        DT::formatStyle("rms.FC",
          ## background = DT::styleColorBar(c(0,3), 'lightblue'),
          background = color_from_middle(fc.rms, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) %>%
        DT::formatStyle(fc.cols,
          ## background = DT::styleColorBar(c(0,3), 'lightblue'),
          background = color_from_middle(F, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )

      if (length(qv.cols) > 0) {
        dt <- dt %>%
          DT::formatSignif(columns = qv.cols, digits = 3)
      }
      dt
    })

    fctable.RENDER_modal <- shiny::reactive({
      dt <- fctable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "datasets",
      func = fctable.RENDER,
      func2 = fctable.RENDER_modal,
      selector = "none"
    )
  }) # end module server
} # end server
