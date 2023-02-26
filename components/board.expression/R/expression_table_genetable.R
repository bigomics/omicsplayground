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
expression_table_genetable_ui <- function(id, width, height) {
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

  genetable_text <- "Table <strong>I</strong> shows the results of the statistical tests. To increase the statistical reliability of the Omics Playground, we perform the DE analysis using four commonly accepted methods in the literature, namely, T-test (standard, Welch), <a href='https://www.ncbi.nlm.nih.gov/pubmed/25605792'> limma</a> (no trend, trend, voom), <a href='https://www.ncbi.nlm.nih.gov/pubmed/19910308'> edgeR</a> (QLF, LRT), and <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4302049'> DESeq2</a> (Wald, LRT), and merge the results.
<br><br>For a selected comparison under the <code>Contrast</code> setting, the results of the selected methods are combined and reported under the table, where <code>meta.q</code> for a gene represents the highest <code>q</code> value among the methods and the number of stars for a gene indicate how many methods identified significant <code>q</code> values (<code>q < 0.05</code>). The table is interactive (scrollable, clickable); users can sort genes by <code>logFC</code>, <code>meta.q</code>, or average expression in either conditions. Users can filter top N = {10} differently expressed genes in the table by clicking the <code>top 10 genes</code> from the table <i>Settings</i>."

  TableModuleUI(
    ns("datasets"),
    info.text = genetable_text,
    width = width,
    height = height,
    options = genetable_opts,
    title = "Differential expression analysis",
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

    table.RENDER <- shiny::reactive({
      res <- res()

      if (is.null(res) || nrow(res) == 0) {
        return(NULL)
      }

      fx.col <- grep("fc|fx|mean.diff|logfc|foldchange", tolower(colnames(res)))[1]
      fx.col
      fx <- res[, fx.col]

      if ("gene_title" %in% colnames(res)) res$gene_title <- shortstring(res$gene_title, 50)
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
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          # paging = TRUE,
          # pageLength = 16, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE,
          scrollY = "20vh",
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
          background = color_from_middle(fx, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    table.RENDER_modal <- shiny::reactive({
      dt <- table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    genetable <- TableModuleServer(
      "datasets",
      func = table.RENDER,
      func2 = table.RENDER_modal,
      selector = "single"
    )

    return(genetable)
  })
}
