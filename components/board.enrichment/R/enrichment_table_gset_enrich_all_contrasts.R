##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

enrichment_table_gset_enrich_all_contrasts_ui <- function(id) {
  ns <- shiny::NS(id)

  tableWidget(ns("fctable"))
}

enrichment_table_gset_enrich_all_contrasts_server <- function(id,
                                                              inputData,
                                                              getFilteredGeneSetTable) {
  moduleServer(id, function(input, output, session) {
    tabH <- 340 ## row height of panels

    fctable.RENDER <- shiny::reactive({
      ngs <- inputData()

      ## get all contrasts
      F <- sapply(ngs$gset.meta$meta, function(x) x[, "meta.fx"])
      colnames(F) <- gsub("_", " ", colnames(F))
      rownames(F) <- rownames(ngs$gset.meta$meta[[1]])
      fc.var <- round(rowMeans(F**2, na.rm = TRUE), digits = 3)
      gs <- substring(rownames(F), 1, 60)
      F1 <- data.frame(geneset = gs, fc.var = fc.var, round(F, digits = 3), check.names = FALSE)

      ## get current filtered geneset and extract names of gene sets
      rpt <- getFilteredGeneSetTable()
      F1 <- F1[intersect(rownames(rpt), rownames(F1)), , drop = FALSE]
      F1$geneset <- wrapHyperLink(F1$geneset, rownames(F1))

      DT::datatable(F1,
        rownames = FALSE, escape = -1,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          scrollX = TRUE,
          scrollY = tabH,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("fc.var",
          background = color_from_middle(fc.var, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) %>%
        DT::formatStyle(colnames(F),
          background = color_from_middle(F[, ], "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    gx_fctable_text <- "The <strong>Enrichment (all)</strong> panel reports the gene set enrichment for all contrasts in the selected dataset."

    gx_fctable_caption <- "<b>Enrichment for all contrasts.</b> Table summarizing the enrichment for all gene sets across all contrasts. The column `fc.var` corresponds to the variance of the gene set across all contrasts."

    shiny::callModule(
      tableModule,
      id = "fctable",
      func = fctable.RENDER,
      title = "Gene set enrichment for all contrasts",
      info.text = gx_fctable_text,
      caption = gx_fctable_caption,
      height = c(295, 750),
      width = c("100%", 1600)
    )
  })
}
