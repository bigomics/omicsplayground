##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_table_gset_enrich_all_contrasts_ui <- function(id, width, height) {
  ns <- shiny::NS(id)

  info_text <- "The <strong>Enrichment (all)</strong> panel reports the gene set enrichment for all contrasts in the selected dataset."

  TableModuleUI(
    ns("datasets"),
    info.text = info_text,
    width = width,
    height = height,
    title = "Gene set enrichment for all contrasts"
  )
}

enrichment_table_gset_enrich_all_contrasts_server <- function(id,
                                                              pgx,
                                                              getFilteredGeneSetTable) {
  moduleServer(id, function(input, output, session) {
    tabH <- 340 ## row height of panels

    fctable.RENDER <- shiny::reactive({

      ## get all contrasts
      F <- sapply(pgx$gset.meta$meta, function(x) x[, "meta.fx"])
      colnames(F) <- gsub("_", " ", colnames(F))
      rownames(F) <- rownames(pgx$gset.meta$meta[[1]])
      fc.var <- round(rowMeans(F**2, na.rm = TRUE), digits = 3)
      gs <- substring(rownames(F), 1, 60)
      F1 <- data.frame(geneset = gs, fc.var = fc.var, round(F, digits = 3), check.names = FALSE)

      ## get current filtered geneset and extract names of gene sets
      rpt <- getFilteredGeneSetTable()
      F1 <- F1[intersect(rownames(rpt), rownames(F1)), , drop = FALSE]
      F1$geneset <- playbase::wrapHyperLink(F1$geneset, rownames(F1))

      DT::datatable(F1,
        rownames = FALSE, escape = -1,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          scrollX = TRUE,
          scrollY = "20vh",
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("fc.var",
          background = playbase::color_from_middle(fc.var, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) %>%
        DT::formatStyle(colnames(F),
          background = playbase::color_from_middle(F[, ], "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
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
  })
}
