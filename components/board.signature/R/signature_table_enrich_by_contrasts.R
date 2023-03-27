##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

signature_table_enrich_by_contrasts_ui <- function(id, width, height) {
  ns <- shiny::NS(id)

  info_text <- "<b>Enrichment by contrast.</b> Enrichment scores of query signature across all contrasts. The table summarizes the enrichment statistics of the gene list in all contrasts using the GSEA algorithm. The NES corresponds to the normalized enrichment score of the GSEA analysis.  "

  TableModuleUI(
    ns("datasets"),
    info.text = info_text,
    width = width,
    height = height,
    title = "Enrichment by contrasts",
    label = "a"
  )
}

signature_table_enrich_by_contrasts_server <- function(id,
                                                       sigCalculateGSEA,
                                                       tabH) {
  moduleServer(id, function(input, output, session) {
    enrichmentContrastTable.RENDER <- shiny::reactive({
      gsea <- sigCalculateGSEA()
      if (is.null(gsea)) {
        return(NULL)
      }

      output <- as.matrix(gsea$output)
      output <- round(output, digits = 4)
      output <- data.frame(contrast = rownames(output), output)
      if (!DEV) {
        output$p <- NULL
        output$rho <- NULL
      }

      color_fx <- as.numeric(output[, "NES"])
      color_fx[is.na(color_fx)] <- 0 ## yikes...
      numeric.cols <- which(sapply(output, is.numeric))

      DT::datatable(output,
        class = "compact cell-border stripe",
        rownames = FALSE,
        extensions = c("Scroller"),
        selection = "single",
        fillContainer = TRUE,
        options = list(
          dom = "lrtip",
          scrollX = TRUE, scrollY = "20vh", scroller = TRUE,
          deferRender = FALSE
        )
      ) %>% ## end of options.list
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("NES",
          background = color_from_middle(color_fx, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    enrichmentContrastTable.RENDER_render <- shiny::reactive({
      dt <- enrichmentContrastTable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    enrichmentContrastTable <- TableModuleServer(
      "datasets",
      func = enrichmentContrastTable.RENDER,
      func2 = enrichmentContrastTable.RENDER_render,
      selector = "single"
    )

    return(enrichmentContrastTable)
  })
}
