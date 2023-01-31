##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

signature_table_enrich_by_contrasts_ui <- function(id) {
  ns <- shiny::NS(id)

  tableWidget(ns("table"))
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

      dbg("enrichmentContrastTable.RENDER: reacted")

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
      numeric.cols

      DT::datatable(output,
        class = "compact cell-border stripe",
        rownames = FALSE,
        extensions = c("Scroller"),
        selection = "single",
        fillContainer = TRUE,
        options = list(
          dom = "lrtip",
          scrollX = TRUE, scrollY = tabH, scroller = TRUE,
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

    info.text1 <- "<b>Enrichment by contrast.</b> Enrichment scores of query signature across all contrasts. The table summarizes the enrichment statistics of the gene list in all contrasts using the GSEA algorithm. The NES corresponds to the normalized enrichment score of the GSEA analysis.  "

    enrichmentContrastTable <- shiny::callModule(
      tableModule,
      id = "table",
      func = enrichmentContrastTable.RENDER,
      info.text = info.text1,
      caption2 = info.text1,
      title = tags$div(
        HTML('<span class="module-label">(a)</span>Enrichment by contrasts')
      ),
      height = c(230, 700)
    )
    return(enrichmentContrastTable)
  })
}
