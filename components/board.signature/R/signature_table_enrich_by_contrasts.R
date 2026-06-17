##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

signature_table_enrich_by_contrasts_ui <- function(
  id,
  title,
  info.text,
  caption,
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
    title = title,
    label = "a"
  )
}

signature_get_enrichment_contrasts <- function(enrichmentContrastTable,
                                               selected_only = FALSE) {
  table_data <- enrichmentContrastTable$data()
  if (is.null(table_data) ||
    nrow(table_data) == 0 ||
    !"contrast" %in% colnames(table_data)) {
    return(NULL)
  }

  rows <- enrichmentContrastTable$rows_selected()
  if ((is.null(rows) || length(rows) == 0) && !selected_only) {
    rows <- enrichmentContrastTable$rows_all()
  }
  if (is.null(rows) || length(rows) == 0) {
    return(NULL)
  }

  rows <- as.integer(rows)
  rows <- rows[!is.na(rows) & rows >= 1 & rows <= nrow(table_data)]
  if (length(rows) == 0) {
    return(NULL)
  }

  unname(as.character(table_data$contrast[rows]))
}

signature_table_enrich_by_contrasts_server <- function(id,
                                                       sigCalculateGSEA) {
  moduleServer(id, function(input, output, session) {
    enrichmentContrastTable.RENDER <- shiny::reactive({
      gsea <- sigCalculateGSEA()
      if (is.null(gsea)) {
        return(NULL)
      }

      output <- as.matrix(gsea$output)
      output <- round(output, digits = 4)
      output <- data.frame(contrast = rownames(output), output)
      output <- output[!grepl("^IA:", output$contrast), , drop = FALSE]
      rownames(output) <- output$contrast
      output$q <- NULL
      output$rho <- NULL

      color_fx <- as.numeric(output[, "NES"])
      color_fx[is.na(color_fx)] <- 0 ## yikes...
      numeric.cols <- which(sapply(output, is.numeric))

      DT::datatable(output,
        class = "compact cell-border stripe",
        rownames = FALSE,
        extensions = c("Scroller"),
        plugins = "scrollResize",
        selection = "single",
        fillContainer = TRUE,
        options = list(
          dom = "lrtip",
          scrollX = TRUE,
          scrollY = "20vh",
          scrollResize = TRUE,
          scroller = TRUE,
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
