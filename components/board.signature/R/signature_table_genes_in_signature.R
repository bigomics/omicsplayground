##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

signature_table_genes_in_signature_ui <- function(
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
    width = width,
    height = height,
    title = title,
    caption, caption,
    label = "b"
  )
}

signature_table_genes_in_signature_server <- function(id,
                                                      organism,
                                                      getEnrichmentGeneTable) {
  moduleServer(id, function(input, output, session) {
    enrichmentGeneTable.RENDER <- shiny::reactive({
      df <- getEnrichmentGeneTable()
      if (is.null(df)) {
        shiny::validate(shiny::need(!is.null(df), "Please select a contrast"))
        return(NULL)
      }

      if (organism %in% c("Human", "human")) {
        df$human_ortholog <- NULL
      }
      if (sum(df$feature %in% df$symbol) > nrow(df) * .8) {
        df$feature <- NULL
      }

      num_cols <- sapply(df, is.numeric)
      color_fx <- df[, num_cols]
      color_fx[is.na(color_fx)] <- 0 ## yikes...
      numeric.cols <- colnames(df)[num_cols]

      DT::datatable(df,
        class = "compact cell-border stripe",
        rownames = FALSE,
        extensions = c("Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = NULL),
        fillContainer = TRUE,
        options = list(
          dom = "lrftip",
          scrollX = TRUE,
          scrollY = "30vh",
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = FALSE
        )
      ) %>% ## end of options.list
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(
          numeric.cols,
          background = color_from_middle(color_fx, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    enrichmentGeneTable.RENDER_modal <- shiny::reactive({
      dt <- enrichmentGeneTable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    enrichmentGeneTable <- TableModuleServer(
      "datasets",
      func = enrichmentGeneTable.RENDER,
      func2 = enrichmentGeneTable.RENDER_modal,
      selector = "single"
    )
    return(enrichmentGeneTable)
  })
}
