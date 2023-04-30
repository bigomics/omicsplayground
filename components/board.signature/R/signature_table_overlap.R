##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

signature_table_overlap_ui <- function(
  id,
  title,
  info.text,
  caption,
  width,
  height) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = "b"
  )
}

signature_table_overlap_server <- function(id,
                                           getOverlapTable,
                                           fullH,
                                           tabH) {
  moduleServer(id, function(input, output, session) {
    overlapTable.RENDER <- shiny::reactive({
      df <- getOverlapTable()
      shiny::req(df)

      df$geneset <- playbase::wrapHyperLink(df$geneset, df$geneset)

      numeric.cols <- which(sapply(df, is.numeric))
      numeric.cols <- intersect(c("p.fisher", "q.fisher"), colnames(df))

      DT::datatable(df,
##      class = "compact cell-border stripe",
        rownames = FALSE, escape = c(-1, -2),
        extensions = c("Scroller"),
        plugins = "scrollResize",
        selection = "none",
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          scrollX = TRUE,
          scrollY = "25vh",
          scrollResize = TRUE,
          scroller = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("score",
          background = playbase::color_from_middle(df$score, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    overlapTable.RENDER_modal <- shiny::reactive({
      dt <- overlapTable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    overlapTable <- TableModuleServer(
      "datasets",
      func = overlapTable.RENDER,
      func2 = overlapTable.RENDER_modal,
      selector = "none"
    )
    return(overlapTable)
  })
}
