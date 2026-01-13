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
  height
) {
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
    table_data <- shiny::reactive({
      df <- getOverlapTable()
      shiny::req(df)
      return(df)
    })

    overlapTable.RENDER <- shiny::reactive({
      df <- table_data()
      numeric.cols <- intersect(c("p.fisher", "q.fisher"), colnames(df))

      geneset_link <- playbase::wrapHyperLink(
        rep_len("<i class='fa-solid fa-arrow-up-right-from-square weblink'></i>", nrow(df)),
        df$geneset
      ) |> HandleNoLinkFound(
        NoLinkString = "<i class='fa-solid fa-arrow-up-right-from-square'></i>",
        SubstituteString = "<i class='fa-solid fa-arrow-up-right-from-square blank_icon'></i>"
      )

      ## strip class, add link
      df$geneset <- sub(".*[:]", "", df$geneset)
      df$geneset <- paste(df$geneset, geneset_link)
      rownames(df) <- 1:nrow(df)

      DT::datatable(
        df,
        rownames = FALSE,
        escape = FALSE,
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
          background = color_from_middle(df$score, "lightblue", "#f5aeae"),
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
      csvFunc = table_data,
      selector = "none"
    )
    return(overlapTable)
  })
}
