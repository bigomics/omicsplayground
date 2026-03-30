##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_table_mgsea_ui <- function(
  id,
  label = "",
  title = "",
  info.text = "",
  caption = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("table"),
    info.text = info.text,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = label
  )
}

mofa_table_mgsea_server <- function(id,
                                    mgsea,
                                    input_k = reactive(1),
                                    top = 20) {
  moduleServer(id, function(input, output, session) {
    csvFunc <- function() {
      mgsea <- mgsea()
      validate(need(!is.null(mgsea), "missing GSEA data."))
      k <- 1
      k <- input_k()
      shiny::req(k, mgsea[[k]])
      df <- mgsea[[k]]
      df
    }

    table.RENDER <- function(full = FALSE) {
      mgsea <- mgsea()
      validate(need(!is.null(mgsea), "missing GSEA data."))
      k <- input_k()
      shiny::req(k, mgsea[[k]])
      df <- mgsea[[k]]
      df <- cbind(pathway = rownames(df), df)
      if (!full) {
        num.cols <- grep("^num", colnames(df), value = TRUE)
        df <- df[, c("pathway", "multi.score", "multi.q", num.cols)]
      }
      numeric.cols <- grep("score|pval|p$|q$|rho", colnames(df))

      DT::datatable(
        df,
        rownames = FALSE, #
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        plugins = "scrollResize",
        options = list(
          dom = "lfrtip", #
          scrollX = TRUE, #
          scrollY = "70vh",
          scroller = TRUE,
          scrollResize = TRUE,
          deferRender = TRUE,
          columnDefs = list(
            list(
              targets = "pathway", ## with no rownames column 1 is column 2
              render = DT::JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 100 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                "}"
              )
            )
          )
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    table.RENDER2 <- function(full = FALSE) {
      table.RENDER(full = TRUE)
    }

    table <- TableModuleServer(
      "table",
      func = table.RENDER,
      func2 = table.RENDER2,
      csvFunc = csvFunc,
      selector = "single"
    )

    return(table)
  })
}
