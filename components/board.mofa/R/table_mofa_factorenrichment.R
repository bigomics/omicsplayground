##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_table_factorenrichment_ui <- function(
  id,
  label = "",
  title = "",
  info.text = "",
  caption = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::checkboxInput(ns("onlyshared"), "only shared pathways", FALSE)
  )

  TableModuleUI(
    ns("table"),
    info.text = info.text,
    options = options,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = label
  )
}

mofa_table_factorenrichment_server <- function(id,
                                               gsea,
                                               selected_factor = reactive(1)) {
  moduleServer(id, function(input, output, session) {
    table.RENDER <- function(full = FALSE) {
      gsea <- gsea()

      validate(need(!is.null(table), "missing GSEA data."))
      k <- 1
      k <- selected_factor() ## which factor/phenotype
      shiny::req(k)

      df <- gsea$table[[k]]
      df <- df[order(-df$NES), ]
      df <- df[, grep("pathway|NES|pval|padj|size|leadingEdge", colnames(df))]
      if (input$onlyshared) {
        size.cols <- grep("^size", colnames(df))
        sel <- apply(df[, size.cols], 1, function(x) all(x > 0))
        df <- df[sel, ]
      }
      df <- data.frame(factor = k, df, check.names = FALSE)
      if (!full) {
        ## filter out pvalue columns
        df <- df[, grep("pval|leadingEdge", colnames(df), invert = TRUE)]
      }

      numeric.cols <- grep("NES|score|pval|padj|rho", colnames(df))

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
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
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
                "return type === 'display' && data.length > 80 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 80) + '...</span>' : data;",
                "}"
              )
            )
          )
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(
          "NES",
          background = color_from_middle(df$NES, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }

    table.RENDER2 <- function(full = FALSE) {
      table.RENDER(full = TRUE)
    }

    table <- TableModuleServer(
      "table",
      func = table.RENDER,
      func2 = table.RENDER2,
      selector = "single"
    )

    return(table)
  })
}
