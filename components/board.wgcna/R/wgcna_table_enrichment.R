##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_table_enrichment_ui <- function(
  id,
  label = "",
  title = "",
  info.text = "",
  caption = "",
  width,
  height
) {
  ns <- shiny::NS(id)

  options <- tagList(
    checkboxInput(ns("showallmodules"), "Show all modules")
  )

  TableModuleUI(
    ns("datasets"),
    title = title,
    options = options,
    info.text = info.text,
    caption = caption,
    width = width,
    height = height,
    label = label
  )
}

wgcna_table_enrichment_server <- function(id,
                                          wgcna,
                                          selected_module
                                          ## enrich_table
) {
  moduleServer(id, function(input, output, session) {
    table_data <- function() {
      gse <- wgcna()$gse
      k <- selected_module()
      if (input$showallmodules) k <- "<all>"
      if (k %in% names(gse)) {
        df <- gse[[k]]
      } else {
        df <- do.call(rbind, gse)
      }
      if (!"score" %in% colnames(df)) {
        df$odd.ratio[is.infinite(df$odd.ratio)] <- 99
        df$score <- df$odd.ratio * -log10(df$p.value)
      }
      cols <- c(
        "module", "geneset", "score", "p.value", "q.value",
        "overlap", "genes"
      )
      cols <- intersect(cols, colnames(df))
      df <- df[, cols]
      df <- df[order(-df$score), ]
      df
    }

    render_table <- function(full = FALSE) {
      df <- table_data()
      if (!full) {
        cols <- c("geneset", "score", "q.value", "overlap", "genes")
        cols <- intersect(cols, colnames(df))
        df <- df[, cols]
      }

      numeric.cols <- grep("score|value|ratio", colnames(df))

      DT::datatable(
        df,
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = "70vh",
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE,
          columnDefs = list(
            list(
              targets = "geneset", ## with no rownames column 1 is column 2
              render = DT::JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 60 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 60) + '...</span>' : data;",
                "}"
              )
            )
          )
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "10px", lineHeight = "70%") %>%
        DT::formatStyle("score",
          background = color_from_middle(df$score, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }

    RENDER <- function() {
      render_table(full = FALSE)
    }

    RENDER_modal <- function() {
      dt <- render_table(full = TRUE)
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    }

    tablemodule <- TableModuleServer(
      "datasets",
      func = RENDER,
      func2 = RENDER_modal,
      csvFunc = table_data,
      selector = "single"
    )

    return(tablemodule)
  })
}
