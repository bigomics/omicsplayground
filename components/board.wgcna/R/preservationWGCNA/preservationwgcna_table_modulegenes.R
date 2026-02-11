##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

preservationWGCNA_table_modulegenes_ui <- function(
  id,
  label = "a",
  title = "Title",
  info.text = "Info",
  caption = "Caption",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("showallmodules"),
      label = "Show all modules",
      value = FALSE
    )
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

preservationWGCNA_table_modulegenes_server <- function(id,
                                                       rwgcna,
                                                       rannot,
                                                       rtrait = reactive(NULL),
                                                       rmodule = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    table_df <- function() {
      pres <- rwgcna()

      trait <- rtrait()
      module <- rmodule()
      annot <- rannot()

      shiny::req(pres)
      shiny::req(trait)
      shiny::req(module)
      shiny::req(annot)

      refname <- names(pres$layers)[1]
      refname

      if (input$showallmodules) module <- NULL
      df <- playbase::wgcna.getGeneStats(
        wgcna = NULL,
        stats = pres$stats,
        trait = trait,
        module = module,
        labels = pres$colors[, refname],
        plot = FALSE
      )

      if (!is.null(annot)) {
        df$title <- playbase::probe2symbol(df$feature, annot, query = "gene_title")
        symbol <- playbase::probe2symbol(df$feature, annot, query = "symbol")
        if (mean(df$feature == symbol, na.rm = TRUE) < 0.2) df$symbol <- symbol
      }

      score.sign <- sign(mean(df$score, na.rm = TRUE))
      df <- df[order(-df$score * score.sign), ]

      return(df)
    }

    render_table <- function(full = TRUE) {
      df <- table_df()

      ## set correct types for filter
      df$module <- factor(df$module)

      score.cols <- grepl("^score", colnames(df)) & !grepl("Pvalue", colnames(df))
      if (!full) {
        cols <- c("module", "feature", "symbol", "title")
        cols <- c(cols, colnames(df)[which(score.cols)])
        cols <- intersect(cols, colnames(df))
        if ("symbol" %in% cols) cols <- setdiff(cols, "title")
        df <- df[, cols]
      }

      if (!input$showallmodules) df$module <- NULL

      ## order name first
      cols <- unique(c("module", "feature", "symbol", "title", colnames(df)))
      cols <- intersect(cols, colnames(df))
      df <- df[, cols]

      ## rename
      colnames(df) <- sub("moduleMembership", "MM", colnames(df))
      colnames(df) <- sub("traitSignificance", "TS", colnames(df))

      numeric.cols <- which(sapply(df, class) == "numeric")
      score.cols <- grepl("^score", colnames(df)) & !grepl("Pvalue", colnames(df))
      score.vals <- df[, score.cols, drop = FALSE]
      short60.cols <- intersect(c("title"), colnames(df))
      short20.cols <- intersect(c("feature", "symbol"), colnames(df))

      DT::datatable(
        df,
        rownames = FALSE, #
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        plugins = c("scrollResize", "ellipsis"),
        # filter = 'top',
        options = list(
          dom = "lfrtip", #
          scrollX = TRUE, #
          scrollY = "70vh",
          scroller = TRUE,
          scrollResize = TRUE,
          deferRender = TRUE,
          columnDefs = list(
            list(
              # targets = c(0,1), ## without rownames column 2 is target 1
              targets = short60.cols,
              render = DT::JS("$.fn.dataTable.render.ellipsis( 60, false )")
            ),
            list(
              # targets = c(0,1), ## without rownames column 2 is target 1
              targets = short20.cols,
              render = DT::JS("$.fn.dataTable.render.ellipsis( 20, false )")
            )
          )
        ) ## end of options
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(
          score.cols,
          background = color_from_middle(score.vals, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }

    table.RENDER <- function() {
      render_table(full = FALSE)
    }

    table.RENDER2 <- function() {
      render_table(full = TRUE)
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
