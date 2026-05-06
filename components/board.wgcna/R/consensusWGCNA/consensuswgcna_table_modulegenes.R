##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

consensusWGCNA_table_modulegenes_ui <- function(
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

consensusWGCNA_table_modulegenes_server <- function(id,
                                                    mwgcna,
                                                    r_annot,
                                                    r_trait = reactive(NULL),
                                                    r_module = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    table_df <- function() {
      cons <- mwgcna()
      trait <- r_trait()
      module <- r_module()
      annot <- r_annot()

      shiny::req(cons)
      shiny::req(trait)
      shiny::req(module)
      shiny::req(annot)

      if (input$showallmodules) module <- NULL

      stats <- playbase::wgcna.getConsensusGeneStats(cons,
        stats = cons$stats,
        trait = trait, module = module
      )
      ## NULL is returned when the selected trait is absent from at least one layer's
      ## traitSignificance (e.g. constant-within-group traits, sample-level IDs).
      ## This is expected: modTraits (trait dropdown source) and per-layer gene stats
      ## are built from different sample subsets, so not every dropdown trait maps
      ## to computable gene statistics in all layers.
      shiny::validate(shiny::need(!is.null(stats), "Selected trait has no gene statistics. Please select a different trait."))

      cm <- Reduce(intersect, lapply(stats, rownames))
      kk <- c(grep("score\\.", colnames(stats[["full"]]), value = TRUE), "consensus")
      if (input$showallmodules) kk <- c("module", kk)
      d1 <- stats[["full"]][cm, c("feature", kk)]
      d2 <- stats[["consensus"]][cm, c("score", "scorePvalue")]
      df <- cbind(d1, d2)
      rm(d1, d2)

      if (!is.null(annot)) {
        df$title <- playbase::probe2symbol(df$feature, annot, query = "gene_title")
        symbol <- playbase::probe2symbol(df$feature, annot, query = "symbol")
        if (mean(df$feature == symbol, na.rm = TRUE) < 0.2) df$symbol <- symbol
      }

      kk <- grep("score\\.", colnames(df), value = TRUE)
      kk <- c("module", "feature", "symbol", "title", kk, "score", "scorePvalue", "consensus")
      df <- df[, c(intersect(kk, colnames(df))), drop = FALSE]

      return(df)
    }

    render_table <- function(full = TRUE) {
      df <- table_df()
      if ("module" %in% colnames(df)) df$module <- factor(df$module)
      numeric.cols <- which(sapply(df, class) == "numeric")
      score.cols <- grepl("^score", colnames(df)) & !grepl("Pvalue", colnames(df))
      score.vals <- df[, score.cols, drop = FALSE]

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
              targets = c(1), ## without rownames column 2 is target 1
              render = DT::JS("$.fn.dataTable.render.ellipsis( 20, false )")
            ),
            list(
              targets = c(2), ## without rownames column 2 is target 1
              render = DT::JS("$.fn.dataTable.render.ellipsis( 40, false )")
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
