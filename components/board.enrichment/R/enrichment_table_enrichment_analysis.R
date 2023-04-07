##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_table_enrichment_analysis_ui <- function(
  id,
  title,
  info.text,
  caption,
  width,
  height) {
  ns <- shiny::NS(id)

  gseatable_opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("gs_showqvalues"), "show indivivual q-values", FALSE),
      "Show all q-values of each individual statistical method in the table.",
      placement = "top", options = list(container = "body")
    )
  )

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    width = width,
    caption = caption,
    height = height,
    options = gseatable_opts,
    title = title,
    label = "I"
  )
}

enrichment_table_enrichment_analysis_server <- function(id,
                                                        getFilteredGeneSetTable) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    gseatable.RENDER <- shiny::reactive({
      rpt <- getFilteredGeneSetTable()
      if (is.null(rpt)) {
        return(NULL)
      }
      if (nrow(rpt) == 0) {
        return(NULL)
      }

      if (!("GS" %in% colnames(rpt))) rpt <- cbind(GS = rownames(rpt), rpt)
      if ("GS" %in% colnames(rpt)) rpt$GS <- playbase::shortstring(rpt$GS, 72)
      if ("size" %in% colnames(rpt)) rpt$size <- as.integer(rpt$size)

      fx <- NULL
      fx.col <- grep("score|fx|fc|sign|NES|logFC", colnames(rpt))[1]
      if (length(fx.col) > 0) fx <- rpt[, fx.col]

      jj <- which(sapply(rpt, is.numeric))
      if (length(jj) > 0) rpt[, jj] <- round(rpt[, jj], digits = 4)
      jj <- which(sapply(rpt, is.character) | sapply(rpt, is.factor))
      if (length(jj) > 0) rpt[, jj] <- apply(rpt[, jj, drop = FALSE], 2, playbase::shortstring, 100)

      if (!input$gs_showqvalues) {
        rpt <- rpt[, grep("^q[.]|^q$", colnames(rpt), invert = TRUE)]
      }

      ## wrap genesets names with known links.
      rpt$GS <- playbase::wrapHyperLink(rpt$GS, rownames(rpt))
      selectmode <- "single"

      is.numcol <- sapply(rpt, is.numeric)
      numcols <- which(is.numcol & !colnames(rpt) %in% c("size"))
      numcols <- colnames(rpt)[numcols]

      colnames(rpt) <- sub("GS", "geneset", colnames(rpt))

      DT::datatable(rpt,
        class = "compact cell-border stripe hover",
        rownames = FALSE,
        escape = c(-1, -5),
        extensions = c("Scroller"),
        fillContainer = TRUE,
        selection = list(mode = selectmode, target = "row", selected = NULL),
        options = list(
          dom = "frtip",
          paging = TRUE,
          pageLength = 15, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE,
          scrollY = "20vh",
          scroller = TRUE,
          deferRender = TRUE,
          search = list(
            regex = TRUE,
            caseInsensitive = TRUE
          )
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numcols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(fx.col,
          background = playbase::color_from_middle(fx, "lightblue", "#f5aeae")
        )
    })

    gseatable.RENDER_modal <- shiny::reactive({
      dt <- gseatable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    gseatable <- TableModuleServer(
      "datasets",
      func = gseatable.RENDER,
      func2 = gseatable.RENDER_modal,
      selector = "single"
    )

    return(gseatable)
  })
}
