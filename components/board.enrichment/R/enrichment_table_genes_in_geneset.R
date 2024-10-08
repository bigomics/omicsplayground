##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_table_genes_in_geneset_ui <- function(
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
    label = "II"
  )
}

enrichment_table_genes_in_geneset_server <- function(id,
                                                     organism,
                                                     geneDetails) {
  moduleServer(id, function(input, output, session) {
    genetable.RENDER <- shiny::reactive({
      rpt <- geneDetails()
      if (is.null(rpt) || nrow(rpt) == 0) {
        shiny::validate(shiny::need(
          nrow(rpt) > 0,
          tspan("Please select a geneset from the table on the left to view genes.", js = FALSE)
        ))
        return(NULL)
      }

      if (organism %in% c("Human", "human")) {
        rpt$human_ortholog <- NULL
      }
      if (sum(rpt$feature %in% rpt$symbol) > nrow(rpt) * .8) {
        rpt$feature <- NULL
      }

      # rpt$gene_title <- NULL  # this is important for metabolomics
      if (!is.null(rpt) && nrow(rpt) > 0) {
        jj <- which(sapply(rpt, is.numeric))
        rpt[, jj] <- round(rpt[, jj], digits = 4)
        jj <- which(sapply(rpt, is.character) | sapply(rpt, is.factor))
        if (length(jj) > 0) rpt[, jj] <- apply(rpt[, jj, drop = FALSE], 2, playbase::shortstring, 60)
      } else {
        rpt <- data.frame("", 0, 0, 0)[0, ]
        colnames(rpt) <- c("gene_name", "fc", "p", "q")
      }

      colnames(rpt) <- sub("^GS$", "gene set", colnames(rpt))
      numeric.cols <- which(sapply(rpt, is.numeric))

      tbl <- DT::datatable(rpt,
        class = "compact cell-border stripe", rownames = FALSE,
        extensions = c("Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          # pageLength = 15, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE,
          scrollY = "calc(45vh - 260px)",
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE,
          search = list(
            regex = TRUE,
            caseInsensitive = TRUE
          )
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4)

      if (nrow(rpt) > 0 && ("fc" %in% colnames(rpt))) {
        fx <- rpt[, "fc"]
        tbl <- tbl %>%
          DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
          DT::formatStyle("fc", background = color_from_middle(fx, "lightblue", "#f5aeae"))
      }
      tbl
    })

    genetable.RENDER_modal <- shiny::reactive({
      dt <- genetable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    genetable <- TableModuleServer(
      "datasets",
      func = genetable.RENDER,
      func2 = genetable.RENDER_modal,
      selector = "single"
    )

    return(genetable)
  })
}
