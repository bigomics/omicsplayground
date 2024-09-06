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

    table_data <- shiny::reactive({
      rpt <- getFilteredGeneSetTable()
      if (!("GS" %in% colnames(rpt))) rpt <- cbind(GS = rownames(rpt), rpt)
      if ("matched.genes" %in% names(rpt)) {
        names(rpt)[names(rpt) == "matched.genes"] <- "size"
      }
      rpt
    })

    gseatable.RENDER <- function() {
      rpt <- table_data()

      if (is.null(rpt)) {
        return(NULL)
      }
      if (nrow(rpt) == 0) {
        return(NULL)
      }

      if ("GS" %in% colnames(rpt)) rpt$GS <- playbase::shortstring(rpt$GS, 72)
      if ("size" %in% colnames(rpt)) rpt$size <- as.integer(rpt$size)

      fx <- NULL
      fx.col <- grep("score|fx|fc|sign|NES|logFC", colnames(rpt))[1]
      if (length(fx.col) > 0) fx <- rpt[, fx.col]

      jj <- which(sapply(rpt, function(col) is.numeric(col) && !is.integer(col)))
      if (length(jj) > 0) rpt[, jj] <- round(rpt[, jj], digits = 4)
      jj <- which(sapply(rpt, is.character) | sapply(rpt, is.factor))
      if (length(jj) > 0) rpt[, jj] <- apply(rpt[, jj, drop = FALSE], 2, playbase::shortstring, 100)
      if (!input$gs_showqvalues) {
        rpt <- rpt[, grep("^q[.]|^q$", colnames(rpt), invert = TRUE)]
      }

      ## wrap genesets names with known links.
      GS_link <- playbase::wrapHyperLink(
        rep_len("<i class='fa-solid fa-arrow-up-right-from-square weblink'></i>", nrow(rpt)),
        rownames(rpt)
      ) |> HandleNoLinkFound(
        NoLinkString = "<i class='fa-solid fa-arrow-up-right-from-square weblink'></i>",
        SubstituteString = "<i class='fa-solid fa-arrow-up-right-from-square blank_icon'></i>"
      )
      rpt$GS <- paste(rpt$GS, "&nbsp;", GS_link)
      colnames(rpt) <- sub("GS", "geneset", colnames(rpt))

      is.numcol <- sapply(rpt, function(col) is.numeric(col) && !is.integer(col))
      numcols <- which(is.numcol & !colnames(rpt) %in% c("size"))
      numcols <- colnames(rpt)[numcols]


      DT::datatable(rpt,
        class = "compact cell-border stripe hover",
        rownames = FALSE,
        escape = c(-1, -2),
        extensions = c("Scroller"),
        plugins = "scrollResize",
        fillContainer = TRUE,
        selection = list(mode = "single", target = "row", selected = 1),
        options = list(
          dom = "frtip",
          paging = TRUE,
          pageLength = 15, #
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
        DT::formatSignif(numcols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(fx.col,
          background = color_from_middle(fx, "lightblue", "#f5aeae")
        )
    }

    gseatable.RENDER_modal <- shiny::reactive({
      dt <- gseatable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    table_data_csv <- function() {
      df <- table_data()
      df$stars <- NULL
      return(df)
    }

    gseatable <- TableModuleServer(
      "datasets",
      func = gseatable.RENDER,
      func2 = gseatable.RENDER_modal,
      csvFunc = table_data_csv,
      selector = "single"
    )

    return(gseatable)
  })
}
