##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_table_gset_enrich_all_contrasts_ui <- function(
  id,
  title,
  info.text,
  caption,
  width,
  height) {
  ns <- shiny::NS(id)

  fctable_opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("showq"), "show q-values", FALSE),
      "Show q-values next to enrichment values.",
      placement = "right", options = list(container = "body")
    )
  )
  
  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    width = width,
    height = height,
    options = fctable_opts,
    title = title,
    caption = caption
  )
}

enrichment_table_gset_enrich_all_contrasts_server <- function(id,
                                                              pgx,
                                                              getFilteredGeneSetTable,
                                                              metaFC,
                                                              metaQ) {
  moduleServer(id, function(input, output, session) {
    tabH <- 340 ## row height of panels

    fctable.RENDER <- shiny::reactive({

      ## get enrichment for all contrasts
      ##F <- sapply(pgx$gset.meta$meta, function(x) x[, "meta.fx"])
      F <- metaFC()
      colnames(F) <- gsub("_", " ", colnames(F))
      ##rownames(F) <- rownames(pgx$gset.meta$meta[[1]])

      ## get enrichment q-value for all contrasts
      ##Q <- sapply(pgx$gset.meta$meta, function(x) x[,"meta.q"])
      Q <- metaQ()
      colnames(Q) <- gsub("_", " ", colnames(Q))
      ##rownames(Q) <- rownames(pgx$gset.meta$meta[[1]])

      ## RMS (non-centered variance)
      fc.rms <- sqrt(F[, 1]**2)
      if (NCOL(F) > 1) {
        fc.rms <- round(sqrt(rowMeans(F**2, na.rm = TRUE)), digits = 3)
      }
      
      ## show Q values?
      show.q <- TRUE
      show.q <- input$showq
      df <- NULL
      gs <- substring(rownames(F), 1, 60)
      if (show.q) {
        F1 <- do.call(cbind, lapply(1:ncol(F), function(i) cbind(F[, i], Q[, i])))
        colnames(F1) <- as.vector(rbind(paste0("ES.", colnames(F)), paste0("q.", colnames(Q))))
        ## colnames(F1) <- sub("q.*","q",colnames(F1))
        df <- data.frame(geneset = gs, rms.ES = fc.rms, F1, check.names = FALSE)
      } else {
        F1 <- F
        colnames(F1) <- paste0("ES.", colnames(F))
        df <- data.frame(geneset = gs, rms.ES = fc.rms, F1, check.names = FALSE)
      }

      ## get current filtered geneset and extract names of gene sets
      res <- getFilteredGeneSetTable()
      df <- df[intersect(rownames(df), rownames(res)), ] ## take intersection of current comparison
      df <- df[order(-df$rms.ES), ]
      colnames(df) <- gsub("_", " ", colnames(df)) ## so it allows wrap line
      colnames(F1) <- gsub("_", " ", colnames(F1)) ## so it allows wrap line
      qv.cols <- grep("^q", colnames(df))
      fc.cols <- setdiff(which(colnames(df) %in% colnames(F1)), qv.cols)
      ## if(length(qv.cols)==0) qv = 0
      df$geneset <- wrapHyperLink(df$geneset, rownames(df))
      
      dt <- DT::datatable(
        df,
        rownames = FALSE,
        escape = -1,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          scrollX = TRUE,
          scrollY = "20vh",
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatSignif(columns = fc.cols, digits = 4) %>%
        DT::formatStyle(
          "rms.ES",
          background = playbase::color_from_middle(fc.rms, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) %>%
        DT::formatStyle(
          fc.cols,
          background = playbase::color_from_middle(F[, ], "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )

      if (length(qv.cols) > 0) {
        dt <- dt %>%
          DT::formatSignif(columns = qv.cols, digits = 4)
      }
      return(dt)
    })

    fctable.RENDER_modal <- shiny::reactive({
      dt <- fctable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      return(dt)
    })

    TableModuleServer(
      "datasets",
      func = fctable.RENDER,
      func2 = fctable.RENDER_modal,
      selector = "none"
    )
  })
}
