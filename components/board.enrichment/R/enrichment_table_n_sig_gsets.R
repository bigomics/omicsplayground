##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_table_n_sig_gsets_ui <- function(
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
    caption = caption,
    title = title
  )
}

enrichment_table_n_sig_gsets_server <- function(id,
                                                pgx,
                                                gs_statmethod) {
  moduleServer(id, function(input, output, session) {
    tabH <- 340 ## row height of panels

    FDRtable.RENDER <- shiny::reactive({
      shiny::req(pgx$X, gs_statmethod())

      meta <- pgx$gset.meta
      test <- gs_statmethod()

      if (length(test) == 1) {
        sig.up <- meta$sig.counts[[test]][["up"]]
        sig.down <- meta$sig.counts[[test]][["down"]]
        rownames(sig.up) <- paste0(rownames(sig.up), "::", test[1])
        rownames(sig.down) <- paste0(rownames(sig.down), "::", test[1])
      } else {
        sig.up <- c()
        sig.down <- c()
        for (i in 1:length(test)) {
          sig1 <- meta$sig.counts[[test[i]]][["up"]]
          sig2 <- meta$sig.counts[[test[i]]][["down"]]
          rownames(sig1) <- paste0(rownames(sig1), "::", test[i])
          rownames(sig2) <- paste0(rownames(sig2), "::", test[i])
          sig.up <- rbind(sig.up, sig1)
          sig.down <- rbind(sig.down, sig2)
        }
      }
      sig.up <- sig.up[order(rownames(sig.up)), , drop = FALSE]
      sig.down <- sig.down[order(rownames(sig.down)), , drop = FALSE]
      pvals <- sort(c(1e-16, 10**seq(-8, -2, 2), 0.05, 0.1, 0.2, 0.5, 1))
      kk <- intersect(colnames(sig.up), pvals)
      sig.up <- sig.up[, match(kk, colnames(sig.up)), drop = FALSE]
      sig.down <- sig.down[, match(kk, colnames(sig.down)), drop = FALSE]

      colnames(sig.up)[1] <- paste("UP   FDR = ", colnames(sig.up)[1])
      colnames(sig.down)[1] <- paste("DOWN   FDR = ", colnames(sig.down)[1])
      colnames(sig.down) <- paste0("  ", colnames(sig.down))
      sigcount <- cbind(sig.down, sig.up[rownames(sig.down), , drop = FALSE])
      dim(sigcount)
      maxsig <- 0.99 * max(sigcount, na.rm = TRUE)

      contr <- sub("::.*", "", rownames(sigcount))
      metd <- sub(".*::", "", rownames(sigcount))
      D <- data.frame(method = metd, contrast = contr, sigcount, check.names = FALSE)

      DT::datatable(D,
        rownames = FALSE,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        plugins = "scrollResize",
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          pageLength = 999, #
          scrollX = TRUE,
          scrollY = "20vh",
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(colnames(sig.up),
          background = DT::styleColorBar(c(0, maxsig), "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) %>%
        DT::formatStyle(colnames(sig.down),
          background = DT::styleColorBar(c(0, maxsig), "lightblue"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    FDRtable.RENDER_modal <- shiny::reactive({
      dt <- FDRtable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "datasets",
      func = FDRtable.RENDER,
      func2 = FDRtable.RENDER_modal,
      selector = "none"
    )
  })
}
