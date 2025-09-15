##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_boxplots_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::selectizeInput(
      ns("selected_pheno"), "Select phenotype",
      choices = NULL, multiple = TRUE
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

mofa_plot_boxplots_server <- function(id,
                                      mofa,
                                      input_factor = reactive(1),
                                      watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    observeEvent(mofa(), {
      data <- mofa()
      phenos <- colnames(data$samples)
      updateSelectInput(
        session,
        "selected_pheno",
        choices = phenos,
        selected = head(phenos, 12)
      )
    })

    plot.RENDER <- function() {
      res <- mofa()
      shiny::req(res)
      resF <- data.frame(res$F)
      samples <- data.frame(res$samples)

      k <- input_factor()
      if (!is.null(k)) {
        shiny::req(k %in% colnames(resF))
      }

      pheno <- input$selected_pheno
      shiny::validate(shiny::need(length(pheno) > 0, "Please select at least one phenotype."))

      cm <- intersect(pheno, colnames(samples))
      pheno <- pheno[pheno %in% cm]
      samples <- samples[, cm, drop = FALSE]

      nph <- length(pheno)
      nr <- max(ceiling(sqrt(nph)), 2)
      nc <- ceiling(nph / nr)
      par(mfrow = c(nr, nc), mar = c(3, 4, 2.8, 0.5))

      i <- 1
      for (i in 1:length(pheno)) {
        jj <- which(colnames(samples) == pheno[i])
        if (!any(jj)) next
        y <- samples[, jj]
        if (all(is.na(y))) next
        f1 <- resF[, k]
        c1 <- (class(y) %in% c("numeric", "integer"))
        c2 <- (all(y %in% c(NA, 0:9)))
        if (c1 && c2) y <- factor(y)
        isfactor <- (class(y) %in% c("character", "factor", "logical"))
        ylab <- paste(k, "score")
        if (isfactor) {
          y <- factor(y)
          boxplot(f1 ~ y, main = "", ylab = ylab, xlab = "", las = 1)
        } else {
          plot(y, f1, main = "", ylab = ylab, xlab = "", las = 1)
        }
        title(pheno[i], cex.main = 1.2)
      }
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 10, pdf.height = 6,
      res = c(85, 120),
      add.watermark = watermark
    )
  })
}
