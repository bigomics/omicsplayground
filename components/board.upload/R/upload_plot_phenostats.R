##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

upload_plot_phenostats_ui <- function(id,
                                      label = "",
                                      height,
                                      width,
                                      title,
                                      caption,
                                      info.text) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    plotlib = "base",
    info.text = info.text,
    caption = caption,
    options = NULL,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

plotPhenoDistribution <- function(pheno) {
  if (nrow(pheno) == 0 || NCOL(pheno) == 0) {
    return(NULL)
  }
  px <- head(colnames(pheno), 20) ## show maximum??
  df <- type.convert(pheno[, px, drop = FALSE], as.is = TRUE)
  vt <- df %>% inspectdf::inspect_types()

  ## discretized continuous variable into 10 bins
  ii <- unlist(vt$col_name[c("numeric", "integer")])
  if (!is.null(ii) && length(ii)) {
    cat("[UploadModule::phenoStats] discretizing variables:", ii, "\n")
    df[, ii] <- apply(df[, ii, drop = FALSE], 2, function(x) {
      if (any(is.infinite(x))) x[which(is.infinite(x))] <- NA
      cut(x, breaks = 10)
    })
  }

  p1 <- df %>%
    inspectdf::inspect_cat() %>%
    inspectdf::show_plot()
  tt2 <- paste(nrow(pheno), "samples x", ncol(pheno), "phenotypes")
  #
  p1 <- p1 + ggplot2::ggtitle("PHENOTYPES", subtitle = tt2) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        size = 12,
        margin = ggplot2::margin(0, 0, 0, 25),
        hjust = 1
      )
    )
  return(p1)
}


upload_plot_phenostats_server <- function(id, checkTables, samplesRT, watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## extract data from pgx object
    plot_data <- shiny::reactive({
      pheno <- samplesRT()
      has.pheno <- !is.null(pheno) && NCOL(pheno) > 0
      check <- checkTables()

      status.ok <- check["samples.csv", "status"]
      status.ds <- tolower(check["samples.csv", "description"])
      error.msg <- paste(
        toupper(status.ok), "\nPlease upload 'samples.csv' (Required):",
        status.ds
      )
      shiny::validate(
        shiny::need(
          status.ok == "OK" && has.pheno,
          error.msg
        )
      )
      pheno <- as.data.frame(pheno, check.names = FALSE, drop = FALSE)
      return(pheno)
    })

    plot.RENDER <- function() {
      pheno <- plot_data()
      p1 <- plotPhenoDistribution(pheno)
      return(p1)
    }

    modal_plot.RENDER <- function() {
      plot.RENDER()
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "base",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 90), ## resolution of plots
      pdf.width = 4, pdf.height = 4,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
