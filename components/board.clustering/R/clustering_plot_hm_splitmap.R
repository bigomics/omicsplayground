##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' Clustering plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
clustering_plot_hm_splitmap_ui <- function(id,
                                           label = "",
                                           height,
                                           width) {
  ns <- shiny::NS(id)

  topmodes <- c("sd", "pca", "specific")

  hm_splitmap_opts <- shiny::tagList(
    # withTooltip( shiny::radioButtons(ns("hm_plottype"), "Plot type:",
    #                                  choices=c("ComplexHeatmap","iHeatmap"),
    #                                  selected="ComplexHeatmap", inline=TRUE, width='100%'),
    #              "Choose plot type: ComplexHeatmap (static) or iHeatmap (interactive)",
    #              placement="right",options = list(container = "body")),
    withTooltip(
      shiny::radioButtons(
        ns("hm_splitby"), "Split samples by:",
        inline = TRUE,
        ## selected="phenotype",
        choices = c("none", "phenotype", "gene")
      ),
      "Split the samples by phenotype or expression level of a gene.",
      placement = "right", options = list(container = "body")
    ),
    shiny::conditionalPanel(
      "input.hm_splitby != 'none'",
      ns = ns,
      withTooltip(shiny::selectInput(ns("hm_splitvar"), NULL, choices = ""),
        "Specify phenotype or gene for splitting the columns of the heatmap.",
        placement = "right", options = list(container = "body")
      ),
    ),
    shiny::fillRow(
      height = 50,
      withTooltip(shiny::selectInput(ns("hm_topmode"), "Top mode:", topmodes, width = "100%"),
        "Specify the criteria for selecting top features to be shown in the heatmap.",
        placement = "right", options = list(container = "body")
      ),
      withTooltip(shiny::selectInput(ns("hm_ntop"), "Top N:", c(50, 150, 500), selected = 50),
        "Select the number of top features in the heatmap.",
        placement = "right", options = list(container = "body")
      ),
      withTooltip(shiny::selectInput(ns("hm_clustk"), "K:", 1:6, selected = 4),
        "Select the number of gene clusters.",
        placement = "right", options = list(container = "body")
      )
    ),
    ## br(),
    withTooltip(
      shiny::radioButtons(
        ns("hm_scale"), "Scale:",
        choices = c("relative", "absolute", "BMC"), inline = TRUE
      ),
      ## ns('hm_scale'), 'Scale:', choices=c('relative','absolute'), inline=TRUE),
      "Show relative (i.e. mean-centered), absolute expression values or batch-mean-centered.",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("hm_legend"), "show legend",
        value = TRUE
      ), "Show or hide the legend.",
      placement = "right", options = list(container = "body")
    ),
    shiny::fillRow(
      height = 50,
      ## shiny::checkboxInput(ns("hm_labRow"),NULL),
      withTooltip(shiny::numericInput(ns("hm_cexRow"), "cexRow:", 1, 0, 1.4, 0.1, width = "100%"),
        "Specify the row label size. Set to 0 to suppress row labels.",
        placement = "right", options = list(container = "body")
      ),
      withTooltip(shiny::numericInput(ns("hm_cexCol"), "cexCol:", 1, 0, 1.4, 0.1, width = "100%"),
        "Specify the column label size. Set to 0 to suppress column labels.",
        placement = "right", options = list(container = "body")
      )
    ),
    shiny::br()
  )

  info_text <- "<b>Clustered heatmap.</b> Heatmap showing gene expression sorted by 2-way hierarchical clustering. Red corresponds to overexpression, blue to underexpression of the gene.  At the same time, gene clusters are functionally annotated in the 'Annotate clusters' panel on the right. Hierarchical clustering can be performed on gene level or gene set level expression in which users have to specify it under the {Level} dropdown list. Under the plot settings, users can split the samples by a phenotype class (e.g., tissue, cell type, or gender) using the {split by} setting. In addition, users can specify the top N = (50, 150, 500) features to be used in the heatmap. The ordering of top features is selected under {top mode}. The criteria to select the top features are: SD - features with the highest standard deviation across all the samples,specific - features that are overexpressed in each phenotype class compared to the rest, or by PCA - by principal components. Users can also choose between 'relative' or 'absolute' expression scale. Under the {cexCol} and {cexRow} settings, it is also possible to adjust the cex for the column and row labels."

  PlotModuleUI(
    ns("pltmod"),
    title = "Clustered Heatmap",
    label = label,
    plotlib = "plotly",
    info.text = info_text,
    options = hm_splitmap_opts,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
  )
}

#' Clustering plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#' @param watermark
#'
#'
#'
#' @export
clustering_plot_hm_splitmap_server <- function(id,
                                               pgx,
                                               getTopMatrix,
                                               hm_level,
                                               watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    fullH <- 850

    ns <- session$ns

    shiny::observeEvent(input$hm_splitby, {
      shiny::req(pgx$X, pgx$samples)

      if (input$hm_splitby == "none") {
        return()
      }
      if (input$hm_splitby == "gene") {
        xgenes <- sort(rownames(pgx$X))
        shiny::updateSelectizeInput(session, "hm_splitvar", choices = xgenes, server = TRUE)
      }
      if (input$hm_splitby == "phenotype") {
        cvar <- sort(pgx.getCategoricalPhenotypes(pgx$samples, min.ncat = 2, max.ncat = 999))
        sel <- cvar[1]
        cvar0 <- grep("^[.]", cvar, value = TRUE, invert = TRUE) ## no estimated vars
        sel <- head(c(
          grep("type|family|class|stat", cvar0, ignore.case = TRUE, value = TRUE),
          cvar0, cvar
        ), 1)
        shiny::updateSelectInput(session, "hm_splitvar", choices = cvar, selected = sel)
      }
    })

    ## update filter choices upon change of data set
    shiny::observe({
      shiny::updateRadioButtons(session, "hm_splitby", selected = "none")
    })

    plot_data_hm1 <- shiny::reactive({
      ## ComplexHeatmap based splitted heatmap ##########

      filt <- getTopMatrix()
      shiny::req(filt)
      ## if(is.null(filt)) return(NULL)

      ## if(input$hm_group) {
      zx <- filt$mat
      annot <- filt$annot
      zx.idx <- filt$idx

      return(list(
        zx = zx,
        annot = annot,
        zx.idx = zx.idx,
        filt = filt
      ))
    })

    hm1_splitmap.RENDER <- function() {
      pd <- plot_data_hm1()

      zx <- pd[["zx"]]
      annot <- pd[["annot"]]
      zx.idx <- pd[["zx.idx"]]
      filt <- pd[["filt"]]

      if (nrow(zx) <= 1) {
        return(NULL)
      }

      show_rownames <- TRUE
      if (nrow(zx) > 100) show_rownames <- FALSE

      cex1 <- ifelse(ncol(zx) > 50, 0.75, 1)
      cex1 <- ifelse(ncol(zx) > 100, 0.5, cex1)
      cex1 <- ifelse(ncol(zx) > 200, 0, cex1)

      scale.mode <- "none"
      if (input$hm_scale == "relative") scale.mode <- "row.center"
      if (input$hm_scale == "BMC") scale.mode <- "row.bmc"
      scale.mode

      ## split genes dimension in 5 groups
      splity <- 5
      splity <- 6
      if (!is.null(zx.idx)) splity <- zx.idx

      ## split samples
      splitx <- NULL
      splitx <- filt$grp

      show_legend <- show_colnames <- TRUE
      show_legend <- input$hm_legend
      if (hm_level() == "geneset" || !is.null(splitx)) show_legend <- FALSE

      annot$group <- NULL ## no group in annotation??
      show_colnames <- (input$hm_cexCol != 0)
      ## if(ncol(zx) > 200) show_colnames <- FALSE ## never...

      if (hm_level() == "gene") {
        ## strip any prefix
        rownames(zx) <- sub(".*:", "", rownames(zx))
      }
      rownames(zx) <- sub("HALLMARK:HALLMARK_", "HALLMARK:", rownames(zx))
      rownames(zx) <- gsub(GSET.PREFIX.REGEX, "", rownames(zx))
      rownames(zx) <- substring(rownames(zx), 1, 50) ## cut long names...
      if (hm_level() == "geneset") rownames(zx) <- tolower(rownames(zx))

      cex2 <- ifelse(nrow(zx) > 60, 0.8, 0.9)
      cex1 <- as.numeric(input$hm_cexCol) * 0.85
      cex2 <- as.numeric(input$hm_cexRow) * 0.75
      cex0 <- ifelse(!is.null(splitx) && length(splitx) <= 10, 1.05, 0.85) ## title

      crot <- 0
      totnchar <- nchar(paste0(unique(splitx), collapse = ""))
      totnchar
      nx <- length(unique(splitx))
      if (!is.null(splitx) & (totnchar > 44 || nx >= 6)) crot <- 90

      nrownames <- 60
      nrownames <- 9999
      if (input$hm_cexRow == 0) nrownames <- 0

      shiny::showNotification("Rendering heatmap...")
      # plt <- grid::grid.grabExpr(
      gx.splitmap(
        zx,
        split = splity, splitx = splitx,
        scale = scale.mode, show_legend = show_legend,
        show_colnames = show_colnames, column_title_rot = crot,
        column_names_rot = 45,
        show_rownames = nrownames, rownames_width = 40,
        softmax = 0,
        ## side.height.fraction=0.03+0.055*NCOL(annot),
        title_cex = cex0, cexCol = cex1, cexRow = cex2,
        col.annot = annot, row.annot = NULL, annot.ht = 2.3,
        key.offset = c(0.89, 1.01),
        main = " ", nmax = -1, mar = c(8, 16)
      )
      p <- grDevices::recordPlot()
      p
      # )
      # browser()
      # plt
    }

    hm2_splitmap.RENDER <- function() {
      ## iHeatmap based splitted heatmap #########

      shiny::req(pgx$genes)

      ## -------------- variable to split samples
      ## scale = ifelse(input$hm_scale=="relative","row.center","none")
      scale <- "none"
      if (input$hm_scale == "relative") scale <- "row.center"
      if (input$hm_scale == "BMC") scale <- "row.bmc"

      plt <- NULL

      filt <- getTopMatrix()
      ## if(is.null(filt)) return(NULL)
      shiny::req(filt)

      ## if(input$hm_group) {
      X <- filt$mat
      annot <- filt$annot
      idx <- filt$idx

      ## sample clustering index
      splitx <- NULL
      splitx <- filt$grp

      ## iheatmapr needs factors for sharing between groups
      annotF <- data.frame(as.list(annot), stringsAsFactors = TRUE)
      rownames(annotF) <- rownames(annot)

      colcex <- as.numeric(input$hm_cexCol)
      rowcex <- as.numeric(input$hm_cexRow)

      tooltips <- NULL
      if (hm_level() == "gene") {
        getInfo <- function(g) {
          aa <- paste0(
            "<b>", pgx$genes[g, "gene_name"], "</b>. ",
            ## pgx$genes[g,"map"],". ",
            pgx$genes[g, "gene_title"], "."
          )
          breakstring2(aa, 50, brk = "<br>")
        }
        tooltips <- sapply(rownames(X), getInfo)
      } else {
        aa <- gsub("_", " ", rownames(X)) ## just geneset names
        tooltips <- breakstring2(aa, 50, brk = "<br>")
      }
      ## genetips = rownames(X)

      shiny::showNotification("Rendering iHeatmap...")

      plt <- pgx.splitHeatmapFromMatrix(
        X = X, annot = annotF, ytips = tooltips,
        idx = idx, splitx = splitx, scale = scale,
        row_annot_width = 0.03, rowcex = rowcex,
        colcex = colcex
      )
      obj2 <- plt %>% iheatmapr::to_plotly_list()
      plt <- plotly::as_widget(obj2)
      return(plt)
    }


    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = hm2_splitmap.RENDER,
      csvFunc = plot_data_hm1,
      res = c(80, 95), ## resolution of plots
      pdf.width = 10, pdf.height = 8,
      add.watermark = watermark
    )

    return(list(
      hm_ntop = shiny::reactive(input$hm_ntop),
      hm_splitvar = shiny::reactive(input$hm_splitvar),
      hm_splitby = shiny::reactive(input$hm_splitby),
      hm_scale = shiny::reactive(input$hm_scale),
      hm_topmode = shiny::reactive(input$hm_topmode),
      hm_clustk = shiny::reactive(input$hm_clustk)
    ))
  }) ## end of moduleServer
}
