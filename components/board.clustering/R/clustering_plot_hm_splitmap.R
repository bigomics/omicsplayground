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

  topmodes <- c("sd","pca","specific")

  hm_splitmap_opts = shiny::tagList(
    withTooltip( shiny::radioButtons(ns("hm_plottype"), "Plot type:",
                                     choices=c("ComplexHeatmap","iHeatmap"),
                                     selected="ComplexHeatmap", inline=TRUE, width='100%'),
                 "Choose plot type: ComplexHeatmap (static) or iHeatmap (interactive)",
                 placement="right",options = list(container = "body")),
    withTooltip( shiny::radioButtons(
      ns("hm_splitby"), "Split samples by:", inline=TRUE,
      ## selected="phenotype",
      choices=c("none","phenotype","gene")),
      "Split the samples by phenotype or expression level of a gene.",
      placement="right",options = list(container = "body")),
    shiny::conditionalPanel(
      "input.hm_splitby != 'none'", ns=ns,
      withTooltip( shiny::selectInput(ns("hm_splitvar"), NULL, choices=""),
                   "Specify phenotype or gene for splitting the columns of the heatmap.",
                   placement="right",options = list(container = "body")),
    ),
    shiny::fillRow(
      height = 50,
      withTooltip( shiny::selectInput(ns('hm_topmode'),'Top mode:',topmodes, width='100%'),
                   "Specify the criteria for selecting top features to be shown in the heatmap.",
                   placement = "right", options = list(container = "body")),
      withTooltip( shiny::selectInput(ns('hm_ntop'),'Top N:',c(50,150,500),selected=50),
                   "Select the number of top features in the heatmap.",
                   placement="right", options = list(container = "body")),
      withTooltip( shiny::selectInput(ns('hm_clustk'),'K:',1:6,selected=4),
                   "Select the number of gene clusters.",
                   placement="right", options = list(container = "body"))
    ),
    ##br(),
    withTooltip( shiny::radioButtons(
      ns('hm_scale'), 'Scale:', choices=c('relative','absolute','BMC'), inline=TRUE),
      ## ns('hm_scale'), 'Scale:', choices=c('relative','absolute'), inline=TRUE),
      "Show relative (i.e. mean-centered), absolute expression values or batch-mean-centered.",
      placement="right", options = list(container = "body")),
    withTooltip( shiny::checkboxInput(
      ns('hm_legend'), 'show legend', value=TRUE), "Show or hide the legend.",
      placement="right", options = list(container = "body")),
    shiny::fillRow(
      height = 50,
      ## shiny::checkboxInput(ns("hm_labRow"),NULL),
      withTooltip( shiny::numericInput(ns("hm_cexRow"), "cexRow:", 1, 0, 1.4, 0.1, width='100%'),
                   "Specify the row label size. Set to 0 to suppress row labels.",
                   placement="right",options = list(container = "body")),
      withTooltip( shiny::numericInput(ns("hm_cexCol"), "cexCol:", 1, 0, 1.4, 0.1, width='100%'),
                   "Specify the column label size. Set to 0 to suppress column labels.",
                   placement="right", options = list(container = "body"))
    ),
    shiny::br()
  )

  info_text <- "Under the <strong>Heatmap</strong> panel, hierarchical clustering can be performed on gene level or gene set level expression in which users have to specify it under the {Level} dropdown list. <p>Under the plot configuration {{Settings}}, users can split the samples by a phenotype class (e.g., tissue, cell type, or gender) using the {split by} setting. In addition, users can specify the top N = (50, 150, 500) features to be used in the heatmap. The ordering of top features is selected under {top mode}. The criteria to select the top features are: <ol><li>SD - features with the highest standard deviation across all the samples, </li><li>specific - features that are overexpressed in each phenotype class compared to the rest, or by </li><li>PCA - by principal components.<br></ol> <br><p>Users can also choose between 'relative' or 'absolute' expression scale. Under the {cexCol} and {cexRow} settings, it is also possible to adjust the cex for the column and row labels."

  PlotModuleUI(
    ns("pltmod"),
    title = "Clustered Heatmap",
    label = label,
    plotlib = "generic", #generic
    outputFunc = plotly::plotlyOutput, #"uiOutput"
    outputFunc2 = plotly::plotlyOutput, #"uiOutput",
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
                                               getTopMatrix,
                                               watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    fullH <- 850

    ns <- session$ns

    # reactive function listening for changes in input
    topmodes <- c("sd","pca","specific")

    plot_data_hm1 <- shiny::reactive({

      ## ComplexHeatmap based splitted heatmap ##########

      filt <- getTopMatrix()
      shiny::req(filt)
      ##if(is.null(filt)) return(NULL)

      ##if(input$hm_group) {
      zx <- filt$mat
      annot = filt$annot
      zx.idx <- filt$idx

      return(list(
        zx = zx,
        annot = annot,
        zx.idx = zx.idx
      ))
    })

    hm1_splitmap.RENDER<- function() {
      pd <- plot_data_hm1()

      zx = pd[["zx"]]
      annot = pd[["annot"]]
      zx.idx = pd[["zx.idx"]]


      if(nrow(zx) <= 1) return(NULL)

      show_rownames = TRUE
      if(nrow(zx) > 100) show_rownames = FALSE

      cex1 = ifelse(ncol(zx)>50,0.75,1)
      cex1 = ifelse(ncol(zx)>100,0.5,cex1)
      cex1 = ifelse(ncol(zx)>200,0,cex1)

      scale.mode = "none"
      if(input$hm_scale=="relative") scale.mode <- "row.center"
      if(input$hm_scale=="BMC") scale.mode <- "row.bmc"
      scale.mode

      ## split genes dimension in 5 groups
      splity = 5
      splity = 6
      if(!is.null(zx.idx)) splity = zx.idx

      ## split samples
      splitx = NULL
      splitx = filt$grp

      show_legend=show_colnames=TRUE
      show_legend <- input$hm_legend
      if(input$hm_level=="geneset" || !is.null(splitx)) show_legend = FALSE

      annot$group = NULL  ## no group in annotation??
      show_colnames <- (input$hm_cexCol != 0)
      ##if(ncol(zx) > 200) show_colnames <- FALSE ## never...

      if(input$hm_level=="gene") {
        ## strip any prefix
        rownames(zx) = sub(".*:","",rownames(zx))
      }
      rownames(zx) <- sub("HALLMARK:HALLMARK_","HALLMARK:",rownames(zx))
      rownames(zx) = gsub(GSET.PREFIX.REGEX,"",rownames(zx))
      rownames(zx) = substring(rownames(zx),1,50)  ## cut long names...
      if(input$hm_level=="geneset")  rownames(zx) <- tolower(rownames(zx))

      cex2 <- ifelse( nrow(zx) > 60, 0.8, 0.9)
      cex1 <- as.numeric(input$hm_cexCol)*0.85
      cex2 <- as.numeric(input$hm_cexRow)*0.75
      cex0 <- ifelse(!is.null(splitx) && length(splitx)<=10, 1.05, 0.85)  ## title

      crot <- 0
      totnchar <- nchar(paste0(unique(splitx),collapse=""))
      totnchar
      nx <- length(unique(splitx))
      if(!is.null(splitx) & (totnchar > 44 || nx>=6) ) crot=90

      nrownames = 60
      nrownames = 9999
      if(input$hm_cexRow==0) nrownames <- 0

      shiny::showNotification('rendering heatmap...')
      plt <- grid::grid.grabExpr(
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
          key.offset = c(0.89,1.01),
          main=" ", nmax = -1, mar = c(8,16)
        )
      )
      plt

    }

    hm2_splitmap.RENDER<- function() {

      ## iHeatmap based splitted heatmap #########

      shiny::req(pgx$genes)

      ## -------------- variable to split samples
      ##scale = ifelse(input$hm_scale=="relative","row.center","none")
      scale = "none"
      if(input$hm_scale=="relative") scale <- "row.center"
      if(input$hm_scale=="BMC") scale <- "row.bmc"

      plt <- NULL

      filt <- getTopMatrix()
      ##if(is.null(filt)) return(NULL)
      shiny::req(filt)

      ##if(input$hm_group) {
      X <- filt$mat
      annot = filt$annot
      idx <- filt$idx

      ## sample clustering index
      splitx <- NULL
      splitx <- filt$grp

      ## iheatmapr needs factors for sharing between groups
      annotF <- data.frame(as.list(annot),stringsAsFactors=TRUE)
      rownames(annotF) = rownames(annot)

      colcex <- as.numeric(input$hm_cexCol)
      rowcex = as.numeric(input$hm_cexRow)

      tooltips = NULL
      if(input$hm_level=="gene") {
        getInfo <- function(g) {
          aa = paste0("<b>",pgx$genes[g,"gene_name"],"</b>. ",
                      ## pgx$genes[g,"map"],". ",
                      pgx$genes[g,"gene_title"],".")
          breakstring2(aa, 50, brk="<br>")
        }
        tooltips = sapply(rownames(X), getInfo)
      } else {
        aa = gsub("_"," ",rownames(X)) ## just geneset names
        tooltips = breakstring2(aa, 50, brk="<br>")
      }
      ##genetips = rownames(X)

      shiny::showNotification('rendering iHeatmap...')

      plt <- pgx.splitHeatmapFromMatrix(
        X=X, annot=annotF, ytips=tooltips,
        idx=idx, splitx=splitx, scale=scale,
        row_annot_width=0.03, rowcex=rowcex,
        colcex=colcex )

      return(plt)

    }


    output$hm1_splitmap <- shiny::renderPlot({
      plt <- hm1_splitmap.RENDER()
      grid::grid.draw(plt, recording=FALSE)
    }, res=90)

    output$hm2_splitmap <- renderIheatmap({
      hm2_splitmap.RENDER()
    })

    hm_splitmap.switchRENDER <- shiny::reactive({
      ##req(input$hm_plottype)
      p = NULL
      if(input$hm_plottype %in% c("ComplexHeatmap","static") ) {
        p = shiny::plotOutput(ns("hm1_splitmap"), height=fullH-80)  ## height defined here!!
      } else {
        p = iheatmaprOutput(ns("hm2_splitmap"), height=fullH-80) ## height defined here!!
      }
      return(p)
    })

    ##output$hm_splitmap_pdf <- shiny::downloadHandler(
    hm_splitmap_downloadPDF <- shiny::downloadHandler(
      filename = "plot.pdf",
      content = function(file) {
        ##PDFFILE = hm_splitmap_module$.tmpfile["pdf"]  ## from above!
        PDFFILE = paste0(gsub("file","plot",tempfile()),".pdf")
        dbg("[ClusteringBoard] hm_splitmap_downloadPDF: exporting SWITCH to PDF...")
        ##showNotification("exporting to PDF")
        ##wd <- input$hm_pdfwidth
        ##ht <- input$hm_pdfheight
        ##wd <- input$pdf_width
        ##ht <- input$pdf_height
        wd <- input[["hm_splitmap-pdf_width"]]  ## ugly!!
        ht <- input[["hm_splitmap-pdf_height"]] ## ugly!!

        if(1 && input$hm_plottype %in% c("ComplexHeatmap","static")) {
          pdf(PDFFILE, width=wd, height=ht)
          grid::grid.draw(hm1_splitmap.RENDER())
          ##print(hm1_splitmap.RENDER())
          ##hm1_splitmap.RENDER()
          dev.off()
        } else {
          save_iheatmap(hm2_splitmap.RENDER(), filename=PDFFILE,
                        vwidth=wd*100, vheight=ht*100)
        }
        if(WATERMARK) {
          dbg("[ClusteringBoard] adding watermark to PDF...\n")
          addWatermark.PDF(PDFFILE) ## from pgx-modules.R
        }
        dbg("[ClusteringBoard] hm_splitmap_downloadPDF: exporting done...")
        file.copy(PDFFILE,file)
      }
    )

    hm_splitmap_downloadPNG <- shiny::downloadHandler(
      filename = "plot.png",
      content = function(file) {
        PNGFILE = paste0(gsub("file","plot",tempfile()),".png")
        dbg("[ClusteringBoard] hm_splitmap_downloadPDF:: exporting SWITCH to PNG...")
        ##showNotification("exporting to PNG")
        wd <- 100*as.integer(input[["hm_splitmap-pdf_width"]])
        ht <- 100*as.integer(input[["hm_splitmap-pdf_height"]])
        if(1 && input$hm_plottype %in% c("ComplexHeatmap","static")) {
          png(PNGFILE, width=wd, height=ht, pointsize=24)
          grid::grid.draw(hm1_splitmap.RENDER())
          ##print(hm1_splitmap.RENDER())  ## should be done inside render for base plot...
          ##hm1_splitmap.RENDER()  ## should be done inside render for base plot...
          ##plot(sin)
          dev.off()
        } else {
          save_iheatmap(hm2_splitmap.RENDER(), filename=PNGFILE,  vwidth=wd, vheight=ht)
        }
        dbg("[ClusteringBoard] hm_splitmap_downloadPNG: exporting done...")
        file.copy(PNGFILE,file)
      }
    )

    # hm_splitmap_downloadHTML <- shiny::downloadHandler(
    #   filename = "plot.html",
    #   content = function(file) {
    #     ##HTMLFILE = hm_splitmap_module$.tmpfile["html"]  ## from above!
    #     HTMLFILE = paste0(gsub("file","plot",tempfile()),".html")
    #     dbg("renderIheatmap:: exporting SWITCH to HTML...")
    #     shiny::withProgress({
    #       ##write("<body>HTML export error</body>", file=HTMLFILE)
    #       p <- hm2_splitmap.RENDER()
    #       shiny::incProgress(0.5)
    #       save_iheatmap(p, filename=HTMLFILE)
    #     }, message="exporting to HTML", value=0 )
    #     dbg("renderIheatmap:: ... exporting done")
    #     file.copy(HTMLFILE,file)
    #   }
    # )

    PlotModuleServer(
      "pltmod",
      plotlib = "generic",
      func = hm_splitmap.switchRENDER,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI,
      func2 = hm_splitmap.switchRENDER,
      res = c(80, 95), ## resolution of plots
      pdf.width = 10, pdf.height = 8,
      download.pdf = hm_splitmap_downloadPDF,
      download.png = hm_splitmap_downloadPNG,
      add.watermark = watermark

    )


    # ## call plotModule
    # hm_splitmap_module <- shiny::callModule(
    #   plotModule,
    #   id = "hm_splitmap",
    #   func = hm_splitmap.switchRENDER, ## ns=ns,
    #   ## func2 = hm_splitmap.switchRENDER, ## ns=ns,
    #   show.maximize = FALSE,
    #   plotlib = "generic",
    #   renderFunc = "renderUI",
    #   outputFunc = "uiOutput",
    #   download.fmt = c("pdf","png"),
    #   options = hm_splitmap_opts,
    #   height = fullH-80, ##???
    #   width = '100%',
    #   pdf.width = 10, pdf.height = 8,
    #   title ="Clustered Heatmap",
    #   info.text = hm_splitmap_text,
    #   info.width = "350px",
    #   ## caption = hm_splitmap_caption,
    #   download.pdf = hm_splitmap_downloadPDF,
    #   download.png = hm_splitmap_downloadPNG,
    #   download.html = hm_splitmap_downloadHTML,
    # #   add.watermark = WATERMARK
    # )
    #


  }) ## end of moduleServer
}
