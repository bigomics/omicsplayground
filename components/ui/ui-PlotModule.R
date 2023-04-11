##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


if(0) {
    id="plot1"
    info.text="Figure"
    title=""
    options = NULL
    label=""
    caption=""
    caption2=info.text ## header=NULL
    plotlib = "base"
    plotlib2 = plotlib
    ##renderFunc=NULL
    outputFunc=shiny::plotOutput
    ##csvFunc=NULL
    ##renderFunc2=NULL
    outputFunc2=shiny::plotOutput
    no.download = FALSE
    download.fmt=c("png","pdf")
    just.info=FALSE
    info.width="300px"
    show.maximize = TRUE
    height = c(640,800)
    width = c("auto",1400)
    pdf.width = 6
    pdf.height = 6
}

PlotModuleUI <- function(id,
                       info.text="Figure",
                       title="",
                       options = NULL,
                       label="",
                       caption="",
                       caption2=caption, ## header=NULL,
                       plotlib = "base",
                       plotlib2 = plotlib,
                       outputFunc = NULL,
                       outputFunc2 = NULL,
                       no.download = FALSE,
                       download.fmt=c("png","pdf"),
                       just.info=FALSE,
                       info.width="300px",
                       show.maximize = TRUE,
                       height = c(400,800),
                       card_footer_height = "3rem",
                       width = c("auto","100%"),
                       pdf.width = 8,
                       pdf.height = 8
                       )
{
    require(magrittr)
    ns <- shiny::NS(id)

    if(is.null(plotlib2)) plotlib2 <- plotlib
    if(length(height)==1) height <- c(height,800)
    if(length(width)==1)  width  <- c(width,"100%")


    ifnotchar.int <- function(s) suppressWarnings(
      ifelse(!is.na(as.integer(s)), paste0(as.integer(s),"px"), s))
    width.1  <- ifnotchar.int(width[1])
    width.2  <- ifnotchar.int(width[2])
    height.1 <- ifnotchar.int(height[1])
    height.2 <- ifnotchar.int(height[2])

    ## OVERRIDE WIDTH: for fullscreen modal always 100%
    width.2 = "100%"   
    
    getOutputFunc <- function(plotlib)
    {
        FUN <- switch(
            plotlib,
            generic = NULL,
            htmlwidget = NULL,
            plotly = plotly::plotlyOutput,
##          echarts4r = echarts4r::echarts4rOutput,
##          scatterD3 = scatterD3::scatterD3Output,
            pairsD3 = pairsD3::pairsD3Output,
            visnetwork = visNetwork::visNetworkOutput,
            ggplot = shiny::plotOutput,
            grid = shiny::plotOutput,
            iheatmapr = iheatmapr::iheatmaprOutput,
            image = shiny::imageOutput,
            base = shiny::plotOutput,
            shiny::plotOutput
        )
        FUN
    }

    if(is.null(plotlib2))   plotlib2 <- plotlib
    if(is.null(outputFunc))  outputFunc  <- getOutputFunc(plotlib)
    if(is.null(outputFunc2)) outputFunc2 <- getOutputFunc(plotlib2)

    ##--------------------------------------------------------------------------------
    ##------------------------ BUTTONS -----------------------------------------------
    ##--------------------------------------------------------------------------------

    ##if(is.null(inputs) || length(inputs)==0 ) inputs <- ""
    options.button <- ""

    if(!just.info && !is.null(options) && length(options)>0) {
        options.button <- DropdownMenu(
          options,
          size = "xs",
          width = "auto",
          icon = shiny::icon("bars"),
          status = "default"
        )
    }

    dload.csv = dload.pdf = dload.png = dload.html = dload.obj = NULL
    if("pdf" %in% download.fmt)   dload.pdf  <- shiny::downloadButton(ns("pdf"), "PDF")
    if("png" %in% download.fmt)   dload.png  <- shiny::downloadButton(ns("png"), "PNG")
    if("html" %in% download.fmt)  dload.html <- shiny::downloadButton(ns("html"), "HTML")
    if("csv"  %in% download.fmt)  dload.csv  <- shiny::downloadButton(ns("csv"), "CSV")
    if("obj"  %in% download.fmt)  dload.obj  <- shiny::downloadButton(ns("obj"), "obj")


    pdf_size = NULL
    if(TRUE || plotlib!="base") {
        pdf_size <- shiny::tagList(
            shiny::fillRow(
                shiny::numericInput(ns("pdf_width"), "Width", pdf.width, 1, 20, 1, width='95%'),
                shiny::numericInput(ns("pdf_height"), "Height", pdf.height, 1, 20, 1, width='100%')
            ),
            shiny::br(),shiny::br(),shiny::br()
        )
    }

    dload.button <- DropdownMenu(
      div(
        style = "width: 150px;",
        shiny::selectInput(
          inputId = ns("downloadOption"),
          label = "Format",
          choices = download.fmt
        ),
        shiny::conditionalPanel(
          condition = "input.downloadOption == 'pdf' || input.downloadOption == 'png'",
          ns = ns,
          shiny::div(
            pdf_size,
            shiny::br()
          )
        ),
        div(
          shiny::downloadButton(
            outputId = ns("download"),
            label = "Download",
          )
        )
      ),
      size = "xs",
      icon = shiny::icon("download"),
      status = "default"
    )

    if(no.download || length(download.fmt)==0 ) dload.button <- ""

    if(!is.null(label) && label!="") label <- paste0("&nbsp;(",label,")")
    label <- shiny::div(class = "plotmodule-title", shiny::HTML(label))

    zoom.button <- NULL
    if(1 && show.maximize) {
        zoom.button <- modalTrigger(ns("zoombutton"),
            ns("plotPopup"),
            icon("window-maximize"),
            class="btn-circle-xs"
        )
    }

    header <- shiny::fillRow(
        flex = c(1,NA,NA,NA,NA),
        class="plotmodule-header",
        shiny::div(class='plotmodule-title', title=title, title),
        DropdownMenu(
            shiny::div(class='plotmodule-info', shiny::HTML(paste0("<b>", as.character(title),".", "</b>", "&nbsp;", as.character(info.text)))),
            width = "250px",
            size = "xs",
            icon = shiny::icon("info"),
            status = "default"
        ),
        options.button,
        shiny::div(class='download-button', title='download', dload.button),
        shiny::div(class='zoom-button', title='zoom', zoom.button)
    )

    ## ------------------------------------------------------------------------
    ## --------------- modal UI (former output$popupfig) ----------------------
    ## ------------------------------------------------------------------------

    popupfigUI <- function() {
        w <- width.2
        h <- height.2

        ## NOTE: this was in the server before and we could ask the
        ## image size. How to do this in the UI part?
        if(FALSE && plotlib2=="image") {
            ## retains aspect ratio
            ##
            img.file <- func()$src
            img.dim <- c(h,w)
            if(grepl("png|PNG",img.file)) img.dim <- dim(png::readPNG(img.file))[1:2]
            if(grepl("jpg|JPG",img.file)) img.dim <- dim(jpeg::readJPEG(img.file))[1:2]
            r <- min( width.2 / img.dim[2], height.2 / img.dim[1])
            h <- img.dim[1]*r
            w <- img.dim[2]*r
        }

        ## render caption2 (for modal)
        if(any(class(caption2)=="reactive")) {
            caption2 <- caption2()
        }
        caption2 <- shiny::div(
          class="caption2 popup-plot-caption",
          shiny::HTML(paste0("<b>", as.character(title),".</b>&nbsp;&nbsp;",
            as.character(caption2)))
        )

        ##        shiny::tagList(
        shiny::div(
          class = "popup-plot-body",
          shiny::div(
            class = "popup-plot",
            tryCatch({
            outputFunc2(ns("renderpopup"), width=w, height=h, inline=FALSE)
            }, error = function(x){ # This try-catch is used for render functions that
                                  # do not have the arguments, such as
                                  # iheatmap plots
            outputFunc2(ns("renderpopup"), width=w, height=h)
            })
          ),
          caption2
        )
    }

    popupfigUI_editor <- function(){
        htmlOutput(ns("editor_frame"))
    }

    ## inline styles (should be in CSS...)
    modaldialog.style <- paste0("#",ns("plotPopup")," .modal-dialog {width:",width.2,";}")
    modalbody.style <- paste0("#",ns("plotPopup")," .modal-body {min-height:",height.2,"; padding:30px 150px;}")
    modalcontent.style <- paste0("#",ns("plotPopup")," .modal-content {width:100vw;}")
    modalfooter.none <- paste0("#",ns("plotPopup")," .modal-footer{display:none;}")

    if(any(class(caption)=="reactive")) {
        caption <- caption()
    }
#    if(any(class(caption)=="character")) {
#        caption <- shiny::HTML(caption)
#        caption <- shiny::span(caption)
#    }
    e <- bslib::card(
      full_screen = FALSE, #full_screen = TRUE breaks reactivity
       style = "overflow: visible;",
       bslib::as.card_item(
        div(header)
       ),
      bslib::card_body_fill( #TODO card_body_fill will be deprecated soon, switch to card_body after dev bslib install
        style = paste0("height: ",height.1,";"),
        outputFunc(ns("renderfigure")), #  %>% shinycssloaders::withSpinner(),
        shiny::div(class="popup-modal",
                    modalUI(
                          id = ns("plotPopup"),
                          title = title,
                          size = "fullscreen",
                          footer = NULL,
                          popupfigUI()
                      )
                    ),
        shiny::div(class="popup-modal",
                    modalUI(
                          id = ns("plotPopup_editor"),
                          title = "Editor",
                          size = "fullscreen",
                          footer = NULL,
                          popupfigUI_editor()
                      )
                    ),
        shiny::tagList(
                    shiny::tags$head(shiny::tags$style(modaldialog.style)),
                    shiny::tags$head(shiny::tags$style(modalbody.style)),
                    shiny::tags$head(shiny::tags$style(modalcontent.style)),
                    shiny::tags$head(shiny::tags$style(modalfooter.none))
                    )
      ),
      bslib::card_body( #TODO probably want to set fillable and fill to FALSE 
        class = "card-footer", # center the content horizontally and vertically
        style = paste0("height:", card_footer_height, ";"), # add left and top margin of 2 pixels
         div(
          class = "caption",
          shiny::HTML(paste0("<b>", as.character(title),".</b>&nbsp;",
            as.character(caption)))
         )
      )
    ) # end of card
    ## e <- htmltools::bindFillRole(e, container = FALSE, item = FALSE, overwrite = TRUE)
    return(e)
}


PlotModuleServer <- function(
         id,
         func,
         func2=NULL,
         plotlib = "base",
         plotlib2 = plotlib,
         renderFunc = NULL,
         renderFunc2 = NULL,
         csvFunc=NULL,
         download.fmt=c("png","pdf"),
         height = c(640,800),
         width = c("auto",1400),
         res=c(100,170),
         download.pdf = NULL,
         download.png = NULL,
         download.html = NULL,
         download.csv = NULL,
         download.obj = NULL,
         pdf.width=8,
         pdf.height=6,
         pdf.pointsize=12,
         add.watermark=FALSE )
{
    moduleServer(
      id,
      function(input, output, session) {

          ns <- session$ns

          ##--------------------------------------------------------------------------------
          ##------------------------ Plotly editor------------------------------------------
          ##--------------------------------------------------------------------------------
          octocat <- list(
            name = "Editor",
            icon = list(
              path = "M410.052,46.398c-0.812-10.885-5.509-21.129-13.226-28.845c-16.089-16.089-41.044-17.965-59.34-4.459l-7.427,5.487C257.281,72.291,191.872,135.46,135.647,206.336c-14.115,17.795-27.792,36.287-40.715,55.015c-0.928-0.042-1.859-0.068-2.795-0.068c-16.279,0-31.583,6.339-43.094,17.851C28.607,299.57,27.77,319.906,26.96,339.572c-0.745,18.1-1.449,35.196-16.99,54.271L0,406.081h15.785c37.145,0,96.119-17.431,119.447-40.759c11.511-11.511,17.85-26.815,17.85-43.094c0-0.941-0.026-1.877-0.068-2.81c18.766-12.941,37.258-26.614,55.01-40.704C278.873,222.52,342.046,157.11,395.787,84.302l5.479-7.419C407.747,68.111,410.867,57.284,410.052,46.398z M124.625,354.715c-16.334,16.334-58.786,31.89-94.095,35.555c10.098-18.012,10.791-34.866,11.417-50.082c0.754-18.326,1.406-34.152,17.702-50.449c8.678-8.678,20.216-13.457,32.488-13.457s23.81,4.779,32.488,13.457s13.457,20.215,13.457,32.487C138.082,334.5,133.303,346.037,124.625,354.715z M135.232,279.133c-6.875-6.875-15.11-11.889-24.091-14.825c10.801-15.429,22.107-30.656,33.724-45.426c12.79,1.717,24.7,7.567,33.871,16.737c9.174,9.174,15.027,21.087,16.745,33.875c-14.743,11.601-29.97,22.905-45.427,33.719C147.116,294.236,142.104,286.006,135.232,279.133z M389.2,67.971l-5.481,7.421c-50.415,68.302-109.268,129.976-175.037,183.518c-3.279-12.747-9.915-24.473-19.34-33.897c-9.421-9.421-21.145-16.055-33.893-19.333C209.017,139.887,270.692,81.036,338.97,30.649l7.427-5.488c12.275-9.062,29.023-7.801,39.823,3c5.177,5.177,8.329,12.05,8.874,19.355C395.641,54.822,393.548,62.086,389.2,67.971z",
              transform = 'scale(0.035)'
            ),
            click = htmlwidgets::JS(paste0(
              "function(gd){$('#", ns("plotPopup_editor"), "').modal('show')}"
            )
            )
          )

          getEditorUrl <- function(session, path_object) {
            cd <- session$clientData
            sprintf(
              "%s/?plotURL=%s//%s:%s%s%s",
              "editor/index.html",
              cd$url_protocol,
              cd$url_hostname,
              cd$url_port,
              cd$url_pathname,
              path_object
            )
          }

          output$editor_frame <- renderUI({
            plot <- func()
            json <- plotly::plotly_json(plot, TRUE) # requires `listviewer` to work properly
            res <- session$registerDataObj(
              "plotly_graph", json$x$data,
              function(data, req) {
                httpResponse(
                  status = 200,
                  content_type = 'application/json',
                  content = data
                )
              }
            )
            url <- getEditorUrl(session, res)
            tags$iframe(src=url, style = "height: 85vh; width: 100%;")
          })

          ##--------------------------------------------------------------------------------
          ##------------------------ FIGURE ------------------------------------------------
          ##--------------------------------------------------------------------------------

          ## these engines cannot (yet) provide html
          if(plotlib %in% c("base")) {
              download.fmt <- setdiff(download.fmt, c("html"))
          }

          do.pdf = "pdf" %in% download.fmt
          do.png = "png" %in% download.fmt
          do.html = "html" %in% download.fmt
          do.obj = "obj" %in% download.fmt

          ##do.csv  = "csv" %in% download.fmt && !is.null(csvFunc)
          do.csv = !is.null(csvFunc)

          PNGFILE=PDFFILE=HTMLFILE=CSVFILE=NULL
          if(do.pdf) PDFFILE = paste0(gsub("file","plot",tempfile()),".pdf")
          if(do.png) PNGFILE = paste0(gsub("file","plot",tempfile()),".png")
          if(do.csv) CSVFILE = paste0(gsub("file","data",tempfile()),".csv")
          HTMLFILE = paste0(gsub("file","plot",tempfile()),".html")  ## tempory for webshot
          HTMLFILE
          unlink(HTMLFILE)

          ##============================================================
          ##=============== Download Handlers ==========================
          ##============================================================
          ## download.pdf = NULL
          ##download.png = download.html = NULL

          if(do.png && is.null(download.png)) {
            download.png <- shiny::downloadHandler(
              filename = "plot.png",
              content = function(file) {

                pdf.width  <- input$pdf_width
                pdf.height <- input$pdf_height
                resx <- 4  ## upresolution

                shiny::withProgress({
                  ## unlink(PNGFILE) ## do not remove!
                  if(plotlib=="plotly") {
                    p <- func()
                    p$width = pdf.width * 80
                    p$height = pdf.height * 80
                    plotlyExport(p, PNGFILE, width=p$width, height=p$height)
                  } else if(plotlib=="iheatmapr") {
                    p <- func()
                    iheatmapr::save_iheatmap(p, vwidth=pdf.width*80,vheight=pdf.height*80,PNGFILE)
                  } else if(plotlib=="visnetwork") {
                    p <- func()
                    dbg("[plotModule] visnetwork download PNG : visSave : HTMLFILE=",HTMLFILE)
                    visNetwork::visSave(p, HTMLFILE)
                    dbg("[plotModule] visnetwork download PNG : webshot : PNGFILE = ",PNGFILE)
                    webshot::webshot(url=HTMLFILE,file=PNGFILE,vwidth=pdf.width*100,vheight=pdf.height*100)
                  } else if(plotlib %in% c("htmlwidget","pairsD3","scatterD3")) {
                    p <- func()
                    htmlwidgets::saveWidget(p, HTMLFILE)
                    webshot::webshot(url=HTMLFILE,file=PNGFILE,vwidth=pdf.width*100,vheight=pdf.height*100)
                  } else if(plotlib %in% c("ggplot","ggplot2")) {
                    ggsave(PNGFILE, plot = func(), dpi=300)
                  } else if(plotlib=="grid") {
                    p <- func()
                    png(PNGFILE, width=pdf.width*100*resx, height=pdf.height*100*resx,
                      pointsize=1.2*pdf.pointsize, res=72*resx)
                    grid::grid.draw(p)
                    dev.off()
                  } else if(plotlib=="image") {
                    p <- func()
                    dbg("[downloadHandler.PNG] copy image ",p$src,"to PNGFILE",PNGFILE)
                    file.copy(p$src, PNGFILE, overwrite=TRUE)
                  } else if(plotlib=="generic") {
                    ## generic function should produce PNG inside plot func()
                    ##
                  } else if(plotlib=="base") {
                    png(PNGFILE, width=pdf.width*100*resx, height=pdf.height*100*resx,
                      pointsize=1.2*pdf.pointsize, res=72*resx)
                    func()
                    dev.off()  ## important!!
                  } else { ## end base
                    png(PNGFILE, pointsize=pdf.pointsize)
                    plot.new()
                    mtext("Error. PNG not available.",line=-8)
                    dev.off()
                  }

                  ## finally copy to final exported file
                  dbg("[downloadHandler.PNG] copy PNGFILE",PNGFILE,"to download file",file )
                  file.copy(PNGFILE, file, overwrite=TRUE)
                  ## ImageMagick or pdftk
                  if(TRUE && add.watermark) {
                    message("[plotModule] adding watermark to PNG...")
                    addWatermark.PNG(file)
                  }
                }, message="Exporting to PNG", value=0.8)
              } ## content
            ) ## PNG downloadHandler
          } ## end if do.png

          if(do.pdf && is.null(download.pdf) ) {
              download.pdf <- shiny::downloadHandler(
                filename = "plot.pdf",
                content = function(file) {
                  pdf.width  <- input$pdf_width
                  pdf.height <- input$pdf_height
                  shiny::withProgress({
                    ## unlink(PDFFILE) ## do not remove!
                    if(plotlib=="plotly") {
                      p <- func()
                      p$width = pdf.width * 80
                      p$height = pdf.height * 80
                      ##err <- try(plotly::export(p, PDFFILE))  ## deprecated
                      ##err <- try(plotly::orca(p, PDFFILE))
                      ##err <- try(ORCA$export(p, PDFFILE, width=p$width, height=p$height))
                      plotlyExport(p, PDFFILE, width=p$width, height=p$height)
                    } else if(plotlib=="iheatmapr") {
                      p <- func()
                      iheatmapr::save_iheatmap(p, vwidth=pdf.width*80,vheight=pdf.height*80,PDFFILE)
                    } else if(plotlib=="visnetwork") {
                      p <- func()
                      dbg("[plotModule] visnetwork :: download PDF : visSave : HTMLFILE=",HTMLFILE)
                      visNetwork::visSave(p, HTMLFILE)
                      dbg("[plotModule] visnetwork :: download PDF : webshot ; PDFFILE=",PDFFILE)
                      webshot::webshot(url=HTMLFILE,file=PDFFILE,vwidth=pdf.width*100,vheight=pdf.height*100)
                    } else if(plotlib %in% c("htmlwidget","pairsD3","scatterD3")) {
                      p <- func()
                      htmlwidgets::saveWidget(p, HTMLFILE)
                      webshot::webshot(url=HTMLFILE, file=PDFFILE, vwidth=pdf.width*100,vheight=pdf.height*100)
                    } else if(plotlib %in% c("ggplot","ggplot2")) {
                      p <- func()
                      pdf(PDFFILE, width=pdf.width, height=pdf.height, pointsize=pdf.pointsize)
                      print(p)
                      dev.off()
                    } else if(plotlib %in% c("grid")) {
                      p <- func()
                      pdf(PDFFILE, width=pdf.width, height=pdf.height, pointsize=pdf.pointsize)
                      grid::grid.draw(p)
                      dev.off()
                    } else if(plotlib=="image") {
                      p <- func()
                      ## p$src  ## PNG image file
                      ## generic function should produce PDF inside plot func()
                      ##
                    } else if(plotlib=="generic") {
                      ## generic function should produce PDF inside plot func()
                      ##
                    } else if(plotlib=="base") {
                      pdf(file=PDFFILE, width=pdf.width, height=pdf.height,
                        pointsize=pdf.pointsize)
                      func()
                      dev.off()  ## important!!
                    } else { ## end base
                      pdf(PDFFILE, pointsize=pdf.pointsize)
                      plot.new()
                      mtext("Error. PDF not available.",line=-8)
                      dev.off()
                    }

                    ## finally copy to final exported file
                    dbg("[downloadHandler.PDF] copy PDFFILE",PDFFILE,"to download file",file )
                    file.copy(PDFFILE, file, overwrite=TRUE)

                    ## ImageMagick or pdftk
                    if(TRUE && add.watermark) {
                      message("[plotModule] adding watermark to PDF...")
                      addWatermark.PDF(file)
                    }
                  }, message="Exporting to PDF", value=0.8)
                                         } ## content
              ) ## PDF downloadHandler
          } ## end if do.pdf

          saveHTML <- function() {
              ## unlink(HTMLFILE) ## do not remove!
              if(plotlib == "plotly" ) {
                  p <- func()
                  htmlwidgets::saveWidget(p, HTMLFILE)
              } else if(plotlib %in% c("htmlwidget","pairsD3","scatterD3") ) {
                  p <- func()
                  htmlwidgets::saveWidget(p, HTMLFILE)
              } else if(plotlib == "iheatmapr") {
                  p <- func()
                  iheatmapr::save_iheatmap(p, HTMLFILE)
              } else if(plotlib == "visnetwork") {
                  p <- func()
                  visNetwork::visSave(p, HTMLFILE)
              } else if(plotlib %in% c("ggplot","ggplot2")) {
                  p <- func()
                  ##ggsave(PDFFILE, width=pdf.width, height=pdf.height)
                  htmlwidgets::saveWidget( plotly::ggplotly(p), file = HTMLFILE);
              } else if(plotlib=="image") {
                  write("<body>image cannot export to HTML</body>",HTMLFILE)
              } else if(plotlib=="generic") {
                  ## generic function should produce PDF inside plot func()
                  ##
              } else if(plotlib=="base") {
                  write("<body>R base plots cannot export to HTML</body>",HTMLFILE)
              } else { ## end base
                  write("<body>HTML export error</body>",file=HTMLFILE)
              }
              return(HTMLFILE)
          }

          if(do.html && is.null(download.html) )  {
            download.html <- shiny::downloadHandler(
              filename = "plot.html",
              content = function(file) {
                shiny::withProgress({
                  ## unlink(HTMLFILE) ## do not remove!
                  if(plotlib == "plotly") {
                    p <- func()
                    htmlwidgets::saveWidget(p, HTMLFILE)
                  } else if(plotlib %in% c("htmlwidget","pairsD3","scatterD3") ) {
                    p <- func()
                    htmlwidgets::saveWidget(p, HTMLFILE)
                  } else if(plotlib == "iheatmapr") {
                    p <- func()
                    iheatmapr::save_iheatmap(p, HTMLFILE)
                  } else if(plotlib == "visnetwork") {
                    p <- func()
                    visNetwork::visSave(p, HTMLFILE)
                  } else if(plotlib %in% c("ggplot","ggplot2")) {
                    p <- func()
                    ##ggsave(PDFFILE, width=pdf.width, height=pdf.height)
                    htmlwidgets::saveWidget( plotly::ggplotly(p), file = HTMLFILE);
                  } else if(plotlib=="generic") {
                    ## generic function should produce PDF inside plot func()
                    ##
                  } else if(plotlib=="image") {
                    write("<body>image cannot be exported to HTML</body>",HTMLFILE)
                  } else if(plotlib=="base") {
                    write("<body>R base plots cannot be exported to HTML</body>",HTMLFILE)
                  } else { ## end base
                    write("<body>HTML export error</body>",file=HTMLFILE)
                  }
                  ## finally copy to fina lexport file
                  file.copy(HTMLFILE, file, overwrite=TRUE)
                }, message="Exporting to HTML", value=0.8)
              } ## end of content
            ) ## end of HTML downloadHandler
          } ## end of do HTML

          if(do.obj)  {
            if(plotlib == "plotly") {
              download.obj <- shiny::downloadHandler(
                filename = "plot.rds",
                content = function(file) {
                  shiny::withProgress({
                    p <- func()
                    ## we need to strip away unnecessary environment to prevent save bloat
                    b <- plotly::plotly_build(p)$x[c("data", "layout", "config")]
                    #b <- plotly_build(p); $x$attr <- NULL; b$x$visdat <- NULL
                    b <- plotly::as_widget(b)   ## from JSON back to R object
                    saveRDS(b, file=file)
                  }, message="saving plot object", value=0.2)
                } ## end of content
              ) ## end of object downloadHandler
            }
          } ## end of do object

          ##if(do.csv && is.null(download.csv) )  {
          if(do.csv)  {
              download.csv <- shiny::downloadHandler(
                filename = "data.csv",
                content = function(file) {
                  shiny::withProgress({
                    data <- csvFunc()
                    if(is.list(data)) data <- data[[1]]
                    write.csv(data, file=file)
                  }, message="Exporting to CSV", value=0.8)
                } ## end of content
              ) ## end of HTML downloadHandler
          } ## end of do HTML

          ##--------------------------------------------------------------------------------
          ##------------------------ OUTPUT ------------------------------------------------
          ##--------------------------------------------------------------------------------
          observeEvent(input$downloadOption, {
            if(input$downloadOption == "png"){
              output$download <- download.png
            }
            if(input$downloadOption == "pdf"){
              output$download <- download.pdf
            }
            if(input$downloadOption == "csv"){
              output$download <- download.csv
            }
            if(input$downloadOption == "html"){
              output$download <- download.html
            }
            if(input$downloadOption == "obj"){
              output$download <- download.obj
            }
          })

          ##--------------------------------------------------------------------------------
          ##---------------------------- UI ------------------------------------------------
          ##--------------------------------------------------------------------------------

          if(is.null(func2)) func2 <- func
          if(is.null(plotlib2)) plotlib2 <- plotlib
          if(length(height)==1) height <- c(height,700)
          if(length(width)==1)  width  <- c(width,1200)
          if(length(res)==1)    res    <- c(res, 1.3*res)

          res.1 <- res[1]
          res.2 <- res[2]

          ## width and height should actually be speficied in UI, not here.
          ifnotchar.int <- function(s) ifelse(grepl("[%]$|auto|vmin|vh|vw|vmax",s),s,as.integer(s))
          ##width = c(800,800)
          ##height = c(400,800)
          width.1  <- ifnotchar.int(width[1])
          width.2  <- ifnotchar.int(width[2])
          height.1 <- ifnotchar.int(height[1])
          height.2 <- ifnotchar.int(height[2])

          ## This sets the correct render and output functions for different
          ##plotting libraries.

          getRenderFunc <- function(plotlib) {
              switch(
                  plotlib,
                  generic = NULL,
                  htmlwidget = NULL,
                  plotly = plotly::renderPlotly,
##                  echarts4r = echarts4r::renderEcharts4r,
##                  scatterD3 = scatterD3::renderScatterD3,
                  pairsD3 = pairsD3::renderPairsD3,
                  visnetwork = visNetwork::renderVisNetwork,
                  ggplot = shiny::renderPlot,
                  grid = function(x) shiny::renderPlot(grid::grid.draw(x, recording=FALSE)),
                  iheatmapr = iheatmapr::renderIheatmap,
                  image = shiny::renderImage,
                  base = shiny::renderPlot,
                  shiny::renderPlot
              )
          }

          if(is.null(renderFunc)) {
              renderFunc <- getRenderFunc(plotlib)
          }

          if(is.null(renderFunc2)) {
              renderFunc2 <- getRenderFunc(plotlib2)
          }

          render = render2 = NULL

          if(1) {
              if(!is.null(func) && plotlib=="base") {
                  render <- shiny::renderPlot({ func() }, res=res.1)
                  ##render <- shiny::renderPlot({ func() }, res=res.1, width=width.1, height=height.1)
              }
              if(!is.null(func2) && plotlib2=="base") {
                  render2 <- shiny::renderPlot({ func2() }, res=res.2)
                  ##render2 <- shiny::renderPlot({ func2() }, res=res.2, width=width.2, height=height.2)
              }
              if(plotlib=="image") {
                  render <- shiny::renderImage( func(), deleteFile=FALSE)
              }
              if(!is.null(func2) && plotlib2=="image") {
                  render2 <- shiny::renderImage(func2(), deleteFile=FALSE)
              }

              ##if(grepl("renderCachedPlot",deparse(substitute(renderFunc)))) {
              if(grepl("cacheKeyExpr",head(renderFunc,1))) {
                  render <- shiny::renderCachedPlot(
                                       func(),
                                       cacheKeyExpr = {list(csvFunc())},
                                       res=res.1
                                   )
              }
              ##if(grepl("renderCachedPlot",deparse(substitute(renderFunc2)))) {
              if(grepl("cacheKeyExpr",head(renderFunc2,1))) {
                  render2 <- shiny::renderCachedPlot(
                                        func2(),
                                        cacheKeyExpr = {list(csvFunc())},
                                        res=res.2
                                    )
              }
          }


          if(is.null(render)) {
            if(plotlib == "plotly"){
              # If the plotting function is `plotly`, add the edit button
              render <- renderFunc({
                # By default remove plotly logo from all plots
                plot <- func() %>%
                  plotly::config(displaylogo = FALSE) %>%
                  plotly::plotly_build()
                #                plot <- plot %>%
                #                  plotly_default()
                # If there is already custom buttons, append the edit one
                # (issue #2210 plotly/plotly.R)
                if(inherits(plot$x$config$modeBarButtons, "list")){
                  for (y in 1:length(plot$x$config$modeBarButtons[[1]])) {
                    if(plot$x$config$modeBarButtons[[1]][[y]] == "toImage"){
                      plot$x$config$modeBarButtons[[1]][[y]] <- NULL
                      break
                    }
                  }
                  plot$x$config$modeBarButtons[[1]] <- append(plot$x$config$modeBarButtons[[1]],
                                                              list(octocat))
                } else { # Otherwise, apply the button regularly
                  plot <- plot %>%
                    plotly::config(modeBarButtonsToAdd = list(octocat),
                                   modeBarButtonsToRemove = c("zoomIn2d", "toImage"))
                }
                plot
              })
            } else {
              render <- renderFunc(func())
            }
          }

          if(is.null(render2) && !is.null(func2)) {
            if(plotlib2 == "plotly"){
              render2 <- renderFunc2({
                # By default remove plotly logo from all plots
                plot <- func2() %>% plotly::config(displaylogo = FALSE) %>%
                  plotly::plotly_build()
                # If there is already custom buttons, append the edit one
                # (issue #2210 plotly/plotly.R)
                if(inherits(plot$x$config$modeBarButtons, "list")){
                  for (y in 1:length(plot$x$config$modeBarButtons[[1]])) {
                    if(plot$x$config$modeBarButtons[[1]][[y]] == "toImage"){
                      plot$x$config$modeBarButtons[[1]][[y]] <- NULL
                      break
                    }
                  }
                  plot$x$config$modeBarButtons[[1]] <- append(plot$x$config$modeBarButtons[[1]],
                                                              list(octocat))
                } else { # Otherwise, apply the button regularly
                  plot <- plot %>%
                    plotly::config(modeBarButtonsToAdd = list(octocat),
                                   modeBarButtonsToRemove = c("zoomIn2d", "toImage"))
                }
                plot
              })
            } else {
              render2 <- renderFunc2(func2())
            }
          }

          output$renderfigure <- render
          output$renderpopup  <- render2

          ##--------------------------------------------------------------------------------
          ##---------------------------- RETURN VALUE --------------------------------------
          ##--------------------------------------------------------------------------------

          list(
              plotfun = func,
              plotfun2 = func2,
              .tmpfiles = c(pdf=PDFFILE, html=HTMLFILE),
              render = render,
              render2 = render2,
              download.pdf = download.pdf,
              download.png = download.png,
              download.html = download.html,
              download.csv = download.csv,
              ##getCaption = caption.fun,
              saveHTML = saveHTML,
              ## outputFunc = outputFunc,
              renderFunc = renderFunc
          )

      }
    )
}



##================================================================================
##====================== END OF FILE +============================================
##================================================================================
