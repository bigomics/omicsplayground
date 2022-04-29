##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
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
                       caption2=info.text, ## header=NULL,
                       plotlib = "base",
                       plotlib2 = "base",
                       outputFunc = NULL,
                       outputFunc2 = NULL,
                       no.download = FALSE,
                       download.fmt=c("png","pdf"), 
                       just.info=FALSE,
                       info.width="300px",
                       show.maximize = TRUE,
                       height = c(400,800),
                       width = c("auto","100%"),
                       pdf.width = 6,
                       pdf.height = 6
                       )
{
    ns <- shiny::NS(id)    

    if(is.null(plotlib2)) plotlib2 <- plotlib    
    if(length(height)==1) height <- c(height,800)
    if(length(width)==1)  width  <- c(width,1200)

    ##ifnotchar.int <- function(s) ifelse(grepl("[%]|auto|vh|vw|vmin|vmax",s),s,as.integer(s))
    ifnotchar.int <- function(s) suppressWarnings(
      ifelse(!is.na(as.integer(s)), paste0(as.integer(s),"px"), s))
    width.1  <- ifnotchar.int(width[1])
    width.2  <- ifnotchar.int(width[2])
    height.1 <- ifnotchar.int(height[1])
    height.2 <- ifnotchar.int(height[2])            
    
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
            iheatmap = iheatmapr::iheatmaprOutput,
            image = shiny::imageOutput,
            base = shiny::plotOutput,
            shiny::plotOutput        
        )
        FUN
    }
    
    if(is.null(outputFunc))  outputFunc  <- getOutputFunc(plotlib)
    if(is.null(outputFunc2)) outputFunc2 <- getOutputFunc(plotlib2)    
    
    ##--------------------------------------------------------------------------------
    ##------------------------ BUTTONS -----------------------------------------------
    ##--------------------------------------------------------------------------------
    
    ##if(is.null(inputs) || length(inputs)==0 ) inputs <- ""
    options.button <- ""
    
    if(!just.info && !is.null(options) && length(options)>0) {
        options.button <- shinyWidgets::dropdownButton(
            ##shiny::tags$h3("Options"),
            options,
            ##shiny::br(),
            ##dload,
            circle = TRUE, size = "xs", ## status = "danger",
            ##icon = shiny::icon("gear"),
            ##icon = shiny::icon("align-justify"),
            icon = shiny::icon("bars"),
            width = "250px",
            inputId = ns("options"),
            tooltip = shinyWidgets::tooltipOptions(title = "Settings", placement = "right")
        )
    }
    dload.csv = dload.pdf = dload.png = dload.html = NULL
    if("pdf" %in% download.fmt)  dload.pdf  <- shiny::downloadButton(ns("pdf"), "PDF")
    if("png" %in% download.fmt)  dload.png  <- shiny::downloadButton(ns("png"), "PNG")
    if("html" %in% download.fmt) dload.html <- shiny::downloadButton(ns("html"), "HTML")
    if("csv" %in% download.fmt)  dload.csv  <- shiny::downloadButton(ns("csv"), "CSV")

    pdf_size = NULL
    if(TRUE || plotlib!="base") {
        pdf_size <- shiny::tagList(
            shiny::fillRow(
                shiny::numericInput(ns("pdf_width"), "PDF width", pdf.width, 1, 20, 1, width='100%'),
                shiny::numericInput(ns("pdf_height"), "height", pdf.height, 1, 20, 1, width='100%')
            ),
            shiny::br(),shiny::br(),shiny::br()
        )
    }
    
    dload.button <- shinyWidgets::dropdownButton(
        dload.pdf,
        dload.png,
        dload.csv,
        dload.html,
        pdf_size,
        circle = TRUE, size = "xs", ## status = "danger",
        icon = shiny::icon("download"), width = "40px", right=FALSE,
        tooltip = shinyWidgets::tooltipOptions(title = "Download", placement = "right")
    )

    if(no.download || length(download.fmt)==0 ) dload.button <- ""

    if(!is.null(label) && label!="") label <- paste0("(",label,")")
    label1 = shiny::HTML(paste0("<span class='module-label'>",label,"</span>"))
        
    ##zoom.button <- shinyWidgets::prettyCheckbox(inputId=ns("zoom"),label=NULL,value=FALSE)
    zoom.button <- NULL
    if(1 && show.maximize) {        
        zoom.button <- modalTrigger(ns("zoombutton"),
            ns("plotPopup"),
            icon("window-maximize"),
            class="btn-circle-xs"
        )
##zoom.button <- shinyBS::tipify(zoom.button, "Maximize", placement="right")  ## weird Cairo error!!!
        ##zoom.button <- with_tippy(zoom.button, "Maximize") ## not consistent...       
    }
    
    buttons <- shiny::fillRow(
        flex = c(NA,1,NA,NA,NA,NA),
        ##flex=c(NA,NA,1),
        label1,
##      shiny::HTML(paste("<center>",title,"</center>")),
        shiny::div(class='plotmodule-title', title=title, title),        
        shinyWidgets::dropdownButton(
            shiny::tags$p(shiny::HTML(info.text)),
            shiny::br(),
            circle = TRUE, size = "xs", ## status = "danger",
            icon = shiny::icon("info"), width = info.width,
            inputId = ns("info"), right=FALSE,
            tooltip = shinyWidgets::tooltipOptions(title = "Info", placement = "right")
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
        if(any(class(caption2)=="reactive")) {
            caption2 <- caption2()
        }
        if(any(class(caption2)=="character")) {
            caption2 <- shiny::HTML(caption2)
            caption2 <- shiny::div(caption2, class="caption2")                         
        }
        shiny::tagList(
                   outputFunc2(ns("renderpopup"), width=w, height=h, inline=FALSE),
                   caption2
        )
    }

    modaldialog.style <- paste0("#",ns("plotPopup")," .modal-dialog {width:",width.2,";}")
    modalbody.style <- paste0("#",ns("plotPopup")," .modal-body {min-height:",height.2,"; padding:60px 300px;}")
    modalcontent.style <- paste0("#",ns("plotPopup")," .modal-content {width:100vw;}")    
    modalfooter.none <- paste0("#",ns("plotPopup")," .modal-footer{display:none;}")
    
    if(any(class(caption)=="reactive")) {
        caption <- caption()
    }
    if(any(class(caption)=="character")) {
        caption <- shiny::HTML(caption)
        caption <- shiny::div(caption, class="caption")                 
    }

    div( class="plotmodule",
        shiny::fillCol(
               flex = c(NA,1,NA,0.001,NA),
               height = height.1,
               div( buttons, class="plotmodule-header"),
               outputFunc(ns("renderfigure"), width=width.1, height=height.1) %>% shinycssloaders::withSpinner(),                              
               caption,
               shiny::div(class="popup-plot",
                          modalUI(
                                ns("plotPopup"),
                                title, 
                                size="fullscreen",
                                popupfigUI()
                            )
                          ),
               shiny::tagList(
                          shiny::tags$head(shiny::tags$style(modaldialog.style)),
                          shiny::tags$head(shiny::tags$style(modalbody.style)),
                          shiny::tags$head(shiny::tags$style(modalcontent.style)),
                          shiny::tags$head(shiny::tags$style(modalfooter.none))
                      )
              )
    )

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
         ##height = c(640,800),
         ##width = c("auto",1400),
         res=c(100,170),
         download.pdf = NULL,
         download.png = NULL,
         download.html = NULL,
         download.csv = NULL,
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
          ##------------------------ FIGURE ------------------------------------------------
          ##--------------------------------------------------------------------------------
          
          ## these engines cannot (yet) provide html
          if(plotlib %in% c("base")) {    
              download.fmt <- setdiff(download.fmt, c("html"))
          }
          
          do.pdf = "pdf" %in% download.fmt
          do.png = "png" %in% download.fmt
          do.html = "html" %in% download.fmt
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
                                                 dbg("[downloadHandler.PNG] export to PNG done!")                    
                                             }, message="exporting to PNG", value=0.8)
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
                                                 message("[plotModule] export to PDF done!")
                                             }, message="exporting to PDF", value=0.8)
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
                  DT::saveWidget( plotly::ggplotly(p), file = HTMLFILE);
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
                                                      DT::saveWidget( plotly::ggplotly(p), file = HTMLFILE);
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
                                              }, message="exporting to HTML", value=0.8)
                                          } ## end of content
                                      ) ## end of HTML downloadHandler
          } ## end of do HTML

          ##if(do.csv && is.null(download.csv) )  {
          if(do.csv)  {
              download.csv <- shiny::downloadHandler(
                                         filename = "data.csv",
                                         content = function(file) {
                                             shiny::withProgress({
                                                 data <- csvFunc()
                                                 if(is.list(data)) data <- data[[1]]
                                                 ##file.copy(CSVFILE, file, overwrite=TRUE)
                                                 write.csv(data, file=file)
                                             }, message="exporting to CSV", value=0.8)
                                         } ## end of content
                                     ) ## end of HTML downloadHandler
          } ## end of do HTML
          
          ##--------------------------------------------------------------------------------
          ##------------------------ OUTPUT ------------------------------------------------
          ##--------------------------------------------------------------------------------

          if(!is.null(download.csv))  output$csv  <- download.csv
          if(!is.null(download.pdf))  output$pdf  <- download.pdf
          if(!is.null(download.png))  output$png  <- download.png
          if(!is.null(download.html)) output$html <- download.html

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
                  visNetwork = visNetwork::renderVisNetwork,
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
              ##render <- eval(parse(text=renderFunc))(func())
              render <- renderFunc(func())              
          }
          if(is.null(render2) && !is.null(func2)) {
              ##render2 <- eval(parse(text=renderFunc2))(func2())
              render2 <- renderFunc2(func2())              
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
              ##buttons = buttons,
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
