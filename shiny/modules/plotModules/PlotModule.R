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
    ##renderFUN=NULL
    outputFUN="shiny::plotOutput"
    ##csvFUN=NULL
    ##renderFUN2=NULL
    outputFUN2="shiny::plotOutput" 
    no.download = FALSE
    download.fmt=c("png","pdf") 
    just.info=FALSE
    info.width="300px"
    show.maximize = TRUE
    height = c(640800)
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
                       plotlib2 = NULL,
                       ##renderFUN=NULL,
                       outputFUN="shiny::plotOutput",
                       ##csvFUN=NULL,
                       ##renderFUN2=NULL,
                       outputFUN2="shiny::plotOutput", 
                       no.download = FALSE,
                       download.fmt=c("png","pdf"), 
                       just.info=FALSE,
                       info.width="300px",
                       show.maximize = TRUE,
                       height = c(640,800),
                       width = c("auto",1400),
                       pdf.width = 6,
                       pdf.height = 6                        
                       ##res=c(72,100)
                       )
{
    ns <- shiny::NS(id)    

    if(is.null(plotlib2)) plotlib2 <- plotlib    
    if(length(height)==1) height <- c(height,700)
    if(length(width)==1)  width  <- c(width,1200)

    ifnotchar.int <- function(s) ifelse(grepl("[%]|auto",s),s,as.integer(s))
    width.1  <- ifnotchar.int(width[1])
    width.2  <- ifnotchar.int(width[2])
    height.1 <- ifnotchar.int(height[1])
    height.2 <- ifnotchar.int(height[2])
    

    getOutputFUN <- function(plotlib)
    {
        switch(
            plotlib,
            generic = NULL,
            htmlwidget = NULL,        
            plotly = "plotly::plotlyOutput",
            echarts4r = "echarts4r::echarts4rOutput",
            scatterD3 = "scatterD3::scatterD3Output",
            pairsD3 = "pairsD3::pairsD3Output",
            visnetwork = "visNetwork::visNetworkOutput",
            ggplot = "shiny::plotOutput",
            grid = "shiny::plotOutput",
            iheatmap = "iheatmapr::iheatmaprOutput",
            image = "shiny::imageOutput",
            base = "shiny::plotOutput",
            "shiny::plotOutput"        
        )
    }
    
    if(is.null(outputFUN))  outputFUN <- getOutputFUN(plotlib)
    if(is.null(outputFUN2)) outputFUN2 <- getOutputFUN(plotlib2)    
    
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
    if(plotlib!="base") {
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
    label1 = shiny::HTML(paste0("<span class='module-label'>",label,"</span>"))
        
    ##zoom.button <- shinyWidgets::prettyCheckbox(inputId=ns("zoom"),label=NULL,value=FALSE)
    zoom.button <- NULL
    if(show.maximize) {        
        zoom.button <- shiny::actionButton(inputId=ns("zoombutton"),label=NULL,
                                    icon=icon("window-maximize"),
                                    class="btn-circle-xs")
        zoom.button <- shinyBS::tipify(zoom.button, "Maximize", placement="right")
    }
    
    ##output$renderbuttons <- shiny::renderUI({
    ## button layout
    buttons <- shiny::fillRow(
        flex = c(NA,NA,NA,NA,NA,1),
        ##flex=c(NA,NA,1),
        label1,
        ##div( class="button-group", style="display: inline-block; float: left;",
        shinyWidgets::dropdownButton(
            shiny::tags$p(shiny::HTML(info.text)),
            shiny::br(),
            circle = TRUE, size = "xs", ## status = "danger",
            icon = shiny::icon("info"), width = info.width,
            inputId = ns("info"), right=FALSE,
            tooltip = shinyWidgets::tooltipOptions(title = "Info", placement = "right")
        ),
        options.button,
        shiny::div(class='download-button', dload.button),
        shiny::div(class='zoom-button', zoom.button),
        ##),
        shiny::HTML(paste("<center>",title,"</center>"))
        ##HTML(paste("<center><strong>",title,"</strong></center>"))
        ##shiny::br()
        ##inputs
        ##selectInput("sel123","number",1:10)
    )
    ##return(ui)
    ##})
    
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
        r <- eval(parse(text=outputFUN2))(ns("renderpopup"), width=w, height=h)
        if(any(class(caption2)=="reactive")) {
            caption2 <- caption2()
        }
        if(any(class(caption2)=="character")) {
            caption2 <- shiny::HTML(caption2)
        }
        shiny::tagList(
                   r,
                   shiny::div( caption2, class="caption2") 
               )
    }

    ##output$widget <- shiny::renderUI({    
    ##cat("[output$widget::renderUI] ns.zoombutton = ",ns("zoombutton"),"\n")
    mtop <- paste0("margin: 400px 20px 20px 20px;")
    modaldialog.style <- paste0("#",ns("plotPopup")," .modal-dialog {width:",width.2+40,"px;}")
    modalbody.style <- paste0("#",ns("plotPopup")," .modal-body {min-height:",height.2+40,"px;}")
    ##modaldialog.style <- paste0("#",ns("plotPopup")," .modal-dialog {width:",width.2,"px;}")
    ##modalbody.style <- paste0("#",ns("plotPopup")," .modal-body {min-height:",height.2,"px;}")
    modalfooter.none <- paste0("#",ns("plotPopup")," .modal-footer{display:none;}")
    
    if(any(class(caption)=="reactive")) {
        caption <- caption()
    }
    if(any(class(caption)=="character")) {
        caption <- shiny::HTML(caption)
    }
    
    shiny::fillCol(
               flex = c(NA,NA,1,NA,NA,0.001),
               height = height.1,
               shiny::tagList(
                          shiny::tags$head(shiny::tags$style(modaldialog.style)),
                          shiny::tags$head(shiny::tags$style(modalbody.style)),            
                          shiny::tags$head(shiny::tags$style(modalfooter.none))
                      ),
               
               buttons,
               ##render,
               eval(parse(text=outputFUN))(ns("renderfigure"), width=width.1, height=height.1),
               shiny::br(),
               shiny::div(caption, class="caption"),          
               shiny::div(class="popup-plot",
                          shinyBS::bsModal(ns("plotPopup"), title, size="l",
                                           ns("zoombutton"),
                                           ##shiny::tagList(shiny::uiOutput(ns("popupfig")))
                                           popupfigUI()
                                           )
                          )
           )

}

PlotModuleServer <- function(
         id,
         func,
         func2=NULL,
         ##info.text="Figure",
         ##title="",
         ##inputs=NULL, 
         ##options = NULL,
         ##label="",
         ##caption="",
         ##caption2=info.text, ## header=NULL,
         plotlib = "base",
         plotlib2 = NULL,
         renderFUN=NULL,
         ##outputFUN=NULL,
         csvFUN=NULL,
         renderFUN2=NULL,
         ##outputFUN2=NULL, 
         ##no.download = FALSE,
         download.fmt=c("png","pdf"), 
         ##just.info=FALSE,
         ##info.width="300px",
         ##show.maximize = TRUE,
         ##height = c(640,800),
         ##width = c("auto",1400),
         res=c(72,100),
         download.pdf = NULL,
         download.png = NULL,
         download.html = NULL,
         download.csv = NULL,
         pdf.width=8,
         pdf.height=6,
         pdf.pointsize=12,
         add.watermark=FALSE ) {
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
          ##do.csv  = "csv" %in% download.fmt && !is.null(csvFUN)
          do.csv = !is.null(csvFUN)
          
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
                                                     p <- func()
                                                     png(PNGFILE, width=pdf.width*100, height=pdf.height*100,
                                                         pointsize=1.2*pdf.pointsize)
                                                     print(p) 
                                                     dev.off() 
                                                 } else if(plotlib=="grid") {
                                                     p <- func()
                                                     png(PNGFILE, width=pdf.width*100, height=pdf.height*100,
                                                         pointsize=1.2*pdf.pointsize)
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
                                                     ##cat("downloadHandler:: exporting to base plot to PNG\n")
                                                     ## NEEEDS FIX!! pdf generating should be done
                                                     ## just here, not anymore in the
                                                     ## renderPlot. But cannot get it to work (IK 19.10.02)
                                                     if(0) {
                                                         dbg("[downloadHandler.PNG] creating new PNG device")
                                                         png(file=PNGFILE, width=80*pdf.width, height=80*pdf.height)
                                                         func()
                                                         ##dev.copy2pdf(file=PDFFILE, width=pdf.width, height=pdf.height)
                                                         ##plot(sin)
                                                         dev.off()  ## important!!
                                                     }
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
                                                     ## NEEEDS FIX!! pdf generating should be done
                                                     ## just here, not anymore in the
                                                     ## renderPlot. But cannot get it to work (IK 19.10.02)
                                                     if(0) {
                                                         pdf(file=PDFFILE, width=pdf.width, height=pdf.height,
                                                             pointsize=pdf.pointsize)
                                                         func()
                                                         ##plot(sin)
                                                         dev.off()  ## important!!
                                                     }
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
                                                 data <- csvFUN()                    
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
          ##    ifnotchar.int <- function(s) ifelse(grepl("[%]|auto",s),s,as.integer(s))
          ##    width.1  <- ifnotchar.int(width[1])
          ##    width.2  <- ifnotchar.int(width[2])
          ##    height.1 <- ifnotchar.int(height[1])
          ##    height.2 <- ifnotchar.int(height[2])
          
          ## This sets the correct render and output functions for different
          ##plotting libraries.

          getRenderFUN <- function(plotlib) {
              switch(
                  plotlib,
                  generic = NULL,
                  htmlwidget = NULL,        
                  plotly = "plotly::renderPlotly",
                  echarts4r = "echarts4r::renderEcharts4r",
                  scatterD3 = "scatterD3::renderScatterD3",
                  pairsD3 = "pairsD3::renderPairsD3",
                  visNetwork = "visNetwork::renderVisNetwork",
                  ggplot = "shiny::renderPlot",
                  grid = "function(x) shiny::renderPlot(grid::grid.draw(x, recording=FALSE))",
                  iheatmapr = "iheatmapr::renderIheatmap",
                  image = "shiny::renderImage",
                  base = "shiny::renderPlot",
                  "shiny::renderPlot"
              )
          }
          
          if(is.null(renderFUN)) {
              renderFUN <- getRenderFUN(plotlib)
          }

          if(is.null(renderFUN2)) {
              if(plotlib2 != plotlib ) {    
                  renderFUN2 <- getRenderFUN(plotlib2)
              } else {
                  renderFUN2 <- renderFUN
              }
          }
          
          if(0) {
              dbg("[plotModule] title=",title)
              dbg("[plotModule] plotlib=",plotlib)
              dbg("[plotModule] plotlib2=",plotlib2)
              dbg("[plotModule] outputFUN=",outputFUN)
              dbg("[plotModule] outputFUN2=",outputFUN2)
              dbg("[plotModule] renderFUN=",renderFUN)
              dbg("[plotModule] renderFUN2=",renderFUN2)
              dbg("[plotModule] is.null.renderFUN=",is.null(renderFUN))
              dbg("[plotModule] is.null.renderFUN2=",is.null(renderFUN2))
          }
          
          ##outputFUN <- sub(".*::","",outputFUN)
          render = render2 = NULL   
          if(1 && plotlib=="base") {
              ## For base plots we need to export PDF/PNG here. For some
              ## reason does not work in the downloadHandler...
              ##
              render <- shiny::renderPlot({
                  ##pdf.width0  <- shiny::isolate(input$pdf_width)
                  ##pdf.height0 <- shiny::isolate(input$pdf_height)                
                  pdf.width0  <- pdf.width
                  pdf.height0 <- pdf.height

                  ## save base render in object
                  ##pdf(NULL)
                  ##dev.control(displaylist="enable")
                  func()
                  ##p1.base <- recordPlot()
                  ##invisible(dev.off())
                  if(1) {
                      suppressWarnings( suppressMessages(
                          if(do.pdf) {
                              ##pdf(file=PDFFILE, width=pdf.width, height=pdf.height, pointsize=pdf.pointsize)
                              dev.print(pdf, file=PDFFILE, width=pdf.width0, height=pdf.height0,
                                        pointsize=pdf.pointsize)
                              ##p1.base
                              ##dev.off()  ## important!!
                          }
                      ))
                  }
                  if(1) {
                      suppressWarnings( suppressMessages(
                          if(do.png) {
                              ##png(file=PNGFILE, width=pdf.width*100, height=pdf.height*100,pointsize=pdf.pointsize)
                              dev.print(png, file=PNGFILE,
                                        width = pdf.width0*72, height = pdf.height0*72,
                                        pointsize=1.4*pdf.pointsize)
                              ##p1.base
                              ##dev.off()  ## important!!
                          }
                      ))
                  }
                  ##p1.base
              }, res=res.1)
          }
          
          if(!is.null(func2) && plotlib2=="base") {
              render2 <- shiny::renderPlot({
                  func2()
              }, res=res.2)
          }
          if(plotlib=="image") {
              render <- shiny::renderImage( func(), deleteFile=FALSE)
          }
          if(!is.null(func2) && plotlib2=="image") {
              render2 <- shiny::renderImage(func2(), deleteFile=FALSE)
          }
          if(grepl("renderCachedPlot",renderFUN)) {
              render <- shiny::renderCachedPlot(
                                   func(), cacheKeyExpr = {list(csvFUN())}, res=res.1
                               )
          }
          if(grepl("renderCachedPlot",renderFUN2)) {
              render2 <- shiny::renderCachedPlot(
                                    func2(), cacheKeyExpr = {list(csvFUN())}, res=res.2
                                )
          }

          if(is.null(render)) {
              render <- eval(parse(text=renderFUN))(func())
          }
          if(is.null(render2) && !is.null(func2)) {
              render2 <- eval(parse(text=renderFUN2))(func2())
          }
          
          output$renderfigure <- render
          output$renderpopup  <- render2

          
          
          ##--------------------------------------------------------------------------------
          ##---------------------------- RETURN VALUE --------------------------------------
          ##--------------------------------------------------------------------------------

          res <- list(
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
              outputFUN = outputFUN,
              renderFUN = renderFUN
          )
          return(res)
      }
    )
}



##================================================================================
##====================== END OF FILE +============================================
##================================================================================
