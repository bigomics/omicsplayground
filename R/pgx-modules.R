########################################################################
## PLOT/TABLE MODULES
########################################################################

WATERMARK=TRUE
##WATERMARK=FALSE

addWatermarkPlotly <- function(p) {
    add_annotations(
        p,
        ##yshift = -100, 
        ##x = 1, y=-0.05, xanchor = "right", 
        x = 0.5, y = 0.5, xanchor = "middle", 
        xref="paper", yref="paper", showarrow=FALSE,
        text = "Omics Playground\nFree version",
        font = list(size=64, family='Arial', color="#AAAAAA22")
        ##font = list(size=12, color="#00000066")
    )
}
    
##================================================================================
##================================================================================
##================================================================================

plotWidget <- function(id) {
    ns <- NS(id)
    uiOutput(ns("widget"))
}
    
plotModule <- function(input, output, session, ## ns=NULL,
                       func, func2=NULL, 
                       info.text="Info text", title="", 
                       inputs=NULL, options = NULL, label="",
                       caption="", caption2="", ## header=NULL,
                       plotlib = "base", renderFunc=NULL, outputFunc=NULL, csvFunc=NULL,
                       no.download = FALSE, download.fmt=c("png","pdf"), 
                       just.info=FALSE, info.width="300px", show.maximize = TRUE,
                       height = c(400,720), width = c("auto",1080), res=c(72,100),
                       download.pdf = NULL, download.png = NULL,
                       download.html = NULL, download.csv = NULL,
                       pdf.width=8, pdf.height=8, pdf.pointsize=12)
{
    ns <- session$ns    
    
    ##--------------------------------------------------------------------------------
    ##------------------------ BUTTONS -----------------------------------------------
    ##--------------------------------------------------------------------------------
    
    if(is.null(inputs) || length(inputs)==0 ) inputs <- ""
    options.button <- ""
    
    if(!just.info && !is.null(options) && length(options)>0) {
        options.button <- dropdownButton(
            ##tags$h3("Options"),
            options,
            ##br(),
            ##dload,
            circle = TRUE, size = "xs", ## status = "danger",
            ##icon = icon("gear"),
            ##icon = icon("align-justify"),
            icon = icon("bars"),
            width = "250px",
            inputId = ns("options"),
            tooltip = tooltipOptions(title = "Settings", placement = "right")
        )
    }
    
    dload.csv = dload.pdf = dload.png = dload.html = NULL
    if("pdf" %in% download.fmt)  dload.pdf  <- downloadButton(ns("pdf"), "PDF")
    if("png" %in% download.fmt)  dload.png  <- downloadButton(ns("png"), "PNG")
    if("html" %in% download.fmt) dload.html <- downloadButton(ns("html"), "HTML")
    ##if("csv" %in% download.fmt)  dload.csv  <- downloadButton(ns("csv"), "CSV")
    if(!is.null(csvFunc))  dload.csv  <- downloadButton(ns("csv"), "CSV")
    
    dload.button <- dropdownButton(
        dload.pdf,
        dload.png,
        dload.csv,
        dload.html,
        circle = TRUE, size = "xs", ## status = "danger",
        icon = icon("download"), width = "40px", right=FALSE,
        tooltip = tooltipOptions(title = "Download", placement = "right")
    )

    if(no.download || length(download.fmt)==0 ) dload.button <- ""
    label1 = HTML(paste0("<span class='module-label'>",label,"</span>"))
    
    require(shinyWidgets)
    ##zoom.button <- prettyCheckbox(inputId=ns("zoom"),label=NULL,value=FALSE)
    zoom.button <- NULL
    if(show.maximize) {
        require(shinyBS)
        zoom.button <- actionButton(inputId=ns("zoombutton"),label=NULL,
                                    icon=icon("window-maximize"),
                                    class="btn-circle-xs")
        zoom.button <- tipify(zoom.button, "Maximize", placement="right")
    }
    
    ##output$renderbuttons <- renderUI({
        ## button layout
    buttons <- fillRow(
        flex=c(NA,NA,NA,NA,NA,1),
        ##flex=c(NA,NA,1),
        label1,
        ##div( class="button-group", style="display: inline-block; float: left;",
        dropdownButton(
            tags$p(HTML(info.text)),
            br(),
            circle = TRUE, size = "xs", ## status = "danger",
            icon = icon("info"), width = info.width,
            inputId = ns("info"), right=FALSE,
            tooltip = tooltipOptions(title = "Info", placement = "right")
        ),
        options.button,
        div(class='download-button', dload.button),
        div(class='zoom-button', zoom.button),
        ##),
        HTML(paste("<center>",title,"</center>"))
        ##HTML(paste("<center><strong>",title,"</strong></center>"))
        ##br()
        ##inputs
        ##selectInput("sel123","number",1:10)
    )
    ##return(ui)
    ##})
    
    ##--------------------------------------------------------------------------------
    ##------------------------ FIGURE ------------------------------------------------
    ##--------------------------------------------------------------------------------
    
    require(shinyWidgets)
    require(shinyBS)
    require(webshot)
    
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
    
    ##outputFunc = renderFunc = NULL
    if(plotlib == "generic") {
        if(is.null(renderFunc)) stop("'generic' class must provide renderFunc")
        if(is.null(outputFunc)) stop("'generic' class must provide outputFunc")
    } else if(plotlib == "htmlwidget") {
        if(is.null(renderFunc)) stop("'htmlwidget' class must provide renderFunc")
        if(is.null(outputFunc)) stop("'htmlwidget' class must provide outputFunc")
    } else if(plotlib == "plotly") {
        ##render <- renderPlotly({ func() })
        require(plotly)
        outputFunc = "plotlyOutput"
        renderFunc = "renderPlotly"
    } else if(plotlib=="scatterD3") {
        require(scatterD3)
        renderFunc="renderScatterD3"
        outputFunc="scatterD3Output"
    } else if(plotlib=="pairsD3") {
        require(pairsD3)
        renderFunc="renderPairsD3"
        outputFunc="pairsD3Output"
    } else if(plotlib == "visnetwork") {
        require(visNetwork)
        ##render <- renderVisNetwork({ func() })
        outputFunc="visNetworkOutput"
        renderFunc = "renderVisNetwork"
    } else if(plotlib %in% c("ggplot","ggplot2")) {
        ##render <- renderPlot({ plot(func())}, res=res)
        outputFunc="plotOutput"
        renderFunc="function(x) renderPlot(plot(x))"
    } else if(plotlib == "iheatmapr") {
        require(iheatmapr)
        ##render <- renderIheatmap({ func() })
        outputFunc="iheatmaprOutput"
        renderFunc="renderIheatmap"
    } else if(plotlib == "image") {
        require(iheatmapr)
        ##render <- renderIheatmap({ func() })
        outputFunc="imageOutput"
        renderFunc="renderImage"
    } else {
        ## Base plotting
        ## render <- renderPlot({
        ##     par(mar=c(0,0,0,0),oma=c(0,0,0,0))
        ##     func()
        ## }, res=res, width=width, height=height)
        outputFunc="plotOutput"
        ##renderFunc="renderPlot"
        ##renderFunc="function(x) renderPlot(x, res=res, width=width, height=height)"
        renderFunc="renderPlot"
    }

    ##============================================================
    ##=============== Download Handlers ==========================
    ##============================================================
    ## download.pdf = NULL
    ##download.png = download.html = NULL

    if(do.png && is.null(download.png)) {
        download.png <- downloadHandler(
            filename = "plot.png",
            content = function(file) {
                withProgress({
                    ## unlink(PNGFILE) ## do not remove!
                    if(plotlib=="plotly") {
                        cat("downloadHandler:: exporting plotly to PNG\n")
                        p <- func()
                        p$width = pdf.width * 100
                        p$height = pdf.height * 100
                        ##is.plotly3d <- class(p)[1]=="plotly" && all( c("x","y","z") %in% names(p$x$attrs[[1]]))
                        err <- try(export(p, PNGFILE))  ## deprecated but still works...
                        if(class(err)=="try-error") {
                            cat("downloadHandler:: export failed, trying webshot...\n")
                            htmlwidgets::saveWidget(p, HTMLFILE) 
                            err2 <- try(webshot(HTMLFILE,vwidth=pdf.width*100,
                                                vheight=pdf.height*100,PNGFILE))
                            if(class(err2)=="try-error") {
                                png(PNGFILE)
                                text(0.5,0.5,"PNG export error (plotly)")
                                dev.off()
                            }
                        }
                    } else if(plotlib=="iheatmapr") {
                        cat("downloadHandler:: exporting iheatmapR to PNG\n")
                        p <- func()
                        save_iheatmap(p, vwidth=pdf.width*80,vheight=pdf.height*80,PNGFILE)
                    } else if(plotlib=="visnetwork") {
                        cat("downloadHandler:: exporting visnetwork to PNG\n")
                        p <- func()
                        visSave(p, HTMLFILE)
                        webshot(HTMLFILE,vwidth=pdf.width*100,vheight=pdf.height*100,PNGFILE)
                    } else if(plotlib %in% c("htmlwidget","pairsD3","scatterD3")) {
                        cat("downloadHandler:: exporting htmlwidget to PNG\n")
                        p <- func()
                        htmlwidgets::saveWidget(p, HTMLFILE)
                        webshot(HTMLFILE, vwidth=pdf.width*100,vheight=pdf.height*100,PNGFILE)
                    } else if(plotlib %in% c("ggplot","ggplot2")) {
                        cat("downloadHandler:: exporting ggplot to PNG\n")
                        p <- func()
                        ##p = addSignature(p)                             
                        ##ggsave(PDFFILE, width=pdf.width, height=pdf.height)
                        png(PNGFILE, width=pdf.width, height=pdf.height, pointsize=pdf.pointsize)
                        print(p) 
                        dev.off() 
                    } else if(plotlib=="generic") {
                        cat("downloadHandler:: generic to PNGFILE = ",PNGFILE,"\n")
                        ## generic function should produce PDF inside plot func()
                        ##
                    } else if(plotlib=="base") {
                        ##cat("downloadHandler:: exporting to base plot to PDF\n")
                        cat("downloadHandler:: exporting base plot to PNG\n")
                        
                        ## NEEEDS FIX!! pdf generating should be done
                        ## just here, not anymore in the
                        ## renderPlot. But cannot get it to work (IK 19.10.02)
                        if(0) {
                            cat("downloadHandler:: creating new PDF device\n")
                            add.ADVERTISEMENT = TRUE
                            add.ADVERTISEMENT = FALSE
                            pdf(file=PDFFILE, width=pdf.width, height=pdf.height, pointsize=pdf.pointsize)
                            func()
                            ##dev.copy2pdf(file=PDFFILE, width=pdf.width, height=pdf.height)
                            plot(sin)
                            ## ##mtext( caption.fun(), outer=TRUE, side=1,line=-3.5, cex=1, xpd=NA)
                            ## if(add.ADVERTISEMENT) {
                            ##     par(mfrow=c(1,1),par=c(1,1,1,1)*0.5)
                            ##     frame()
                            ##     motto <- pgx.randomSlogan(b=30)
                            ##     mex = min(pdf.width,pdf.height)/8.0
                            ##     text(0.5, 0.95, "Created with BigOmics Playground", cex=1.2*mex)
                            ##     text(0.5, 0.70, motto, cex=2*mex, font=4)
                            ##     text(0.5, 0.25, "BigOmics Analytics", cex=1.5*mex, font=2)
                            ##     text(0.5, 0.15, "Self-service bioinformatics solutions", cex=1.2*mex)
                            ##     text(0.5, 0.05, "www.bigomics.ch", cex=1.2*mex, font=3, xpd=NA)
                            ## }
                            dev.off()  ## important!!
                        }
                    } else { ## end base                
                        png(PNGFILE, pointsize=pdf.pointsize)
                        plot.new()
                        mtext("Error. PNG not available.",line=-8)
                        dev.off()
                    }
                    
                    ## finally copy to final exported file
                    file.copy(PNGFILE,file)
                }, message="exporting to PNG", value=0.8)
            } ## content 
        ) ## PNG downloadHandler
    } ## end if do.png

    if(do.pdf && is.null(download.pdf) ) {        
        download.pdf <- downloadHandler(
            filename = "plot.pdf",
            content = function(file) {
                withProgress({
                    ## unlink(PDFFILE) ## do not remove!
                    if(plotlib=="plotly") {
                        cat("downloadHandler:: exporting plotly to PDF\n")
                        p <- func()
                        p$width = pdf.width * 100
                        p$height = pdf.height * 100
                        ##is.plotly3d <- class(p)[1]=="plotly" && all( c("x","y","z") %in% names(p$x$attrs[[1]]))
                        err <- try(export(p, PDFFILE))  ## deprecated but still works...
                        if(class(err)=="try-error") {
                            cat("downloadHandler:: export failed, trying webshot...\n")
                            htmlwidgets::saveWidget(p, HTMLFILE) 
                            err2 <- try(webshot(HTMLFILE,vwidth=pdf.width*100,
                                                vheight=pdf.height*100,PDFFILE))
                            if(class(err2)=="try-error") {
                                pdf(PDFFILE)
                                text(0.5,0.5,"PDF export error (plotly)")
                                dev.off()
                            }
                        }
                    } else if(plotlib=="iheatmapr") {
                        cat("downloadHandler:: exporting iheatmapR to PDF\n")
                        p <- func()
                        save_iheatmap(p, vwidth=pdf.width*80,vheight=pdf.height*80,PDFFILE)
                    } else if(plotlib=="visnetwork") {
                        cat("downloadHandler:: exporting visnetwork to PDF\n")
                        p <- func()
                        visSave(p, HTMLFILE)
                        webshot(HTMLFILE,vwidth=pdf.width*100,vheight=pdf.height*100,PDFFILE)
                    } else if(plotlib %in% c("htmlwidget","pairsD3","scatterD3")) {
                        cat("downloadHandler:: exporting htmlwidget to PDF\n")
                        p <- func()
                        htmlwidgets::saveWidget(p, HTMLFILE)
                        webshot(HTMLFILE, vwidth=pdf.width*100,vheight=pdf.height*100,PDFFILE)
                    } else if(plotlib %in% c("ggplot","ggplot2")) {
                        cat("downloadHandler:: exporting ggplot to PDF\n")
                        p <- func()
                        ##p = addSignature(p)                             
                        ##ggsave(PDFFILE, width=pdf.width, height=pdf.height)
                        pdf(PDFFILE, width=pdf.width, height=pdf.height, pointsize=pdf.pointsize)
                        print(p) 
                        dev.off() 
                    } else if(plotlib=="generic") {
                        cat("downloadHandler:: generic to PDFFILE = ",PDFFILE,"\n")
                        ## generic function should produce PDF inside plot func()
                        ##
                    } else if(plotlib=="base") {
                        ##cat("downloadHandler:: exporting to base plot to PDF\n")
                        cat("downloadHandler:: exporting base plot to PDF\n")
                        ## NEEEDS FIX!! pdf generating should be done
                        ## just here, not anymore in the
                        ## renderPlot. But cannot get it to work (IK 19.10.02)
                        if(0) {
                            cat("downloadHandler:: creating new PDF device\n")
                            add.ADVERTISEMENT = TRUE
                            add.ADVERTISEMENT = FALSE
                            pdf(file=PDFFILE, width=pdf.width, height=pdf.height, pointsize=pdf.pointsize)
                            func()
                            ##dev.copy2pdf(file=PDFFILE, width=pdf.width, height=pdf.height)
                            plot(sin)
                            ## ##mtext( caption.fun(), outer=TRUE, side=1,line=-3.5, cex=1, xpd=NA)
                            ## if(add.ADVERTISEMENT) {
                            ##     par(mfrow=c(1,1),par=c(1,1,1,1)*0.5)
                            ##     frame()
                            ##     motto <- pgx.randomSlogan(b=30)
                            ##     mex = min(pdf.width,pdf.height)/8.0
                            ##     text(0.5, 0.95, "Created with BigOmics Playground", cex=1.2*mex)
                            ##     text(0.5, 0.70, motto, cex=2*mex, font=4)
                            ##     text(0.5, 0.25, "BigOmics Analytics", cex=1.5*mex, font=2)
                            ##     text(0.5, 0.15, "Self-service bioinformatics solutions", cex=1.2*mex)
                            ##     text(0.5, 0.05, "www.bigomics.ch", cex=1.2*mex, font=3, xpd=NA)
                            ## }
                            dev.off()  ## important!!
                        }
                    } else { ## end base                
                        pdf(PDFFILE, pointsize=pdf.pointsize)
                        plot.new()
                        mtext("Error. PDF not available.",line=-8)
                        dev.off()
                    }

                    cat("file.exists(PDFFILE) = ",file.exists(PDFFILE),"\n")

                    ## ImageMagick or pdftk
                    if(WATERMARK) {
                        cat("adding watermark to PDF...\n")
                        ##PDFFILE="~/Downloads/plot.pdf"
                        tmp1 <- paste0(tempfile(),".pdf")
                        file.copy(PDFFILE,tmp1)
                        cmd = paste("pdftk",tmp1,"stamp ../lib/watermark.pdf output",PDFFILE)
                        cat("cmd = ",cmd,"\n")
                        system(cmd)
                    }
                    
                    ## finally copy to final exported file
                    file.copy(PDFFILE,file)
                    
                }, message="exporting to PDF", value=0.8)
            } ## content 
        ) ## PDF downloadHandler
    } ## end if do.pdf
    
    saveHTML <- function() {
        ## unlink(HTMLFILE) ## do not remove!
        if(plotlib == "plotly" ) {
            cat("downloadHandler:: exporting plotly to HTML\n")
            p <- func()
            if(WATERMARK) p <- addWatermarkPlotly(p)
            htmlwidgets::saveWidget(p, HTMLFILE) 
        } else if(plotlib %in% c("htmlwidget","pairsD3","scatterD3") ) {
            cat("downloadHandler:: exporting htmlwidget to HTML\n")
            p <- func()
            htmlwidgets::saveWidget(p, HTMLFILE) 
        } else if(plotlib == "iheatmapr") {
            p <- func()
            save_iheatmap(p, HTMLFILE)
        } else if(plotlib == "visnetwork") {
            p <- func()
            visSave(p, HTMLFILE)
        } else if(plotlib %in% c("ggplot","ggplot2")) {
            p <- func()
            ##ggsave(PDFFILE, width=pdf.width, height=pdf.height)
            saveWidget( ggplotly(p), file = HTMLFILE);
        } else if(plotlib=="generic") {
            cat("downloadHandler:: generic to HTMLFILE = ",HTMLFILE,"\n")
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
        download.html <- downloadHandler(
            filename = "plot.html",
            content = function(file) {
                withProgress({
                    ## unlink(HTMLFILE) ## do not remove!
                    if(plotlib == "plotly" ) {
                        cat("downloadHandler:: exporting plotly to HTML\n")
                        p <- func()
                        if(WATERMARK) p <- addWatermarkPlotly(p)
                        htmlwidgets::saveWidget(p, HTMLFILE) 
                    } else if(plotlib %in% c("htmlwidget","pairsD3","scatterD3") ) {
                        cat("downloadHandler:: exporting htmlwidget to HTML\n")
                        p <- func()
                        htmlwidgets::saveWidget(p, HTMLFILE) 
                    } else if(plotlib == "iheatmapr") {
                        p <- func()
                        save_iheatmap(p, HTMLFILE)
                    } else if(plotlib == "visnetwork") {
                        p <- func()
                        visSave(p, HTMLFILE)
                    } else if(plotlib %in% c("ggplot","ggplot2")) {
                        p <- func()
                        ##ggsave(PDFFILE, width=pdf.width, height=pdf.height)
                        saveWidget( ggplotly(p), file = HTMLFILE);
                    } else if(plotlib=="generic") {
                        cat("downloadHandler:: generic to HTMLFILE = ",HTMLFILE,"\n")
                        ## generic function should produce PDF inside plot func()
                        ##
                    } else if(plotlib=="base") {
                        write("<body>R base plots cannot export to HTML</body>",HTMLFILE)
                    } else { ## end base
                        write("<body>HTML export error</body>",file=HTMLFILE)
                    }
                    ## finally copy to fina lexport file
                    file.copy(HTMLFILE,file)
                }, message="exporting to HTML", value=0.8)
            } ## end of content
        ) ## end of HTML downloadHandler
    } ## end of do HTML

    ##if(do.csv && is.null(download.csv) )  {
    if(do.csv)  {
        download.csv <- downloadHandler(
            filename = "data.csv",
            content = function(file) {
                withProgress({
                    data <- csvFunc()                    
                    ##file.copy(CSVFILE,file)
                    write.csv(data, file=file)
                }, message="exporting to CSV", value=0.8)
            } ## end of content
        ) ## end of HTML downloadHandler
    } ## end of do HTML
    

    ##--------------------------------------------------------------------------------
    ##------------------------ OUTPUT ------------------------------------------------
    ##--------------------------------------------------------------------------------

    ## if("download.csv" %in% names(figure))  output$csv <- download.csv
    ## if("download.pdf" %in% names(figure))  output$pdf <- download.pdf
    ## if("download.png" %in% names(figure))  output$png  <- download.png
    ## if("download.html" %in% names(figure)) output$html <- download.html

    if(!is.null(download.csv))  output$csv  <- download.csv
    if(!is.null(download.pdf))  output$pdf  <- download.pdf
    if(!is.null(download.png))  output$png  <- download.png
    if(!is.null(download.html)) output$html <- download.html

    ## all widget IDs
    ##if(is.null(ns)) ns=function(x){x}
    ## buttons <- list(
    ##     info =  ns("info"),
    ##     options = ns("options"),
    ##     html =  ns("html"),
    ##     pdf = ns("pdf"),
    ##     png = ns("png"),
    ##     csv = ns("csv"),
    ##     zoom =  ns("zoom")
    ## )

    ##--------------------------------------------------------------------------------
    ##---------------------------- UI ------------------------------------------------
    ##--------------------------------------------------------------------------------

    ##func2 <- func  ## must be pryr  %<a-% to work!!!
    if(is.null(func2)) func2 <- func
    if(length(height)==1) height <- c(height,700)
    if(length(width)==1)  width  <- c(width,1200)
    if(length(res)==1)    res    <- c(res, 1.3*res)    

    res.1 <- res[1]
    res.2 <- res[2]
    ifnotchar.int <- function(s) ifelse(grepl("[%]|auto",s),s,as.integer(s))
    width.1 <- ifnotchar.int(width[1])
    width.2 <- ifnotchar.int(width[2])
    height.1 <- ifnotchar.int(height[1])
    height.2 <- ifnotchar.int(height[2])
    
    ##outputFunc <- sub(".*::","",outputFunc)
    render = render2 = NULL   
    if(1 && plotlib=="base") {
        ## For base plots we need to export PDF/PNG here. For some
        ## reason does not work in the downloadHandler...
        ##
        render <- renderPlot({
            func()
            if(1) {
                suppressWarnings( suppressMessages(
                    if(do.pdf) {
                        dev.print(pdf, file=PDFFILE, width=pdf.width, height=pdf.height,
                                  pointsize=pdf.pointsize)
                        ##dev.copy2pdf(file=PDFFILE, width=pdf.width, height=pdf.height)
                        ##dev.off()  ## important!!
                    }
                ))            
                suppressWarnings( suppressMessages(
                    if(do.png) {
                        dev.print(png, file=PNGFILE, width=pdf.width*100, height=pdf.height*100,
                                  pointsize=pdf.pointsize)
                        ##dev.copy2pdf(file=PDFFILE, width=pdf.width, height=pdf.height)
                        ##dev.off()  ## important!!
                    }
                ))
            }            
            ##}, res=res[1], width=width[1], height=height[1] )
            ##}, res=res[1], width=width.1, height=height.1 )
            ##}, width=width.1, height=height.1, res=res.1 )
            }, res=res.1)
        ##})
        if(!is.null(func2)) {
            ##render2 <- renderPlot(func2(), res=res[1], width=width[2], height=height[2])
            ##render2 <- renderPlot(func2(), res=res[1], width=width.2, height=height.2)
            ##render2 <- renderPlot(func2(), width=width.2, height=height.2, res=res.2)
            render2 <- renderPlot(func2(), res=res.2)            
            ##render2 <- renderPlot(func2())
        }
    } else if(plotlib=="image") {
        render <- renderImage( func(), deleteFile=FALSE)
        if(!is.null(func2)) {
            render2 <- renderImage(func2(), deleteFile=FALSE)
        }
    } else {
        render <- eval(parse(text=renderFunc))(func())
        if(!is.null(func2)) {
            render2 <- eval(parse(text=renderFunc))(func2())
        }
    }
    
    output$renderfigure <- render
    output$renderpopup  <- render2

    output$popupfig <- renderUI({
        w <- width.2
        h <- height.2
        if(plotlib=="image") {
            ## retains aspect ratio
            ##
            img.file <- func()$src
            require(png)
            require(jpeg)
            img.dim <- c(h,w)
            if(grepl("png|PNG",img.file)) img.dim <- dim(readPNG(img.file))[1:2]
            if(grepl("jpg|JPG",img.file)) img.dim <- dim(readJPEG(img.file))[1:2]
            r <- min( width.2 / img.dim[2], height.2 / img.dim[1])
            h <- img.dim[1]*r
            w <- img.dim[2]*r
        } 
        if(!is.null(render2)) {
            r <- eval(parse(text=outputFunc))(ns("renderpopup"), width=w, height=h)
        } else {
            r <- eval(parse(text=outputFunc))(ns("renderfigure"), width=w, height=h)
        }
  
        if(any(class(caption2)=="reactive")) {
            caption2 <- caption2()
        }
        if(any(class(caption2)=="character")) {
            caption2 <- HTML(caption2)
        }
        tagList(
            r,
            div( caption2, class="caption2")          
        )
    })
        
    output$widget <- renderUI({
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
            caption <- HTML(caption)
        }

        fillCol(
            flex = c(NA,NA,1,NA,NA,0.001),
            height = height.1,
            tagList(
                tags$head(tags$style(modaldialog.style)),
                tags$head(tags$style(modalbody.style)),            
                tags$head(tags$style(modalfooter.none))
                ##tags$head(tags$style(".caption { margin: 20px 5px;}"))
            ),
            ##uiOutput(ns("renderbuttons")),
            buttons,
            ##render,
            eval(parse(text=outputFunc))(ns("renderfigure"), width=width.1, height=height.1),
            br(),
            div( caption, class="caption"),          
            div(class="popup-plot",
                bsModal(ns("plotPopup"), title, size="l",
                        ns("zoombutton"),
                        ##tagList(renderPlot(plot(sin)))
                        tagList(uiOutput(ns("popupfig")))
                        )
                )
        )
    })
    outputOptions(output, "widget", suspendWhenHidden=FALSE) ## important!!!
    
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
        outputFunc = outputFunc,
        renderFunc = renderFunc
    )
    return(res)
}

##================================================================================
##================================================================================
##================================================================================


tableWidget <- function(id) {
    ns <- NS(id)
    uiOutput(ns("widget"))
}

tableModule <- function(input, output, session, 
                        func, func2=NULL, info.text="Info text",
                        title=NULL, label=NULL, caption=NULL, 
                        server=TRUE, filename="data.csv", ##inputs=NULL, 
                        ##no.download = FALSE, just.info=FALSE,
                        width="100%", height="auto",
                        options = NULL, info.width="300px"
                        )
{
    require(shinyWidgets)
    require(shinyBS)

    ns <- session$ns
    
    if(any(class(caption)=="reactive")) {
        caption <- caption()
    }
    if(class(caption)=="character") {
        caption <- HTML(caption)
    }
    
    options.button <- ""    
    if(!is.null(options) && length(options)>0) {
        options.button <- dropdownButton(
            options,
            ##br(),
            ##dload,
            circle = TRUE, size = "xs", ## status = "danger",
            ## icon = icon("gear"),
            icon = icon("bars"),
            width = "250px",
            inputId = ns("options"),
            tooltip = tooltipOptions(title = "Settings", placement = "right")
        )
    }
    label1 = HTML(paste0("<span class='module-label'>",label,"</span>"))

    require(shinyWidgets)
    zoom.button <- actionButton(inputId=ns("zoombutton"),label=NULL,
                                icon=icon("window-maximize"),
                                class="btn-circle-xs")    
    
    buttons <- fillRow(
        ##flex=c(NA,NA,NA,NA,1),
        flex=c(NA,NA,NA,NA,NA,1),
        ##actionLink(options_id, label=NULL, icon = icon("info")),
        label1,
        dropdownButton(
            tags$p(HTML(info.text)),
            br(),
            circle = TRUE, size = "xs", ## status = "danger",
            icon = icon("info"), width = info.width,
            inputId = ns("info"), right=FALSE,
            tooltip = tooltipOptions(title = "Info", placement = "right")
        ),
        options.button,
        div(class='download-button',
            dropdownButton(
                downloadButton(ns("csv"), "CSV"),
                circle = TRUE, size = "xs", ## status = "danger",
                icon = icon("download"), width = "80px", right=FALSE,
                tooltip = tooltipOptions(title = "Download", placement = "right")
            )),
        zoom.button,
        HTML(paste("<center>",title,"</center>"))
        ##HTML(paste("<center><strong>",title,"</strong></center>"))
        ## HTML(paste("<center>",title,"</center>"))
        ##inputs
        ##selectInput("sel123","number",1:10)
    )
    
    CSVFILE = paste0(gsub("file","data",tempfile()),".csv")
    CSVFILE    
    
    ## render2 <- renderPlot({plot_array[[3]]()}, res=res)
    download.csv <- downloadHandler(
        filename = filename,
        content = function(file) {
            dt <- func()
            write.csv(dt$x$data, file=CSVFILE, row.names=FALSE)
            file.copy(CSVFILE,file)
        }
    )
    output$csv <- download.csv
    
    if(is.null(func2)) func2 <- func
    if(length(height)==1) height <- c(height,700)
    if(length(width)==1)  width  <- c(width,1280)
    ifnotchar.int <- function(s) ifelse(grepl("[%]|auto",s),s,as.integer(s))
    width.1 <- ifnotchar.int(width[1])
    width.2 <- ifnotchar.int(width[2])
    height.1 <- ifnotchar.int(height[1])
    height.2 <- ifnotchar.int(height[2])

    output$datatable <- renderDataTable({
        dt <- func()
        ##dt <- dt %>% DT::formatStyle(0, target='row', fontSize='11px', lineHeight='80%')
        dt
    })
    output$datatable2 <- renderDataTable({
        func2() %>% DT::formatStyle(0, target='row', fontSize='12px', lineHeight='90%')
    })
    
    output$popuptable <- renderUI({
        dataTableOutput(ns("datatable2"), width=width.2-40, height=height.2-40)
    })
    
    output$widget <- renderUI({

        ##modaldialog.style <- paste0("#",ns("tablePopup")," .modal-dialog {width:",width.2+40,"px;}")
        ##modalbody.style <- paste0("#",ns("tablePopup")," .modal-body {min-height:",height.2+40,"px;}")
        modaldialog.style <- paste0("#",ns("tablePopup")," .modal-dialog {width:",width.2,"px;}")
        modalbody.style <- paste0("#",ns("tablePopup")," .modal-body {min-height:",height.2,"px;}")
        modalfooter.none <- paste0("#",ns("tablePopup")," .modal-footer{display:none;}")
        div.caption <- NULL
        if(!is.null(caption)) div.caption <- div(caption, class="table-caption")
        
        fillCol(
            flex = c(NA,NA,1,NA),
            ##tags$head(tags$style(".popup-table .modal-dialog {width:1200px;}")),
            ##tags$head(tags$style(".popup-table .modal-body {min-height:700px;}")),
            ##tags$head(tags$style(".popup-table .modal-footer {display:none;}")),            
            tags$head(tags$style(modaldialog.style)),
            tags$head(tags$style(modalbody.style)),
            tags$head(tags$style(modalfooter.none)),
            buttons,
            ##uiOutput(ns("renderbuttons")),
            div.caption,
            dataTableOutput(ns("datatable"), width=width.1, height=height.1),
            ##DT::renderDataTable(func()),
            ##verbatimTextOutput(ns("zoomstatus"))
            div(class="popup-table",
                bsModal(ns("tablePopup"), title, size="l",
                        ns("zoombutton"),
                        ##tagList(renderPlot(plot(sin)))
                        tagList(uiOutput(ns("popuptable")))
                        )
                )
            
        )
    })
    ##outputOptions(output, "widget", suspendWhenHidden=FALSE) ## important!!!

    module <- list(
        rows_selected = reactive(input$datatable_rows_selected),
        rows_all = reactive(input$datatable_rows_all)
    )
    return(module)
}





########################################################################
########################################################################
########################################################################
