########################################################################
## PLOT/TABLE MODULES
########################################################################

##if(!exists("WATERMARK")) WATERMARK=FALSE
##WATERMARK=TRUE
##WATERMARK=FALSE

SLOGAN = c(
    "'I Feel Empowered' - with OmicsPlayground",
    "'I Love Data' - with OmicsPlayground",
    "'I Get Insights' - with OmicsPlayground",
    "'Data, Knowledge, Insight' - with OmicsPlayground",
    "'So Much More' - with OmicsPlayground",
    "'Do-it-myself' - with OmicsPlayground",
    "'Yes-I-Can' - with OmicsPlayground",
    "'Dig Deeper' - with OmicsPlayground",
    "'Never Stop Exploring' - with OmicsPlayground",
    "'Take control' - with OmicsPlayground",    
    "'My Eureka! moments' - with OmicsPlayground",
    "'I See Clearly Now' - with OmicsPlayground"
    )

addWatermark.base <- function(line=-1.5, cex=0.8, col="#6699ff88") {    
    txt = "'Feel empowered' --- with Omics Playground"
    txt = sample(SLOGAN,1)
    ##mtext(txt,outer=TRUE, side=1,line=line, cex=cex, col=col, xpd=NA)
    title(sub=txt, adj=1, line=3.5, font=2, cex.sub=cex, col.sub=col)
}
addWatermark.Plotly <- function(p) {
    add_annotations(
        p,
        ##yshift = -100, 
        ##x = 1, y=-0.1, xanchor = "right",
        ##x = 0.5, y=-0.1, xanchor = "middle", 
        x = 0.5, y = 0.5, xanchor = "middle", 
        xref="paper", yref="paper", showarrow=FALSE,
        ##text = "'Feel empowered' --- with Omics Playground",
        ##text = "'Feel empowered'\nwith Omics Playground",
        text = sub(" - ","\n",sample(SLOGAN,1)),
        font = list(size=50, family='Arial', color="#AAAAAA22")
        ##font = list(size=18, family='Arial', color="#AAAAAA22")
    )
}

colBL="#00448855"
colRD="#88004455"
addWatermark.PDF.SAVE <- function(file, col="#88006655") {
    if(system("which pdftk",ignore.stdout=TRUE)==1) return ## if no pdftk installed...
    tmp1 <- paste0(gsub("file","plot",tempfile()),".pdf")
    tmp2 <- paste0(gsub("file","plot",tempfile()),".pdf")
    txt <- sample(SLOGAN,1)
    ##txt = sub("with Omics","with\nOmics",sub(" - ","\n",txt))
    pdf(tmp1,w=8,h=8)
    frame()
    ##text(0.5,0.5,txt,font=2,cex=4,col=col)
    legend("bottom", txt, cex=1.75, text.font=2, text.col=col,
           xjust = 0.5, yjust = 0.5, adj=-0.012,
           inset=c(0,-0.15), xpd=TRUE,
           bg="#88006611", border=NULL, box.col = NA,
           text.width = 2.05*strwidth(txt) )
    ##text(0.5,1.0,txt,font=2,cex=1.5,col=col,xpd=FALSE)
    ##text(1.08,0.5,txt,srt=90,font=1,cex=1.3,col=col,xpd=TRUE)
    ##text(0.5,-0.2,txt,srt=0,font=2,cex=1,col=col,xpd=TRUE)
    dev.off()
    cmd <- paste("pdftk",file,"stamp",tmp1,"output",tmp2) ## NEED pdftk installed!!!
    cmd
    system(cmd)
    file.copy(tmp2,file,overwrite=TRUE)
    unlink(tmp1)
    unlink(tmp2)
}
addWatermark.PDF <- function(file) {
    if(system("which pdftk",ignore.stdout=TRUE)==1) return ## if no pdftk installed...
    mark <- file.path(FILES,"watermark.pdf")
    tmp <- paste0(gsub("file","plot",tempfile()),".pdf")
    cmd <- paste("pdftk",file,"stamp",mark,"output",tmp) ## NEED pdftk installed!!!
    cmd
    system(cmd)
    file.copy(tmp,file,overwrite=TRUE)
    unlink(tmp)
}

if(0) {
    file = "/home/kwee/Downloads/plot.pdf"
    addWatermark.PDF(file) 
}
## prepare ORCA server
library(plotly)
if(!exists("ORCA") || !ORCA$process$is_alive()) {
    ORCA <- orca_serve(
        more_args="--enable-webgl",
        env = c(DISPLAY=":0")
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
    require(plotly)
    
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

                        p <- func()
                        p$width = pdf.width * 80
                        p$height = pdf.height * 80

                        cat("[pgx-modules::plotModule] exporting plotly to PNG with ORCA\n")
                        ##err <- try(export(p, PNGFILE))  ## deprecated 
                        ##err <- try(orca(p, PNGFILE))
                        err <- try(ORCA$export(p, PNGFILE, width=p$width, height=p$height)) 
                        if(class(err)!="try-error") {
                            cat("[pgx-modules::plotModule] OK!\n")
                        } else {
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
                            cat("downloadHandler:: creating new PNG device\n")
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
                    file.copy(PNGFILE,file)
                    cat("[pgx-modules::plotModule] export to PNG done!\n")                    
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
                        ## if(WATERMARK) p <- addWatermark.Plotly(p)
                        p$width = pdf.width * 80
                        p$height = pdf.height * 80

                        cat("[pgx-modules::plotModule] exporting plotly to PDF with ORCA\n")
                        ##err <- try(export(p, PDFFILE))  ## deprecated
                        ##err <- try(orca(p, PDFFILE))
                        err <- try(ORCA$export(p, PDFFILE, width=p$width, height=p$height)) 
                        if(class(err)!="try-error") {
                            cat("[pgx-modules::plotModule] OK!\n")                            
                        } else {                            
                            cat("[pgx-modules::plotModule] export failed, trying webshot...\n")
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
                            pdf(file=PDFFILE, width=pdf.width, height=pdf.height, pointsize=pdf.pointsize)
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
                    file.copy(PDFFILE,file)

                    ## ImageMagick or pdftk
                    if(TRUE && WATERMARK) {
                        cat("[pgx-modules::plotModule] adding watermark to PDF...\n")
                        ##addWatermark.PDF(PDFFILE)
                        addWatermark.PDF(file) 
                    }
                    cat("[pgx-modules::plotModule] export to PDF done!\n")


                }, message="exporting to PDF", value=0.8)
            } ## content 
        ) ## PDF downloadHandler
    } ## end if do.pdf
    
    saveHTML <- function() {
        ## unlink(HTMLFILE) ## do not remove!
        if(plotlib == "plotly" ) {
            cat("downloadHandler:: exporting plotly to HTML\n")
            p <- func()
            if(WATERMARK) p <- addWatermark.Plotly(p)
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
                        if(WATERMARK) p <- addWatermark.Plotly(p)
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

    if(!is.null(download.csv))  output$csv  <- download.csv
    if(!is.null(download.pdf))  output$pdf  <- download.pdf
    if(!is.null(download.png))  output$png  <- download.png
    if(!is.null(download.html)) output$html <- download.html

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
                        dev.print(pdf, file=PDFFILE, width=pdf.width, height=pdf.height,
                                  pointsize=pdf.pointsize)
                        ##p1.base
                        ##if(WATERMARK) addWatermark.base()
                        ##dev.off()  ## important!!
                    }
                ))
            }
            if(1) {
                suppressWarnings( suppressMessages(
                    if(do.png) {
                        ##png(file=PNGFILE, width=pdf.width*100, height=pdf.height*100,pointsize=pdf.pointsize)
                        dev.print(png, file=PNGFILE, width=pdf.width*100, height=pdf.height*100,
                                  pointsize=pdf.pointsize)
                        ##p1.base
                        ##dev.off()  ## important!!
                    }
                ))
            }
            ##p1.base
        }, res=res.1)

        if(!is.null(func2)) {
            render2 <- renderPlot({
                func2()
                ##if(WATERMARK) addWatermark.base()                
            }, res=res.2)            
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
