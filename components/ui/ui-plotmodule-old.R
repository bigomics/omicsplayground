##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

########################################################################
## PLOT/TABLE MODULES
########################################################################

##if(!exists("WATERMARK")) WATERMARK=FALSE
##WATERMARK=TRUE
##WATERMARK=FALSE

colBL="#00448855"
colRD="#88004455"
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

addWatermark.PNG <- function(file) {
    if(system("which convert",ignore.stdout=TRUE)==1) return ## if no pdftk installed...
    tmp <- paste0(gsub("file","plot",tempfile()),".png")
    cmd = "convert plot.png -font Helvetica -pointsize 13 -extent 100%x105% -draw \"gravity south fill #80000080 text 0,4 'Created using the OmicsPlayground. Developed by BigOmics Analytics in Switzerland.' \"  plot_wmark.png"
    cmd <- sub("plot.png",file,cmd)
    cmd <- sub("plot_wmark.png",tmp,cmd)
    system(cmd)
    file.copy(tmp,file,overwrite=TRUE)
    unlink(tmp)
}

##================================================================================
##=================== Plotly Export function =====================================
##================================================================================

##;format="pdf";width=height=800;scale=1;file="plot.pdf";server=NULL
##p=plot_ly(1:10)
plotlyExport <- function(p, file = "plot.pdf", format = tools::file_ext(file),
                         scale = NULL, width = NULL, height = NULL, server=NULL)
{
    is.docker <- file.exists("/.dockerenv")
    is.docker
    export.ok <- FALSE

    if(class(p)[1] != "plotly") {
        message("[plotlyExport] ERROR : not a plotly object")
        return(NULL)
    }
    ## remove old
    unlink(file,force=TRUE)

    ## See if Kaleido is available
    if(1 && !export.ok) {
        ## https://github.com/plotly/plotly.R/issues/2179
        reticulate::py_run_string("import sys")
        err <- try(suppressMessages(plotly::save_image(p, file=file, width=width, height=height, scale=3)))
        export.ok <- class(err)!="try-error"
        if(export.ok) message("[plotlyExport] --> exported with plotly::save_image() (kaleido)")
        export.ok <- TRUE
    }
    if(1 && !export.ok) {
        ## works only for non-GL plots
        err <- try(plotly::export(p, file, width=width, height=height))
        export.ok <- class(err)!="try-error"
        if(export.ok) message("[plotlyExport] --> exported with plotly::export() (deprecated)")
    }
    if(0 && !export.ok) {
        tmp = paste0(tempfile(),".html")
        htmlwidgets::saveWidget(p, tmp)
        err <- try(webshot::webshot(url=tmp,file=file,vwidth=width*100, vheight=height*100))
        export.ok <- class(err)!="try-error"
        if(export.ok) message("[plotlyExport] --> exported with webshot::webshot()")
    }
    if(!export.ok) {
        message("[plotlyExport] WARNING: export failed!")
        if(format=="png") png(file)
        if(format=="pdf") pdf(file)
        par(mfrow=c(1,1));frame()
        text(0.5,0.5,"Plotly export error",cex=2)
        dev.off()
    }

    message("[plotlyExport] file.exists(file)=",file.exists(file))
    export.ok <- export.ok && file.exists(file)
    return(export.ok)
}

##================================================================================
##================================================================================
##================================================================================

plotWidget <- function(id,...) {
    ns <- shiny::NS(id)
    shiny::uiOutput(ns("widget"),...)
}

plotModule <- function(input, output, session,
                       func, func2=NULL,
                       info.text="Figure", title="",
                       inputs=NULL, options = NULL, label="",
                       caption="", caption2=info.text, ## header=NULL,
                       plotlib = "base", plotlib2 = NULL,
                       renderFunc=NULL, outputFunc=NULL, csvFunc=NULL,
                       renderFunc2=NULL, outputFunc2=NULL,
                       no.download = FALSE, download.fmt=c("png","pdf"),
                       just.info=FALSE, info.width="300px", show.maximize = TRUE,
                       height = c(640,800), width = c("auto",1400), res=c(72,100),
                       download.pdf = NULL, download.png = NULL,
                       download.html = NULL, download.csv = NULL,
                       pdf.width=8, pdf.height=6, pdf.pointsize=12,
                       add.watermark=FALSE )
{
    ns <- session$ns

    ##--------------------------------------------------------------------------------
    ##------------------------ BUTTONS -----------------------------------------------
    ##--------------------------------------------------------------------------------

    if(is.null(inputs) || length(inputs)==0 ) inputs <- ""
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
    ##if("csv" %in% download.fmt)  dload.csv  <- shiny::downloadButton(ns("csv"), "CSV")
    if(!is.null(csvFunc))  dload.csv  <- shiny::downloadButton(ns("csv"), "CSV")

    pdf_size = NULL
    if(plotlib!="base") {
        pdf_size <- shiny::tagList(
            shiny::fillRow(
                shiny::numericInput(ns("pdf_width"), "width", pdf.width, 1, 20, 1, width='100%'),
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
    ##label1 = shiny::HTML(paste0("<span class='module-label'>",label,"</span>"))
    title1 <- title
    if(label!="") {
        title1 = shiny::HTML(paste0(title," (",label,")"))
    }

    ##zoom.button <- shinyWidgets::prettyCheckbox(inputId=ns("zoom"),label=NULL,value=FALSE)
    zoom.button <- NULL
    if(show.maximize) {
        zoom.button <- modalTrigger(ns("zoombutton"), ns("plotpopup"),
            icon("window-maximize"), class="btn-circle-xs")
        zoom.button <- withTooltip(zoom.button, "Maximize plot", placement="right")
    }


    ##output$renderbuttons <- shiny::renderUI({
    ## button layout
    header <- shiny::fillRow(
        flex = c(1,NA,NA,NA,NA),
        ##label1,
        ##shiny::HTML(title1),
        shiny::div(class='plotmodule-title', title=title, title1),
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
        shiny::div(class='zoom-button', zoom.button)
        ##),
    )
    ##return(ui)
    ##})

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
                    data <- csvFunc()
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
    ifnotchar.int <- function(s) ifelse(grepl("[%]|auto",s),s,as.integer(s))
    width.1  <- ifnotchar.int(width[1])
    width.2  <- ifnotchar.int(width[2])
    height.1 <- ifnotchar.int(height[1])
    height.2 <- ifnotchar.int(height[2])

    ## This sets the correct render and output functions for different
    ##plotting libraries.
    getOutputRenderFunc <- function(plotlib, outputFunc, renderFunc)
    {
        if(plotlib == "generic") {
            if(is.null(renderFunc)) stop("'generic' class must provide renderFunc")
            if(is.null(outputFunc)) stop("'generic' class must provide outputFunc")
        } else if(plotlib == "htmlwidget") {
            if(is.null(renderFunc)) stop("'htmlwidget' class must provide renderFunc")
            if(is.null(outputFunc)) stop("'htmlwidget' class must provide outputFunc")
        } else if(plotlib == "plotly") {
            ##render <- plotly::renderPlotly({ func() })

            if(is.null(outputFunc)) outputFunc = "plotly::plotlyOutput"
            if(is.null(renderFunc)) renderFunc = "plotly::renderPlotly"
        } else if(plotlib == "echarts") {
            ##render <- plotly::renderPlotly({ func() })

            if(is.null(outputFunc)) outputFunc = "echarts4r::echarts4rOutput"
            if(is.null(renderFunc)) renderFunc = "echarts4r::renderEcharts4r"
        } else if(plotlib=="scatterD3") {

            if(is.null(renderFunc)) renderFunc="scatterD3::renderScatterD3"
            if(is.null(outputFunc)) outputFunc="scatterD3::scatterD3Output"
        } else if(plotlib=="pairsD3") {
            if(is.null(renderFunc)) renderFunc="pairsD3::renderPairsD3"
            if(is.null(outputFunc)) outputFunc="pairsD3::pairsD3Output"
        } else if(plotlib == "visnetwork") {
            if(is.null(outputFunc)) outputFunc="visNetwork::visNetworkOutput"
            if(is.null(renderFunc)) renderFunc = "visNetwork::renderVisNetwork"
        } else if(plotlib %in% c("ggplot","ggplot2")) {
            if(is.null(outputFunc)) outputFunc="shiny::plotOutput"
            if(is.null(renderFunc)) renderFunc="function(x) shiny::renderPlot(plot(x))"
        } else if(plotlib %in% c("grid")) {
            if(is.null(outputFunc)) outputFunc="shiny::plotOutput"
            if(is.null(renderFunc)) renderFunc="function(x) shiny::renderPlot(grid::grid.draw(x, recording=FALSE))"
        } else if(plotlib == "iheatmapr") {
            if(is.null(outputFunc)) outputFunc="iheatmapr::iheatmaprOutput"
            if(is.null(renderFunc)) renderFunc="iheatmapr::renderIheatmap"
        } else if(plotlib == "image") {
            if(is.null(outputFunc)) outputFunc="shiny::imageOutput"
            if(is.null(renderFunc)) renderFunc="shiny::renderImage"
        } else {
            ## Base plotting
            ## render <- shiny::renderPlot({
            ##     par(mar=c(0,0,0,0),oma=c(0,0,0,0))
            ##     func()
            ## }, res=res, width=width, height=height)
            if(is.null(outputFunc)) outputFunc="shiny::plotOutput"
            ##renderFunc="renderPlot"
            ##renderFunc="function(x) shiny::renderPlot(x, res=res, width=width, height=height)"
            if(is.null(renderFunc)) renderFunc="shiny::renderPlot"
        }
        list( outputFunc=outputFunc, renderFunc=renderFunc )
    }

    if(is.null(outputFunc) || is.null(renderFunc)) {
        out1 <- getOutputRenderFunc(plotlib, outputFunc, renderFunc)
        outputFunc = out1$outputFunc
        renderFunc = out1$renderFunc
    }

    if((is.null(outputFunc2) || is.null(renderFunc2)) && plotlib2 != plotlib ) {
        out2 <- getOutputRenderFunc(plotlib2, outputFunc2, renderFunc2)
        outputFunc2 = out2$outputFunc
        renderFunc2 = out2$renderFunc
    }

    if((is.null(outputFunc2) || is.null(renderFunc2)) && plotlib2 == plotlib ) {
        if(is.null(outputFunc2) && !is.null(outputFunc)) outputFunc2 = outputFunc
        if(is.null(renderFunc2) && !is.null(renderFunc)) renderFunc2 = renderFunc
    }

    if(0) {
        dbg("[plotModule] title=",title)
        dbg("[plotModule] plotlib=",plotlib)
        dbg("[plotModule] plotlib2=",plotlib2)
        dbg("[plotModule] outputFunc=",outputFunc)
        dbg("[plotModule] outputFunc2=",outputFunc2)
        dbg("[plotModule] renderFunc=",renderFunc)
        dbg("[plotModule] renderFunc2=",renderFunc2)
        dbg("[plotModule] is.null.renderFunc=",is.null(renderFunc))
        dbg("[plotModule] is.null.renderFunc2=",is.null(renderFunc2))
    }

    ##outputFunc <- sub(".*::","",outputFunc)
    render = render2 = NULL
    if(1 && plotlib=="base") {
        ## For base plots we need to export PDF/PNG here. For some
        ## reason does not work in the downloadHandler...
        ##
        render <- shiny::renderPlot({
            pdf.width0  <- pdf.width
            pdf.height0 <- pdf.height

            ## save base render in object
            func()
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
    if(grepl("renderCachedPlot",renderFunc)) {
        render <- shiny::renderCachedPlot(
            func(), cacheKeyExpr = {list(csvFunc())}, res=res.1
        )
    }
    if(grepl("renderCachedPlot",renderFunc2)) {
        render2 <- shiny::renderCachedPlot(
            func2(), cacheKeyExpr = {list(csvFunc())}, res=res.2
        )
    }

    if(is.null(render)) {
        render <- eval(parse(text=renderFunc))(func())
    }
    if(is.null(render2) && !is.null(func2)) {
        render2 <- eval(parse(text=renderFunc2))(func2())
    }

    output$renderfigure <- render
    output$renderpopup  <- render2

    output$popupfig <- shiny::renderUI({
        w <- width.2
        h <- height.2
        if(plotlib2=="image") {
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
        r <- eval(parse(text=outputFunc2))(ns("renderpopup"), width=w, height=h)

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
    })

    output$widget <- shiny::renderUI({

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

        div(
            class = "plotmodule",
            shiny::fillCol(
                flex = c(NA,NA,1,NA,NA,0.001),
                height = height.1,
                shiny::tagList(
                    shiny::tags$head(shiny::tags$style(modaldialog.style)),
                    shiny::tags$head(shiny::tags$style(modalbody.style)),
                    shiny::tags$head(shiny::tags$style(modalfooter.none))
                ),
                div(class="plotmodule-header", header),
                ##render,
                eval(parse(text=outputFunc))(ns("renderfigure"), width=width.1, height=height.1),
                shiny::br(),
                shiny::div(caption, class="caption"),
                shiny::div(class="popup-plot",
                    modalUI(ns("plotPopup"), title, size="fullscreen",
                        shiny::uiOutput(ns("popupfig"))
                    )
                )
            )
        )

    })
    shiny::outputOptions(output, "widget", suspendWhenHidden=FALSE) ## important!!!

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
        ##getCaption = caption.fun,
        saveHTML = saveHTML,
        outputFunc = outputFunc,
        renderFunc = renderFunc
    )
    return(res)
}
