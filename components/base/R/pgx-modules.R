##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
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
##============ Welcome to the ORCA hell .... =====================================
##================================================================================

## prepare ORCA server if we are not in a Docker


initOrca <- function(launch=TRUE) {

    responding.local = FALSE
    responding.docker = FALSE
    responding.process = FALSE
    orca.server = NULL
    
    ## See if an orca-server docker is running --------------------------------------
    srv.list <- c("http://orca-server:9091","http://localhost:9091")
    srv.responding <- sapply(srv.list, function(s)
        class(try(httr::POST(s, body=plotly:::to_JSON("")),silent=TRUE))!="try-error")
    srv.responding
    orca.server <- names(srv.responding)[which(srv.responding)]
    if(is.na(orca.server[1])) orca.server <- NULL
    responding.docker  <- !is.null(orca.server)
    ##message("docker ORCA response = ",class(res.docker))
    message("[initOrca] dockerized ORCA is responding = ",responding.docker)
    message("[initOrca] dockerized ORCA server = ", paste(orca.server,collapse=" "))
    
    ## If not, we try to connect to local orca-server at port 5151 -----------------
    if(!responding.docker) {
        res.local <- try(httr::POST("http://localhost:5151", body=plotly:::to_JSON("")),silent=TRUE)
        responding.local <- class(res.local)=="response"
        message("[initOrca] local ORCA server at http://localhost:5151")
        message("[initOrca] local ORCA server is responding = ",responding.local)
        orca.server = "http://localhost:5151"
    }
    
    ## If a orca docker or server is not running, we try to lauch a local orca process
    ORCA = NULL
    if(launch && !responding.docker && !responding.local) {
        message("[initOrca] Trying to launch local orca-server...")
        assignInNamespace("correct_orca", function() return(TRUE), ns="plotly")
        ##ORCA <- plotly::orca_serve(port=5151)
        ORCA <- plotly::orca_serve(port=5151, keep_alive=TRUE, more_args="--enable-webgl")
        ORCA$process$is_alive()
        for(i in 1:10) {
            res.process <- try(httr::POST("http://localhost:5151", body=plotly:::to_JSON("")),silent=TRUE)
            responding.process <- class(res.process)=="response"
            ##message("local ORCA is responding = ",responding.local)
            message("trying to reach local orca...")
            if(responding.process) {
                break
            }
            Sys.sleep(2)
        }
        res.process  <- try(httr::POST("http://localhost:5151", body=plotly:::to_JSON("")),silent=TRUE)
        has.response <- class(res.process)=="response"
        responding.process <- ORCA$process$is_alive() && has.response
        message("[initOrca] ORCA process is alive = ",ORCA$process$is_alive())
        message("[initOrca] ORCA process has response = ",has.response)
        message("[initOrca] ORCA process is responding = ",responding.process)
        orca.server = ORCA
    }
    
    if(responding.docker) {
        message("[initOrca] --> using dockerized ORCA server at ",orca.server)
    }

    if(responding.local) {
        message("[initOrca] --> using local ORCA server at http://localhost:5151")
    }

    if(responding.process) {
        message("[initOrca] --> using local ORCA process from Plotly")
    }

    ##orca.available <- responding.docker || responding.local || responding.process
    orca.available <- !is.null(orca.server)
    if(!orca.available) {
        message("[initOrca] WARNING!!! Could not connect to ORCA server")
    }    

    return(orca.server)
}


##;format="pdf";width=height=800;scale=1;file="plot.pdf";server=NULL
plotlyExport <- function(p, file = "plot.pdf", format = tools::file_ext(file), 
                         scale = NULL, width = NULL, height = NULL, server=NULL)
{
    is.docker <- file.exists("/.dockerenv")
    has.orca.bin <- file.exists("/usr/local/bin/orca")
    has.orca.export <- has.orca.bin && exists("ORCA") && "export" %in% names(ORCA)
    is.docker
    has.orca.bin
    has.orca.export
    export.ok <- FALSE
    
    message("[plotlyExport] class(p) = ",class(p)[1])
    if(class(p)[1] != "plotly") {
        message("[plotlyExport] ERROR : not a plotly object")
        return(NULL)
    }

    message("[plotlyExport] file = ",file)
    message("[plotlyExport] format = ",format)
    message("[plotlyExport] is.docker = ",is.docker)
    message("[plotlyExport] exists(ORCA) = ",exists("ORCA"))
    message("[plotlyExport] class(ORCA) = ",class(ORCA))
    ##message("[plotlyExport] ORCA$port = ",ORCA$port)
    ##message("[plotlyExport] ORCA is alive = ",ORCA$process$is_alive())
    
    ## remove old
    unlink(file,force=TRUE)

    ## See if Kaleido is available
    if(1 && !export.ok) {
        err <- try( plotly::save_image(p, file=file, width=width, height=height) )
        export.ok <- class(err)!="try-error"
        if(export.ok) message("[plotlyExport] --> exported with plotly::save_image() (kaleido)")
        export.ok <- TRUE
    }
    
    ## See if any ORCA server is responding (docker or already local)
    if(1 && !export.ok) {
        global.srv  <- exists("ORCA") && class(ORCA)[1]=="character"
        global.srv
        if(is.null(server) && !global.srv) {
            server <- c("http://orca-server:9091","http://localhost:9091","http://localhost:5151")
        }
        if(is.null(server) && global.srv) {
            server <- ORCA
        }
        srv.responding <- sapply(server, function(s)
            class(try(httr::POST(s, body=plotly:::to_JSON("")),silent=TRUE))!="try-error")
        srv.responding
        server <- names(srv.responding)[which(srv.responding)][1]
        if(length(server)==0 || is.na(server[1])) server <- NULL
        server
        message("[plotlyExport] using orca server = ",server)
        
        if(!is.null(server)) {
            bod <- list(figure = plotly::plotly_build(p)$x[c("data", "layout")],
                        format = format, width = width, height = height, 
                        scale = scale)
            res <- httr::POST(server, body = plotly:::to_JSON(bod))
            httr::stop_for_status(res)
            httr::warn_for_status(res)
            con <- httr::content(res, as = "raw")
            writeBin(con, file)
            message("[plotlyExport] --> exported with ORCA server on ",server)
            export.ok <- TRUE
        }
    }
    
    ## if(0 && has.orca.bin && !export.ok) {
    ##     ## BUG: currently plotly::orca() cannot write to /tmp folder
    ##     err <- try(plotly::orca(p, file=file, format=format, width=width, height=height))
    ##     export.ok <- class(err)!="try-error"
    ##     if(export.ok) message("[plotlyExport] --> exported with plotly::orca()")
    ## }
    
    if(1 && has.orca.export && !export.ok) {
        err <- try(ORCA$export(p, file=file, format=format, width=width, height=height))
        export.ok <- class(err)!="try-error"
        if(export.ok) message("[plotlyExport] --> exported with ORCA$export()")
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
        zoom.button <- modalTrigger(ns("zoombutton"),
            ns("plotpopup"),
                                    icon("window-maximize"),
                                    class="btn-circle-xs")
        zoom.button <- withTooltip(zoom.button, "Maximize", placement="right")
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
            eval(parse(text=outputFunc))(ns("renderfigure"), width=width.1, height=height.1),
            shiny::br(),
            shiny::div(caption, class="caption"),          
            shiny::div(class="popup-plot",
                modalUI(ns("plotPopup"), title, size="fullscreen",
                        shiny::uiOutput(ns("popupfig"))
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
        ##buttons = buttons,
        ##getCaption = caption.fun,
        saveHTML = saveHTML,
        outputFunc = outputFunc,
        renderFunc = renderFunc
    )
    return(res)
}

##================================================================================
##======================= TABLE MODULE ===========================================
##================================================================================


tableWidget <- function(id) {
    ns <- shiny::NS(id)
    shiny::uiOutput(ns("widget"))
}

tableModule <- function(input, output, session, 
                        func, func2=NULL, info.text="Info text",
                        title=NULL, label=NULL, server=TRUE, 
                        caption=NULL, caption2=caption,
                        csvFunc=NULL, filename="data.csv", ##inputs=NULL, 
                        ##no.download = FALSE, just.info=FALSE,
                        width=c("100%","100%"), height=c("auto","auto"),
                        options = NULL, info.width="300px"
                        )
{
    ns <- session$ns
    
    if(any(class(caption)=="reactive")) {
        caption <- caption()
    }
    if(class(caption)=="character") {
        caption <- shiny::HTML(caption)
    }
    
    options.button <- ""    
    if(!is.null(options) && length(options)>0) {
        options.button <- shinyWidgets::dropdownButton(
            options,
            ##shiny::br(),
            ##dload,
            circle = TRUE, size = "xs", ## status = "danger",
            ## icon = shiny::icon("gear"),
            icon = shiny::icon("bars"),
            width = "250px",
            inputId = ns("options"),
            tooltip = shinyWidgets::tooltipOptions(title = "Settings", placement = "right")
        )
    }

    if(!is.null(label) && label!="") label <- paste0("(",label,")")
    label1 = shiny::HTML(paste0("<span class='module-label'>",label,"</span>"))
    
    zoom.button <- modalTrigger(
        ns("zoombutton"),
        ns("tablePopup"),
        icon("window-maximize"),
        class="btn-circle-xs"
    )
    
    buttons <- shiny::fillRow(
        ##flex=c(NA,NA,NA,NA,1),
        flex=c(NA,1,NA,NA,NA,NA),
        label1,
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
        shiny::div(class='download-button',
                   shinyWidgets::dropdownButton(
                                     shiny::downloadButton(ns("csv"), "CSV"),
                                     circle = TRUE, size = "xs", ## status = "danger",
                                     icon = shiny::icon("download"), width = "80px", right=FALSE,
                                     tooltip = shinyWidgets::tooltipOptions(title = "Download",
                                                                            placement = "right")
                                 )),
        zoom.button
    )
    
    CSVFILE = paste0(gsub("file","data",tempfile()),".csv")
    CSVFILE    
    
    ## render2 <- shiny::renderPlot({plot_array[[3]]()}, res=res)
    download.csv <- shiny::downloadHandler(
        filename = filename,
        content = function(file) {
          if(!is.null(csvFunc)) {
            dt <- csvFunc()
          } else {
            dt <- func()$x$data
          }
          ##write.csv(dt, file=CSVFILE, row.names=FALSE)
          ##file.copy(CSVFILE, file, overwrite=TRUE)
          write.csv(dt, file=file, row.names=FALSE)
        }
    )
    output$csv <- download.csv
    
    if(is.null(func2)) func2 <- func
    if(length(height)==1) height <- c(height,height)
    if(length(width)==1)  width  <- c(width,width)
    ##ifnotchar.int <- function(s) ifelse(grepl("[%]$|auto|vmin|vh|vw|vmax",s),s,as.integer(s))
    ifnotchar.int <- function(s) suppressWarnings(
      ifelse(!is.na(as.integer(s)), paste0(as.integer(s),"px"), s))
    width.1  <- ifnotchar.int(width[1])
    width.2  <- ifnotchar.int(width[2])
    height.1 <- ifnotchar.int(height[1])
    height.2 <- ifnotchar.int(height[2])

    output$datatable <- DT::renderDT({
         func()
    })
    output$datatable2 <- DT::renderDT({
         func2() 
    })
    
    output$popuptable <- shiny::renderUI({
        if(any(class(caption2)=="reactive")) {
            caption2 <- caption2()
        }
        if(any(class(caption2)=="character")) {
            caption2 <- shiny::HTML(caption2)
        }
        shiny::tagList(
            shiny::div( caption2, class="caption2"),            
            DT::DTOutput(ns("datatable2"), width=width.2, height=height.2)
        )
    })
    
    output$widget <- shiny::renderUI({

        modaldialog.style <- paste0("#",ns("tablePopup")," .modal-dialog {width:",width.2,";}")
        modalbody.style   <- paste0("#",ns("tablePopup")," .modal-body {min-height:",height.2,";}")
        modalfooter.none  <- paste0("#",ns("tablePopup")," .modal-footer{display:none;}")
        div.caption <- NULL
        if(!is.null(caption)) div.caption <- shiny::div(caption, class="table-caption")

        div(class="tablewidget",
            shiny::fillCol(
              flex = c(NA,NA,1,NA),
              shiny::tags$head(shiny::tags$style(modaldialog.style)),
              shiny::tags$head(shiny::tags$style(modalbody.style)),
              shiny::tags$head(shiny::tags$style(modalfooter.none)),
              div(buttons, class="tablewidget-header"),
              div.caption,
              DT::DTOutput(ns("datatable"), width=width.1, height=height.1),
              shiny::div(class="popup-table",
                         modalUI(
                           id = ns("tablePopup"),
                           title = title,
                           size = "fullscreen",
                           shiny::uiOutput(ns("popuptable"))
                         ))            
            ))
        })
    
    module <- list(
        ##data = func,
        data = shiny::reactive(func()$x$data),
        rows_current = shiny::reactive(input$datatable_rows_current),
        rows_selected = shiny::reactive(input$datatable_rows_selected),        
        rows_all = shiny::reactive(input$datatable_rows_all),
        rownames_current = shiny::reactive({
            rns <- rownames(func()$x$data)
            if(is.null(rns)) rns <- 1:nrow(func()$x$data)
            rns[input$datatable_rows_current]
        }),
        rownames_selected = shiny::reactive({
            rns <- rownames(func()$x$data)
            if(is.null(rns)) rns <- 1:nrow(func()$x$data)
            rns[input$datatable_rows_selected]
        }),        
        rownames_all = shiny::reactive({
            rns <- rownames(func()$x$data)
            if(is.null(rns)) rns <- 1:nrow(func()$x$data)
            rns[input$datatable_rows_all]
        })
    )
    return(module)
}
