########################################################################
## PLOT/TABLE MODULES
########################################################################

ADDSIGNATURE=TRUE
ADDSIGNATURE=FALSE

## ----------------- create widget
moduleWidget <- function(module, outputFunc="plotOutput", ns=NULL,
                         height="100%", width="100%")
{
    ##module.id = ns(module$id)
    module.id = module$id
    if(!is.null(ns)) module.id <- ns(module.id)    
    if(!is.null(module$outputFunc)) outputFunc = module$outputFunc
    if(!is.character(outputFunc)) outputFunc = as.character(quote(outputFunc))
    outputFunc2 <- paste0(outputFunc,"('",module.id,"', height='",
                          height,"', width='",width,"')")
    p <- fillCol(
        flex = c(NA,1,NA),
        module$button,
        eval(parse(text=outputFunc2)),
        div(HTML(module$getCaption()),class="caption")
    )
    ## p <- box(p)
    p
}

attachModule <- function(output, module, ns=NULL) {
    ## Attach render functions to shiny output. Works for both
    ## tableModule and plotModule.
    ##
    ##
    if(!"id" %in% names(module)) stop("module must have id")
    id <- module$id
    ##if(!is.null(ns)) id <- ns(module$id) ## NEED RETHINK!!!
    output[[id]] <- module$render
    if("csv" %in% names(module)) output[[paste0(id,"_csv")]] <- module$csv
    if("pdf" %in% names(module)) output[[paste0(id,"_pdf")]] <- module$pdf
    if("png" %in% names(module)) output[[paste0(id,"_png")]] <- module$png
    if("html" %in% names(module)) output[[paste0(id,"_html")]] <- module$html
    return(output)
}
##attachModule(enrich_fctable_module, "enrich_fctable")

plotModuleButtons <- function(id, text="Help text", title="", ns=NULL,
                              download.fmt = c("png","pdf"), info.width="300px",
                              no.download = FALSE, just.info=FALSE, 
                              options=NULL, inputs=NULL, label="")
{
    require(shinyWidgets)
    require(shinyBS)

    ## all widget IDs
    if(is.null(ns)) ns=function(x){x}
    info_id    <- ns(paste0(id,"_info"))
    options_id <- ns(paste0(id,"_options"))
    html_id <- ns(paste0(id,"_html"))
    pdf_id <- ns(paste0(id,"_pdf"))
    png_id <- ns(paste0(id,"_png"))
    csv_id <- ns(paste0(id,"_csv"))

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
            inputId = options_id,
            tooltip = tooltipOptions(title = "Settings", placement = "right")
        )
    }

    dload.csv = dload.pdf = dload.png = dload.html = NULL
    if("pdf" %in% download.fmt)  dload.pdf  <- downloadButton(pdf_id, "PDF")
    if("png" %in% download.fmt)  dload.png  <- downloadButton(png_id, "PNG")
    if("csv" %in% download.fmt)  dload.csv  <- downloadButton(csv_id, "CSV")
    if("html" %in% download.fmt) dload.html <- downloadButton(html_id, "HTML")
    
    dload.button <- dropdownButton(
        dload.pdf,
        dload.png,
        dload.csv,
        dload.html,
        circle = TRUE, size = "xs", ## status = "danger",
        icon = icon("download"), width = "40px",
        tooltip = tooltipOptions(title = "Download", placement = "right")
    )
    if(no.download || length(download.fmt)==0 ) dload.button <- ""
    label1 = HTML(paste0("<span class='module-label'>",label,"</span>"))
    
    ## button layout
    btn <- fillRow(
        flex=c(NA,NA,NA,NA,1),
        label1,
        dropdownButton(
            tags$p(HTML(text)),
            br(),
            circle = TRUE, size = "xs", ## status = "danger",
            icon = icon("info"), width = info.width,
            inputId = info_id,
            tooltip = tooltipOptions(title = "Info", placement = "right")
        ),
        options.button,
        div(class='download-button', dload.button),
        ##HTML(paste("<center><strong>",title,"</strong></center>"))
        HTML(paste("<center>",title,"</center>"))
        ##br()
        ##inputs
        ##selectInput("sel123","number",1:10)
    )
    return(btn)
}

plotModule <- function(id, func, info.text="Info text", title="", ns=NULL,
                       inputs=NULL, options = NULL, label="", caption="", 
                       plotlib = "base", renderFunc=NULL, outputFunc=NULL,
                       no.download = FALSE, download.fmt=c("png","pdf","html"), 
                       just.info=FALSE, server=TRUE, info.width="300px",
                       width = "auto", height = "auto",  ## for renderPlot
                       pdf.width=8, pdf.height=8, pdf.pointsize=12, res=72)
{
    require(shinyWidgets)
    require(shinyBS)
    require(webshot)

    ## no NS here!!
    ## info_id <- paste0(id,"_info")
    ## options_id <- paste0(id,"_options")
    ## pdf_id <- paste0(id,"_pdf")
    ## png_id <- paste0(id,"_png")
    ## html_id <- ns(paste0(id,"_html"))
    
    ## these engines cannot (yet) provide html
    if(plotlib %in% c("base")) {    
        download.fmt <- setdiff(download.fmt, c("html"))
    }
    
    buttons = plotModuleButtons(
        id=id, text=info.text, title=title,
        ns=ns, ## pass through any namespace prefix
        just.info=just.info, label = label,
        inputs=inputs, options=options,
        info.width=info.width,
        no.download=no.download,
        download.fmt=download.fmt)
    ##toggleDropdownButton(info_id)
    
    
    if(any(class(caption)=="reactive")) {
        caption.fun <- caption
    } else {
        caption.fun <- function() { caption }
    }

    do.pdf = "pdf" %in% download.fmt
    do.png = "png" %in% download.fmt
    do.html = "html" %in% download.fmt

    PNGFILE=PDFFILE=HTMLFILE=NULL
    if(do.pdf) PDFFILE = paste0(gsub("file","plot",tempfile()),".pdf")
    if(do.png) PNGFILE = paste0(gsub("file","plot",tempfile()),".png")
    HTMLFILE = paste0(gsub("file","plot",tempfile()),".html")  ## tempory for webshot
    ##HTMLFILE = "/tmp/plot.html"
    HTMLFILE
    unlink(HTMLFILE)
    
    if(plotlib == "generic") {
        if(is.null(renderFunc)) error("'generic' class must provide renderFunc")
        renderFUN <- eval(parse(text=renderFunc))
        render <- renderFUN({
            func()
        })
    } else if(plotlib == "plotly") {
        render <- renderPlotly({
            p <- func()            
            p
        })
        outputFunc = "plotly::plotlyOutput"
    } else if(plotlib %in% c("htmlwidget","pairsD3","scatterD3")) {
        require(htmlwidgets)
        if(plotlib=="htmlwidget" && is.null(renderFunc)) {
            error("'htmlwidget' class must provide renderFunc")
        }
        if(plotlib=="scatterD3") {
            require(scatterD3)
            renderFunc="scatterD3::renderScatterD3"
            outputFunc="scatterD3::scatterD3Output"
        }
        if(plotlib=="pairsD3") {
            require(pairsD3)
            renderFunc="pairsD3::renderPairsD3"
            outputFunc="pairsD3::pairsD3Output"
        }
        renderFUN <- eval(parse(text=renderFunc))
        render <- renderFUN({
            p <- func()
            return(p)
        })
    } else if(plotlib == "visnetwork") {
        require(visNetwork)
        render <- renderVisNetwork({
            p <- func()
            return(p)
        })
        outputFunc="visNetwork::visNetworkOutput"
    } else if(plotlib %in% c("ggplot","ggplot2")) {
        render <- renderPlot({
            plot(func())
        }, res=res)
        outputFunc="plotOutput"
    } else if(plotlib == "iheatmapr") {
        require(iheatmapr)
        render <- renderIheatmap({
            p <- func()            
            p
        })
        outputFunc="iheatmapr::iheatmaprOutput"
    } else {
        ##------------------------------------------------------------
        ## Base plotting
        ##------------------------------------------------------------
        render <- renderPlot({
            par(mar=c(0,0,0,0),oma=c(0,0,0,0))
            func()

            ## NEEEDS FIX!! pdf generating should be done just in
            ## renderPlot, not here. But cannot get it to work (IK
            ## 19.10.02)
            suppressWarnings( suppressMessages(
                if(do.pdf) {
                    if(ADDSIGNATURE) {
                        mtext("created with Omics Playground",
                              1,line=-1,outer=TRUE,adj=0.98,padj=0,cex=0.6,col="#44444444")
                    }
                    dev.print(pdf, file=PDFFILE, width=pdf.width, height=pdf.height,
                              pointsize=pdf.pointsize)
                    ##dev.copy2pdf(file=PDFFILE, width=pdf.width, height=pdf.height)
                    ##dev.off()  ## important!!
                }
            ))
            
            suppressWarnings( suppressMessages(
                if(do.png) {
                    if(ADDSIGNATURE) {
                        mtext("created with Omics Playground",
                              1,line=-1,outer=TRUE,adj=0.98,padj=0,cex=0.6,col="#44444444")
                    }
                    dev.print(png, file=PNGFILE, width=pdf.width*100, height=pdf.height*100,
                              pointsize=pdf.pointsize)
                    ##dev.copy2pdf(file=PDFFILE, width=pdf.width, height=pdf.height)
                    ##dev.off()  ## important!!
                }
            ))

        }, res=res, width=width, height=height )
        outputFunc="plotOutput"
    }

    ## render2 <- renderPlot({plot_array[[3]]()}, res=res)
    download.pdf = download.png = download.html = NULL

    ##============================================================
    ##=============== Download Handlers ==========================
    ##============================================================
    addSignaturePlotly <- function(p) {
        add_annotations(
            p,
            x = 1, y=-0.05,
            ##yshift = -100, 
            xref="paper", yref="paper",
            text = "created with Omics Playground",
            xanchor = "right", showarrow=FALSE,
            font = list(size=6, color="#44444444")
        )
    }

    if(do.png) {
        download.png <- downloadHandler(
            filename = "plot.png",
            content = function(file) {
                withProgress({
                    ## unlink(PNGFILE) ## do not remove!
                    if(plotlib=="plotly") {
                        cat("downloadHandler:: exporting plotly to PNG\n")
                        p <- func()
                        if(ADDSIGNATURE) p = addSignaturePlotly(p)                             
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
                        if(ADDSIGNATURE) p = addSignaturePlotly(p)                             
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

    if(do.pdf) {        
        download.pdf <- downloadHandler(
            filename = "plot.pdf",
            content = function(file) {
                withProgress({
                    ## unlink(PDFFILE) ## do not remove!
                    if(plotlib=="plotly") {
                        cat("downloadHandler:: exporting plotly to PDF\n")
                        p <- func()
                        if(ADDSIGNATURE) p = addSignaturePlotly(p)                             
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
                        if(ADDSIGNATURE) p = addSignaturePlotly(p)                             
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
            if(ADDSIGNATURE) p <- addSignaturePlotly(p)
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
    
    if(do.html)  {
        download.html <- downloadHandler(
            filename = "plot.html",
            content = function(file) {
                withProgress({
                    ## unlink(HTMLFILE) ## do not remove!
                    if(plotlib == "plotly" ) {
                        cat("downloadHandler:: exporting plotly to HTML\n")
                        p <- func()
                        if(ADDSIGNATURE) p <- addSignaturePlotly(p)
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
    
    
    ## ------------- create a widget????
    ##widget <- function(w,h) box(fillCol(flex=c(NA,1),button,render),width=w,height=h)
    
    module <- list(
        id = id,
        .func = func,
        .tmpfiles = c(pdf=PDFFILE, html=HTMLFILE),
        render = render,
        ##render2 = render2,
        pdf = download.pdf,
        png = download.png,
        html = download.html,
        buttons = buttons,
        getCaption = caption.fun,
        saveHTML = saveHTML,
        outputFunc = outputFunc
    )
    ## attr(module, "class") <- "ShinyModule"
    return(module)
}


tableModule <- function(id, func, info.text="Info text",
                        title="", label="", caption="", ns=NULL,
                        ##inputs=NULL, 
                        options = NULL, info.width="300px",
                        ##no.download = FALSE, just.info=FALSE,
                        server=TRUE)
{
    require(shinyWidgets)
    require(shinyBS)

    ## all callback ids
    if(is.null(ns)) ns=function(x)x
    info_id <- paste0(id,"_info") ## NS?
    options_id <- paste0(id,"_options")  ## NS?
    pdf_id <- ns(paste0(id,"_pdf"))
    csv_id <- ns(paste0(id,"_csv"))    
    label1 = HTML(paste0("<span class='module-label'>",label,"</span>"))

    if(any(class(caption)=="reactive")) {
        caption.fun <- caption
    } else {
        caption.fun <- function() { caption }
    }
    
    options.button <- ""    
    if(!is.null(options) && length(options)>0) {
        options.button <- dropdownButton(
            ##tags$h3("Options"),
            options,
            ##br(),
            ##dload,
            circle = TRUE, size = "xs", ## status = "danger",
            ## icon = icon("gear"),
            icon = icon("bars"),
            width = "250px",
            inputId = options_id,
            tooltip = tooltipOptions(title = "Settings", placement = "right")
        )
    }
    
    buttons <- fillRow(
        flex=c(NA,NA,NA,NA,1),
        ##actionLink(options_id, label=NULL, icon = icon("info")),
        label1,
        dropdownButton(
            tags$p(HTML(info.text)),
            br(),
            circle = TRUE, size = "xs", ## status = "danger",
            icon = icon("info"), width = info.width,
            inputId = info_id,
            tooltip = tooltipOptions(title = "Info", placement = "right")
        ),
        options.button,
        div(class='download-button', dropdownButton(
            downloadButton(csv_id, "CSV"),
            circle = TRUE, size = "xs", ## status = "danger",
            icon = icon("download"), width = "80px",
            tooltip = tooltipOptions(title = "Download", placement = "right")
        )),
        ##HTML(paste("<center><strong>",title,"</strong></center>"))
        HTML(paste("<center>",title,"</center>"))
        ##inputs
        ##selectInput("sel123","number",1:10)
    )
    
    CSVFILE = paste0(gsub("file","data",tempfile()),".csv")
    CSVFILE
    
    render <- renderDataTable({
        dt <- func()
        write.csv(dt$x$data, file=CSVFILE, row.names=FALSE)
        dt
    }, server=server)
    
    ## render2 <- renderPlot({plot_array[[3]]()}, res=res)
    download.csv <- downloadHandler(
        filename = "data.csv",
        content = function(file) file.copy(CSVFILE,file)        
    )
    ##box <- fillCol(flex=c(NA,1), button, render)
    
    module <- list(
        id = id,
        .func = func,
        .tmpfiles = c(csv=CSVFILE),
        getCaption = caption.fun,
        render = render,
        ##render2 = render2,
        csv = download.csv,
        buttons = buttons
    )
    ## attr(module, "class") <- "ShinyModule"
    return(module)
}


pgx.randomSlogan <- function(b=NULL) {
    motto.list <- c("Skip the queue. Take the fast lane",
                    "Analyze your omics data 10x faster, 100x smarter",
                    "Self-service analytics for Big Omics data",
                    "Do-it-yourself. Yes you can",
                    "Fasten your seat belts. Hi-speed data analytics",
                    "Zero coding required",
                    "Make Data Analytics Great Again",
                    "Smart tools for biodata analytics",
                    "No need to write pesky R scripts",
                    "Analytics anywhere. Anytime",
                    "Analyze with confidence. Be a rockstar",
                    "Our platform. Your solution",
                    "Integrate more. Dig deeper",
                    "Visual analytics. See and understand",
                    "Explore omics data freely",
                    "Where do you want to go today?",
                    "Analytics the way you want",
                    "Your analysis doesn't go on coffee breaks",
                    "Fast track your bioinformatics needs"
                    )
    s <- sample(motto.list,1)
    if(!is.null(b)) s <- breakstring2(s,n=b)
    s
}

########################################################################
########################################################################
########################################################################
