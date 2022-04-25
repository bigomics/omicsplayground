##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

IntersectionBoard <- function(id, inputData, selected_gxmethods, selected_gsetmethods)
{
  moduleServer(id, function(input, output, session)
  {
    ns <- session$ns ## NAMESPACE
    fullH = 800       # row height of panel 

    infotext =
        "The <strong>Intersection analysis module</strong> enables users to compare multiple contrasts by intersecting the genes of profiles. The main goal is to identify contrasts showing similar profiles.

<br><br>For the selected contrasts, the platform provides volcano plots and pairwise correlation plots between the profiles in the <strong>Pairs</strong> panel. Simultaneously, a Venn diagram with the number of intersecting genes between the profiles is plotted in <strong>Venn diagram</strong> panel. Details of intersecting genes are also reported in an interactive table. A more detailed scatter plot of two profiles is possible under the <strong>Two-pairs</strong> panel. Users can check the pairwise correlations of the contrasts under the <b>Contrast heatmap</b> panel. Alternatively, the <strong>Connectivity Map (CMap)</strong> shows the similarity of the contrasts profiles as a t-SNE plot.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=5' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>

"

    ## delayed input
    input_comparisons <- shiny::reactive({
        input$comparisons
    }) %>% shiny::debounce(500)    
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    shiny::observeEvent( input$info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Intersection Analysis Board</strong>"),
            shiny::HTML(infotext),
            easyClose = TRUE, size="l" ))
    })

    ## update choices upon change of data set 
    shiny::observe({
        ngs <- inputData()
        ##req(ngs)
        if(is.null(ngs)) return(NULL)
        comparisons <- colnames(ngs$model.parameters$contr.matrix)
        comparisons <- sort(comparisons)
        shiny::updateSelectInput(session, "comparisons", choices=comparisons,
                          selected=head(comparisons,3))
    })

    ## update choices upon change of feature level
    ##observeEvent( input$level, {
    shiny::observe({
        ngs <- inputData()
        ## shiny::req(ngs,input$level)
        if(is.null(ngs)) return(NULL)
        shiny::req(input$level)
        ##flt.choices = names(ngs$families)
        if(input$level=="geneset") {
            ft <- names(COLLECTIONS)
            nn <- sapply(COLLECTIONS, function(x) sum(x %in% rownames(ngs$gsetX)))
            ft <- ft[nn >= 10]
        } else {
            ## gene level
            ft <- pgx.getFamilies(ngs,nmin=10,extended=FALSE)
        }
        ft <- sort(ft)
        ## if(input$level=="gene") ft = sort(c("<custom>",ft))
        ## ft = sort(c("<custom>",ft))        
        shiny::updateSelectInput(session, "filter", choices=ft, selected="<all>")
    })
    
    ## shiny::observe({
    ##     splom.sel <- plotly::event_data("plotly_selected", source="splom")
    ##     sel.keys <- as.character(splom.sel$key)
    ##     if(0 && length(sel.keys)>0) {
    ##         shiny::updateSelectInput(session, "filter", selected="<custom>")
    ##         sel.keys = paste(sel.keys, collapse=" ")
    ##         shiny::updateTextAreaInput(session, "customlist", value=sel.keys)
    ##     }
    ## })
    
    
    ##================================================================================
    ##========================= REACTIVE FUNCTIONS ===================================
    ##================================================================================
    
    getFoldChangeMatrix <- shiny::reactive({
        ## 
        ## Get full foldchange matrix from ngs object.
        ##
        ##
        ##
        ##dbg("<intersectionBoard:getFoldChangeMatrix> reacted\n")
        fc0 = NULL
        qv0 = NULL
        ngs <- inputData()
        alertDataLoaded(session,ngs)        
        shiny::req(ngs)
        
        sel = names(ngs$gset.meta$meta)
        ##sel = input_comparisons()
        ##sel = intersect(sel, names(ngs$gset.meta$meta))
        ##if(length(sel)==0) return(NULL)
        
        if(input$level=="geneset") {
            gsetmethods <- c("gsva","camera","fgsea")
            gsetmethods <- selected_gsetmethods()
            if(length(gsetmethods)<1 || gsetmethods[1]=="") return(NULL)
            
            ##fc0 = sapply(ngs$gset.meta$meta[sel], function(x)
            ##    rowMeans(unclass(x$fc)[,gsetmethods,drop=FALSE]))
            fc0 = sapply(ngs$gset.meta$meta[sel], function(x) x$meta.fx)
            rownames(fc0) <- rownames(ngs$gset.meta$meta[[1]])
            qv0 = sapply(ngs$gset.meta$meta[sel], function(x)
                apply(unclass(x$q)[,gsetmethods,drop=FALSE],1,max))
            rownames(qv0) <- rownames(ngs$gset.meta$meta[[1]])

            ## apply user selected filter
            gsets = rownames(fc0)
            if(input$filter == "<custom>") {
                gsets = strsplit( input$customlist, split="[, ;]")[[1]]
                if(length(gsets)>0) {
                    gsets = intersect(rownames(ngs$gsetX), gsets)
                }
            } else if(input$filter != "<all>") {
                gsets = unique(unlist(COLLECTIONS[input$filter]))
            }
            gsets = intersect(gsets, rownames(fc0))
            fc1 = fc0[gsets,,drop=FALSE]
            qv1 = qv0[gsets,,drop=FALSE]
        } else {
            ## Gene
            ##
            gxmethods <- "trend.limma"
            gxmethods <- c("trend.limma","edger.qlf","deseq2.wald")
            gxmethods <- selected_gxmethods()  ## reactive object from EXPRESSION section

            dbg("[IntersectionBoard:getFoldChangeMatrix] gxmethods = ",gxmethods)

            mq1 <- ngs$gx.meta$meta[[1]]$meta.q
            dbg("[IntersectionBoard:getFoldChangeMatrix] head.meta.q = ",head(mq1))
            
            if(length(gxmethods)<1 || gxmethods[1]=="") return(NULL)
            
            fc0 = sapply(ngs$gx.meta$meta[sel], function(x) x$meta.fx)
            rownames(fc0) <- rownames(ngs$gx.meta$meta[[1]])
            qv0 = sapply(ngs$gx.meta$meta[sel], function(x)
                apply(unclass(x$q)[,gxmethods,drop=FALSE],1,max))
            rownames(qv0) <- rownames(ngs$gx.meta$meta[[1]])
            dim(fc0)
            dim(qv0)

            ## filter with active filter
            sel.probes = rownames(fc0) ## default to all probes
            if(input$filter == "<custom>") {
                genes = strsplit( input$customlist, split="[, ;]")[[1]]
                if(length(genes)>0) {
                    sel.probes = filterProbes(ngs$genes, genes)
                }
            } else  if(input$filter != "<all>") {
                ##gset <- GSETS[[input$filter]]
                gset.genes <- unlist(getGSETS(input$filter))
                sel.probes = filterProbes(ngs$genes, gset.genes)
            }
            sel.probes = intersect(sel.probes, rownames(fc0))
            fc1 = fc0[sel.probes,,drop=FALSE]
            qv1 = qv0[sel.probes,,drop=FALSE]
        }    
        fc1 <- fc1[,!duplicated(colnames(fc1)),drop=FALSE]
        qv1 <- qv1[,!duplicated(colnames(qv1)),drop=FALSE]
        
        res = list(fc=fc1, qv=qv1, fc.full=fc0, qv.full=qv0)
        return(res) 
    })

    
    getActiveFoldChangeMatrix <- shiny::reactive({
        
        res = getFoldChangeMatrix()
        ##if(is.null(res)) return(NULL)
        shiny::req(res)            

        ## match with selected/active contrasts
        ## comp = head(colnames(res$fc),3)
        comp = input_comparisons()
        kk = match(comp, colnames(res$fc))
        if(length(kk)==0) return(NULL)
        if(length(kk)==1) kk = c(kk,kk)
        res$fc = res$fc[,kk,drop=FALSE]
        res$qv = res$qv[,kk,drop=FALSE]        
        res$fc.full = res$fc.full[,kk,drop=FALSE]
        res$qv.full = res$qv.full[,kk,drop=FALSE]        

        return(res)
    })
    
    
    getSignificanceCalls <- shiny::reactive({
        ## Gets the matrix of significance calls.
        ##
        ngs <- inputData()

        sel = head(names(ngs$gset.meta$meta),7)
        sel = input_comparisons()
        sel = intersect(sel, names(ngs$gset.meta$meta))
        if(length(sel)==0) return(NULL)        

        res <- getFoldChangeMatrix()    
        fc <- res$fc[,sel,drop=FALSE]
        qv <- res$qv[,sel,drop=FALSE]

        fdr=0.05;lfc=0.2
        fdr = as.numeric(input$fdr)
        lfc = as.numeric(input$lfc)
        dt = sign(fc) * (qv <= fdr & abs(fc) >= lfc)
        dt[is.na(dt)] = 0
        ## add label of venn intersection region
        dt.labels = LETTERS[1:ncol(dt)]

        venn.intersection = apply( 1*(dt!=0), 1, function(x)
            paste(dt.labels[which(x==1)],collapse=""))
        dt = data.frame( intersection=venn.intersection, dt, check.names=FALSE )

        ## update filter
        choices <- c("<all>",sort(unique(venn.intersection)))
        selected <- venn.intersection[which.max(nchar(venn.intersection))]
        dbg("[getSignificanceCalls] choiced = ",choices)
        shiny::updateSelectInput(
            session, "venntable_intersection",
            choices = choices,
            ##selected = "<all>"
            selected = selected
        )
        
        return(dt)
    })


    getSignificantFoldChangeMatrix <- shiny::reactive({
        ##
        ## Filters FC matrix with significance and user-defined
        ## intersection region.
        dt <- getSignificanceCalls()
        shiny::req(dt)
        
        isect = input$intersection
        fc0 = getFoldChangeMatrix()$fc  
        if( length(isect) == 0) {
            fc1 = fc0
        } else {
            ## only genes at least significant in one group
            jj = which(rowSums(dt[,2:ncol(dt),drop=FALSE]!=0)>0)
            if(length(jj)==0) return(NULL)
            dt = dt[jj,,drop=FALSE]

            ## check same sign
            if(input$include=="up/down") {
                kk = 1 + match(c("B","C"),LETTERS[1:10])
                kk = 1 + match(isect,LETTERS[1:10])
                kk <- intersect(kk, 1:ncol(dt))
                
                dt1 = dt[,kk,drop=FALSE]    
                jj = which( rowMeans(sign(dt1)== +1)==1 |
                            (rowMeans(sign(dt1)== -1)==1) )    
                dt = dt[jj,,drop=FALSE]    
                remove(dt1)
            }
            
            ## only genes in the selected intersection
            intersection="ABC"
            intersection = paste0(input$intersection,collapse="")
            dt = dt[which(dt$intersection == intersection),,drop=FALSE]

        }
        
        ## filtered by family/collection
        fc1 = fc0[intersect(rownames(dt),rownames(fc0)),,drop=FALSE]
        if(nrow(dt)==1) {
            fc1 = matrix(fc1,nrow=1)
            rownames(fc1) <- rownames(dt)
            colnames(fc1) <- colnames(fc0)
        }

        ## filtered by SPLOM selection
        splom.sel <- plotly::event_data("plotly_selected", source="splom")
        sel.keys <- as.character(splom.sel$key)
        if(1 && length(sel.keys)>0) {
            sel <- intersect(sel.keys, rownames(fc1))
            fc1 = fc1[sel,,drop=FALSE]
        }
        
        ## only active/selected comparisons
        sel = colnames(dt)[-1]
        kk = match(sel,gsub(" \\(-\\)","",colnames(fc1)))
        fc1 = fc1[,kk,drop=FALSE]

        ## order
        fc1 <- fc1[order(-rowMeans(fc1)),,drop=FALSE]        
        fc1 <- round(fc1, digits=3)
        colnames(fc1) = LETTERS[1:ncol(fc1)]
        ##fc0 = data.frame(fc0)
        
        ## add intersection code
        sel <- match(rownames(fc1),rownames(dt))
        fc1 <- data.frame(intersection=dt$intersection[sel], fc=fc1)
        
        ## filter on user selection
        vv <- input$venntable_intersection
        if(vv != "<all>") {
            sel <- which(fc1$intersection == vv)
            fc1 <- fc1[sel,,drop=FALSE]
        }
        return(fc1)
    })

    
    getCurrentSig <- shiny::reactive({
        ## Switch between FC profile or NMF vectors
        ##
        ##
        ngs <- inputData()
        shiny::req(ngs)
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        
        ##------------ UMAP clustering (genes) -----------------
        dbg("[getCurrentSig] calculating UMAP for genes...")
        progress$inc(0.33,"calculating UMAP for genes...")
        if("cluster.genes" %in% names(ngs)) {
            pos <- ngs$cluster.genes$pos[['umap2d']]
        } else {
            X1 <- ngs$X
            X1 <- (X1 - rowMeans(X1)) / mean(apply(X1,1,sd,na.rm=TRUE))
            pos <- pgx.clusterBigMatrix(
                t(X1), methods='umap', dims=2, reduce.sd=-1)[[1]]
            pos <- pos.compact(pos)
        }

        ##------------ UMAP clustering (genesets) -----------------
        dbg("[getCurrentSig] calculating UMAP for genesets...")
        progress$inc(0.33,"calculating UMAP for genesets...")
        if("cluster.gsets" %in% names(ngs)) {
            gsea.pos <- ngs$cluster.gsets$pos[['umap2d']]
        } else {
            X2 <- ngs$gsetX
            X2 <- (X2 - rowMeans(X2)) / mean(apply(X2,1,sd,na.rm=TRUE))
            gsea.pos <- pgx.clusterBigMatrix(
                t(X2), methods='umap', dims=2, reduce.sd=-1)[[1]]
            gsea.pos <- pos.compact(gsea.pos)
            dim(gsea.pos)
        }

        ##------------ get signature matrices -----------------
        F <- pgx.getMetaMatrix(ngs, level='gene')
        G <- pgx.getMetaMatrix(ngs, level='geneset')
        ##f.score <- F$fc * -log10(F$qv)
        ##g.score <- G$fc * -log10(G$qv)
        f.score <- F$fc * (1 - F$qv)**4  ## q-weighted FC
        g.score <- G$fc * (1 - G$qv)**4
        
        ii <- intersect(rownames(pos),rownames(f.score))            
        sig  <- f.score[ii,,drop=FALSE]
        pos  <- pos[ii,]                        
        ii <- order(-rowMeans(sig))
        sig <- sig[ii,,drop=FALSE]
        pos <- pos[ii,]
        
        ii <- intersect(rownames(gsea.pos),rownames(g.score))            
        gsea  <- g.score[ii,,drop=FALSE]
        gsea.pos  <- gsea.pos[ii,]                        
        ii <- order(-rowMeans(gsea))
        gsea <- gsea[ii,,drop=FALSE]
        gsea.pos <- gsea.pos[ii,]
        
        out <- list(sig=sig, pos=pos, gsea=gsea, gsea.pos=gsea.pos)

        dbg("[getCurrentSig] done")                        
        progress$close()
        out
    })

    
    ##======================================================================
    ## PLOTTING MODULES
    ##======================================================================
    
    ##======================================================================
    ## Scatterplot matrix in plotly
    ##
    ## From: https://plot.ly/r/splom/
    ##======================================================================
    
    scatterPlotMatrix.PLOT <- shiny::reactive({

        dbg("[IntersectionBoard::scatterPlotMatrix.PLOT]  reacted\n")





        ##res = pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        res = getActiveFoldChangeMatrix()
        fc0 = res$fc.full
        fc1 = res$fc
        
        dim(fc0)
        if(is.null(res)) return(NULL)
        
        ##fc0 <- fc0[order(-rowSums(fc0**2)),]
        fc0 <- fc0[order(-apply(abs(fc0),1,max)),]
        fc0 <- fc0[order(-rowMeans(abs(fc0**2))),]
        
        ## selected genes
        sel.genes = grep("^CD",rownames(fc0),value=TRUE)
        sel.genes = head(rownames(fc0),100)  ## top100
        sel.genes = intersect(rownames(fc0),rownames(fc1))
        sel.genes = intersect(sel.genes, rownames(fc0))
        head(sel.genes)

        ## subsample for speed: take top1000 + 1000
        df <- data.frame(fc0)
        if(1) {
            ntop = 99999
            ##ntop <- input$splom_ntop        
            jj <- match(sel.genes, rownames(df))
            jj <- c(jj, 1:min(ntop,nrow(df)))
            if(nrow(df)>ntop) {
                nremain = setdiff(1:nrow(df),jj)
                jj = c(jj, sample(nremain,min(1000,length(nremain))))  ## add 1000 random
            }
            jj <- unique(jj)
            ##df <- data.frame(head(fc0,ntop))
            df <- data.frame(df[jj,])
        }
        dim(df)
        
        ## resort selection so that selected genes are drawn last to avoid
        ## covering them up.
        is.sel = (rownames(df) %in% sel.genes)
        if(input$splom_highlight) {
            df.color = c("#00000033","#0066FF")[1 + is.sel]
            df.color = c("#AAAAAA","#1e60BB")[1 + is.sel]
            df.color = c("#CCCCCC22","#1e60BB88")[1 + is.sel]
        } else {
            df.color = rep("#00000088",nrow(df))
            df.color = rep("#1e60BB88",nrow(df))
        }
        
        ## Labels for top 50 
        label.text = label.text0 = head(rownames(df)[which(is.sel)],50)
        label.text <- sub(".*[:]","",label.text)  ## strip prefix??
        label.text <- shortstring(label.text,30)
        if(sum(is.na(label.text))) label.text[is.na(label.text)] <- ""
        
        ## reorder so the selected genes don't get overlapped
        jj <- order(is.sel)
        df <- df[jj,]
        df.color <- df.color[jj]
        sel1 <- match(label.text0, rownames(df))  ## index for labeled

        ## Tooltip text for all 
        tt <- rownames(df)  ## strip prefix
        ## tt <- sub("","",tt)  ## strip prefix??
        if(input$level == "gene") {
            ngs <- inputData()
            g <- rownames(df)
            tt <- paste0("<b>",g,"</b> ", ngs$genes[g,"gene_title"])
        }
        tt <- gsub("_"," ",tt)
        tt <- sapply(tt,breakstring2,50,brk="<br>")

        ## plotly
        ##   
        axis = list( showline=TRUE,
                    zeroline=TRUE,
                    gridcolor='#dddf',
                    ticklen=4
                    )
        
        if(ncol(df)<=2) {

            ## ----------------------------------------------------
            ## Single pairs plot
            ## ----------------------------------------------------
            
            rho = cor(df[,1], df[,2])        
            annot.rho <- list(
                text = paste("r=",round(rho,4)),
                font = list(size=14),
                align = 'left',
                showarrow = FALSE,
                xref = 'paper',
                yref = 'paper',
                x = 0.03,
                y = 0.97,
                borderpad = 8, 
                bordercolor = 'black',
                borderwidth = 0.6)

            p <- plotly::plot_ly( data=df[,1:2], x = df[,1], y= df[,2],
                         type = 'scattergl', mode = 'markers',
                         source='splom', key = rownames(df),
                         ##type = 'scatter', mode = 'markers',
                         text = tt,
                         hovertemplate = paste0("<br>%{text}<br>x: %{x}<br>y: %{y}<extra></extra>"),
                         marker = list(
                             color = df.color,
                             size = 8,
                             line = list(
                                 width = 0.3,
                                 ##color = 'rgb(230,230,230)'
                                 color = 'rgb(0,0,0)'
                             ))
                         ) %>%
                plotly::add_annotations(
                    x = df[sel1,1],
                    y = df[sel1,2],
                    text = as.character(label.text),
                    ##text = rep("mylabel",length(sel1)),
                    ##xanchor = 'left',
                    xanchor = 'center',
                    yanchor = 'top',
                    font = list(size=14),
                    xref = "x",
                    yref = "y",
                    showarrow = FALSE,
                    ax = 20,
                    ay = -40
                ) %>%
                plotly::layout(
                    ## title= 'Scatterplot',
                    annotations = annot.rho,
                    hovermode = 'closest',
                    dragmode= 'select',
                    ##plot_bgcolor='rgba(240,240,240, 0.95)',
                    ## template = "plotly_dark",
                    xaxis = c(title = paste(colnames(df)[1],"   (logFC)"), axis),
                    yaxis = c(title = paste(colnames(df)[2],"   (logFC)"), axis)
                ) 
            
        } else  {

            ## ----------------------------------------------------
            ## Scatter pairs matrix plot
            ## ----------------------------------------------------

            dimensions = lapply(colnames(df), function(a) list(label=a, values = df[,a]))

            ## compute correlations 
            rho = cor(df)
            rho.text = paste("r=",as.vector(round(rho,digits=3)))
            n = ncol(df)

            ## annotation positions (approximated by eye...)
            xann = 1.02*(as.vector(mapply(rep,seq(0,0.98,1/n),n)) + 0.05*1/n)
            ##xann = as.vector(mapply(rep,seq(0,1,1/(n-1)),n))
            yann = 1.08*(as.vector(rep(seq(1,0.02,-1/n),n)) - 0.15*1/n - 0.04)
            ##yann = as.vector(rep(seq(1,0.0,-1/(n-1)),n))
            
            p <- plotly::plot_ly(df, source="splom", key=rownames(df) ) %>%
                plotly::add_trace(
                    type = 'splom',
                    dimensions = dimensions,
                    text = tt,
                    hovertemplate = paste0("<br>%{text}<br>x: %{x}<br>y: %{y}<extra></extra>"),
                    marker = list(
                        color = df.color,
                        ## colorscale = pl_colorscale,
                        size = 5,
                        line = list(
                            width = 0.3,
                            ##color = 'rgb(230,230,230)'
                            color = 'rgb(0,0,0)'
                        )
                    )
                ) %>%
                plotly::add_annotations(
                    x = xann,
                    y = yann,
                    text = rho.text,
                    font = list(size=11),
                    xanchor = "left",
                    align = 'left',
                    showarrow = FALSE,
                    xref = 'paper',
                    yref = 'paper',
                    borderpad = 3, 
                    bordercolor = 'black',
                    borderwidth = 0.6) %>%
                plotly::layout(
                    ##title= 'Scatterplot matrix',
                    hovermode='closest',
                    dragmode= 'select',
                    ##annotations = annot,
                    ## plot_bgcolor='rgba(240,240,240, 0.95)',
                    ## template = "plotly_dark",
                    xaxis = c( domain=NULL, axis),
                    yaxis = c( domain=NULL, axis),
                    xaxis2=axis, xaxis3=axis, xaxis4=axis, xaxis5=axis, xaxis6=axis, xaxis7=axis,
                    yaxis2=axis, yaxis3=axis, yaxis4=axis, yaxis5=axis, yaxis6=axis, yaxis7=axis
                )
            ## %>% plotly::style(diagonal = list(visible = F))

        }
        
        p <- p %>%
            plotly::layout(margin = list(80,80,80,80) )  ## l,r,b,t

        p <- p %>%
            ##config(displayModeBar = FALSE) %>% ## disable buttons
            plotly::config(modeBarButtonsToRemove = setdiff(all.plotly.buttons,"toImage") ) %>%
            plotly::config(toImageButtonOptions = list(format='svg',
                                               height=800, width=800, scale=1.1)) %>%
            plotly::config(displaylogo = FALSE) %>% 
            plotly::event_register('plotly_selected') 

        dbg("scatterPlotMatrix:: done\n")
        p    
    })

    scatterPlotMatrix.opts = shiny::tagList(
        withTooltip( shiny::checkboxInput(ns("splom_highlight"),"Highlight genes",TRUE),
               "Enable highlighting genes on the plots. Users can highlight points by selecting them with the mouse, using the box selection or the lasso selection tool.")
        ##tipify( shiny::selectInput(ns("splom_ntop"),"Number of top genes",c(100,500,1000,2500,999999),selected=1000),
        ##        "Number of top genes ")
    )

    scatterPlotMatrix_info = "For the selected contrasts, the <strong>Pairs</strong> panel provides pairwise scatterplots for the differential expression profiles corresponding to multiple contrasts. The main purpose of this panel is to identify similarity or dissimilarity between selected contrasts. When K >= 3 contrasts are selected, the figure shows a KxK scatterplot matrix. When K <= 2, the Pairs panel provides an interactive pairwise scatterplots for the differential expression profiles of the two selected contrasts. The pairs plot is interactive and shows information of each gene with a mouse hover-over. Users can also select a number points by selecting points with the mouse, using the box selection or the lasso selection tool. Note that the selected genes will appear in input panel on the left sidebar as '<custom>' selection."

    shiny::callModule(
        plotModule,
        id = "scatterPlotMatrix", 
        func = scatterPlotMatrix.PLOT,
        plotlib = "plotly",
        title = "SCATTERPLOT PAIRS", label="a",
        options = scatterPlotMatrix.opts,
        ## download.fmt = c("pdf","html"),  ## scatterGL does not work for PDF
        ## download.fmt = c("html"),
        pdf.width=8, pdf.height=8,
        height = c(740,750), width = c('100%',1000),
        res=95,
        info.text = scatterPlotMatrix_info,
        ##caption = scatterPlotMatrix_caption,
        add.watermark = WATERMARK
    )
    ##output <- attachModule(output, scatterPlotMatrix_module)

    
    ##======================================================================
    ## Venn diagram
    ##======================================================================

    venndiagram.RENDER <- shiny::reactive({
        
        dt = getSignificanceCalls()
        if(is.null(dt) || nrow(dt)==0) return(NULL)
        
        dt1 = dt[,2:ncol(dt),drop=FALSE]
        label = LETTERS[1:ncol(dt1)]
        colnames(dt1) = label
        include = "both"
        if(input$include=="up/down") {
            include = c("up","down")
        }    
        
        par(mfrow=c(1,1), mar=c(1,1,3,1)*0, bty="n")
        par(oma=c(0.0,0,0,0))
        if(ncol(dt1)==1) {
            frame()
            text(0.5, 0.5, "Error: Venn diagram needs at least two groups")
            return(NULL)            
        } 
        if(ncol(dt1)>7) {
            frame()
            text(0.5, 0.5, "Error: too many groups for Venn diagram")
            return(NULL)
        } 

        p <- NULL
        if(0) {
            ## dt1 = dt1[,1:min(5,ncol(dt1))]
            limma::vennDiagram(
                       dt1,  main=NULL, cex.main=0.2, cex=1.2, mar=c(0,0,2,0),
                       include=include, bty="n", fg=grey(0.7),
                       circle.col=c("turquoise", "salmon","lightgreen","orange") )
            tt = paste(label,"=",colnames(dt)[-1])
            legend("topleft", legend=tt, bty='n', cex=0.9, y.intersp=0.95,
                       inset=c(0.04,-0.01), xpd=TRUE)       
            
        } else {

            ##colnames(dt1) <- colnames(dt)[-1]
            
            x <- apply(dt1, 2, function(x) rownames(dt1)[which(x!=0)])
            ##venntable <- ggVennDiagram::process_region_data(Venn(x))
            xlen <- sapply(x,length) 
            dbg("[venndiagram.RENDER] list.len = ",xlen)
            
            if(length(x)==0 || all(xlen==0)) {
                frame()
                text(0.5, 0.5, "Error: no valid genes. Please adjust thresholds.",
                     col='red')
                return(NULL)
            }            
            
            p <- ggVennDiagram::ggVennDiagram(
                x,
                label = "count",
                edge_size = 0.4
            ) +
                ggplot2::scale_fill_gradient(low="grey90",high = "red") +
                ggplot2::theme( legend.position = "none", 
                      plot.margin = ggplot2::unit(c(1,1,1,1)*0.3, "cm"))
            
            ## legend
            tt = paste(label,"=",colnames(dt)[-1])
            n1  <- ceiling(length(tt)/2)
            tt1 <- tt[1:n1]
            tt2 <- tt[(n1+1):length(tt)]
            if(length(tt2) < length(tt1)) tt2 <- c(tt2,"   ")
            tt1 <- paste(tt1, collapse='\n')
            tt2 <- paste(tt2, collapse='\n')            
            
            xlim <- ggplot2::ggplot_build(p)$layout$panel_scales_x[[1]]$range$range
            ylim <- ggplot2::ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
            x1 = xlim[1] - 0.1*diff(xlim)
            x2 = xlim[1] + 0.6*diff(xlim)
            y1 = ylim[2] + 0.12*diff(xlim)
            
            p <- p +
                ggplot2::annotate("text", x = x1, y = y1, hjust="left",
                         label = tt1, size=4, lineheight=0.83) +
                ggplot2::annotate("text", x = x2, y = y1, hjust="left",
                             label = tt2, size=4, lineheight=0.83) +
                ggplot2::coord_sf(clip="off")

            grid::grid.draw(p)
        }        
        ##p
    })

    FDR.VALUES2 <- c(1e-9,1e-6,1e-3,0.01,0.05,0.1,0.2,0.5,1)
    venndiagram.opts = shiny::tagList(
        ##    shiny::checkboxGroupInput(ns('intersection'),NULL, choices=c("A","B","C"), inline=TRUE )    
        shiny::fillRow(
            flex=c(1,1), ## height=80,       
            withTooltip( shiny::selectInput(ns("fdr"),"FDR", choices=FDR.VALUES2, selected=0.20),
                   "Threshold for false discovery rate",
                   placement="right", options = list(container = "body")),
            withTooltip( shiny::selectInput(ns("lfc"),"logFC threshold",
                                choices = c(0,0.1,0.2,0.5,1,2,5),
                                selected = 0.2),
                   "Threshold for fold-change (log2 scale)",
                   placement="right", options = list(container = "body"))
        ),
        shiny::br(),br(),br(),br(),
        shiny::radioButtons(ns('include'),'Counting:', choices=c("both","up/down"), inline=TRUE)
    )
   
    shiny::callModule(
        plotModule,
        id = "venndiagram", 
        func = venndiagram.RENDER,
        func2 = venndiagram.RENDER,
        ## plotlib = "ggplot",
        title = "VENN DIAGRAM", label="b",
        ##caption = venntable_buttons,
        info.text = "The Venn diagram visualizes the number of intersecting genes between the profiles. The list of intersecting genes with further details is also reported in an interactive table below, where users can select and remove a particular contrasts from the intersection analysis.",
        options = venndiagram.opts,
        pdf.width = 8, pdf.height = 8,
        height = c(400,700),
        width = c('100%',900),
        res = c(72,90),
        add.watermark = WATERMARK
    )

    ##output$intersection_table <- DT::renderDataTable({
    venntable.RENDER <- shiny::reactive({
        ngs <- inputData()
        shiny::req(ngs)
        
        ## get foldchanges
        fc0 = getSignificantFoldChangeMatrix()  ## isolate??        
        if(is.null(fc0) || nrow(fc0)==0) return(NULL)
                
        ## add gene name/title
        if(input$level == "gene") {
            gene = as.character(ngs$genes[rownames(fc0),"gene_name"])
            gene.tt = substring( GENE.TITLE[gene],1,50)
            gene.tt = as.character(gene.tt)
            ##fc0 = data.frame( name=name, title=gene.tt, fc0)
            fc0 = data.frame(name=gene, fc0, check.names=FALSE)
        } else {
            name = substring(rownames(fc0),1,50)
            name[is.na(name)] = "NA"
            fc0 = data.frame(name=name, fc0, check.names=FALSE)
        }
        
        df = data.frame( fc0, check.names=FALSE)
        nsc = setdiff(1:ncol(df),2)
        ##dt <- dt[rownames(fc0),]    
        ##D <- cbind(intersection=dt$intersection, D)
        DT::datatable(df,
                      class='compact cell-border stripe',
                      rownames=FALSE,
                      extensions = c('Scroller'), selection='none',
                      fillContainer = TRUE,
                      options=list(
                          ##dom = 'lfrtip',
                          dom = 'tip',                          
                          ## buttons = c('copy','csv','pdf'),
                          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          ## columnDefs = list(list(targets=nsc, searchable = FALSE)),
                          scrollX = TRUE, 
                          ##scrollY = 150,
                          scrollY = '70vh',
                          scroller = TRUE,
                          deferRender = TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  
    })

    venntable_buttons <- shiny::inputPanel(
        shiny::div(shiny::checkboxGroupInput(
            ns('intersection'), NULL,
            choices = c("A","B","C"), ## selected=c("A","B","C"),
            inline=TRUE ),
            style="font-size:0.85em; margin-top:-10px;")
    )

    venntable_opts <- shiny::tagList(
        shiny::selectInput(ns('venntable_intersection'),'Filter intersection:', choices=NULL)        
    )
    
    shiny::callModule(
        tableModule,
        id = "venntable", 
        func = venntable.RENDER,
        ##caption = venntable_buttons,
        options = venntable_opts,
        title = "INTERSECTION",
        label = "c",
        info.text = "Table of genes in selected intersection.", 
        ## info.width = "400px",
        height = c(260,750),
        width = c('auto',1200)
    )
    
    
    ##================================================================================
    ## Single-pair scatter plot
    ##================================================================================

    pairsPlot.PLOT <- shiny::reactive({



        
        ##res = pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        res = getActiveFoldChangeMatrix()
        ##if(is.null(res)) return(NULL)
        shiny::req(res)    
        fc0 = res$fc.full[,1:2]  
        fc1 = res$fc[,1:2]  
        
        ##fc0 <- fc0[order(-rowSums(fc0**2)),]
        fc0 <- fc0[order(-apply(abs(fc0),1,max)),]
        fc0 <- fc0[order(-rowMeans(fc0**2)),]
        fm0 <- rowMeans(fc0**2)
        
        ## selected genes (will be labeled)
        sel.genes <- head(rownames(fc0),50)
        pairs.sel <- plotly::event_data("plotly_selected", source="pairs")
        sel.keys <- as.character(pairs.sel$key)
        if(length(sel.keys) > 0) {
            sel.genes <- sel.keys
        }
        ## sel.genes <- head(sel.genes,50)
                
        ## subsample for speed: take top1000 + 1000
        df <- data.frame(fc0)
        if(0) {
            ntop = 99999
            ##ntop <- input$splom_ntop                    
            jj <- match(c(sel.genes,hi.genes), rownames(df))
            jj <- c(jj, 1:min(ntop,nrow(df)))
            if(nrow(df)>ntop) {
                nremain = setdiff(1:nrow(df),jj)
                jj = c(jj, sample(nremain,min(1000,length(nremain))))  ## add 1000 random
            }
            jj <- unique(jj)
            ##df <- data.frame(head(fc0,ntop))
            df <- data.frame(df[jj,])
        }
        dim(df)

        ## resort selection so that selected genes are drawn last to avoid
        ## covering them up.
        is.hi  = (rownames(df) %in% hi.genes)
        is.sel = (rownames(df) %in% sel.genes)

        df.color = c("#00000033","#0066FF")[1 + is.hi]
        df.color = c("#AAAAAA","#1e60BB")[1 + is.hi]
        df.color = c("#AAAAAA55","#1e60BB99")[1 + is.hi]
        df.color = c("#BBBBBB22","#2e80BB55","#1e60BBAA")[1 + is.hi + (is.hi & is.sel)]
        df.color = c("#BBBBBB55","#1e60BBAA")[1 + is.hi]
        ##df.color = c("#AAAAAA44","#1e60BB66","#1e60BB99","#1e60BBAA")[1 + 1*is.hi + 2*is.sel]
       
        ## Labels for top 50
        ##label.text = label.text0 = head(rownames(df)[which(is.sel)],50)
        which.lab <- which(is.sel & is.hi)
        which.lab <- head(which.lab[order(-fm0[rownames(df)[which.lab]])],20)
        label.text = label.text0 = rownames(df)[which.lab]
        label.text <- sub("","",label.text)  ## strip prefix??
        label.text <- shortstring(label.text,30)
        if(sum(is.na(label.text))) label.text[is.na(label.text)] <- ""
        width1 <- rep(0,length(is.sel))
        if(input$pairsplot_showselected) {
            width1 <- 0.3 + 0.5*is.sel
        }

        ## reorder so the selected genes don't get overlapped
        jj <- order(is.sel|is.hi + is.sel&is.hi)
        df <- df[jj,]
        df.color <- df.color[jj]
        width1 <- width1[jj]

        sel1 <- NULL
        if(input$pairsplot_labelgenes) {
            sel1 <- match(label.text0, rownames(df))  ## index for labeled
        }
        
        ## Tooltip text for all 
        tt <- rownames(df)  ## strip prefix
        ## tt <- sub("","",tt)  ## strip prefix??
        if(input$level == "gene") {
            ngs <- inputData()
            g <- rownames(df)
            tt <- paste0("<b>",g,"</b> ", ngs$genes[g,"gene_title"])
        }
        tt <- gsub("_"," ",tt)
        tt <- sapply(tt,breakstring2,50,brk="<br>")

        ## plotly
        ##   
        axis = list( showline=TRUE, zeroline=TRUE,
                    gridcolor='#dddf', ticklen=4 )
        
        rho = cor(df[,1], df[,2])        
        annot.rho <- list(
            text = paste("r=",round(rho,4)),
            font = list(size=14),
            align = 'left',
            showarrow = FALSE,
            xref = 'paper',
            yref = 'paper',
            x = 0.03,
            y = 0.97,
            borderpad = 8, 
            bordercolor = 'black',
            borderwidth = 0.6)
        
        p <- plotly::plot_ly( data=df[,1:2], x = df[,1], y= df[,2],
                     type = 'scattergl', mode = 'markers',
                     source='pairs', key = rownames(df),
                     ##type = 'scatter', mode = 'markers',
                     text = tt,
                     hovertemplate = paste0("<br>%{text}<br>x: %{x}<br>y: %{y}<extra></extra>"),
                     marker = list(
                         color = df.color,
                         size = 8,
                         line = list(
                             width = width1,
                             ##color = 'rgb(230,230,230)'
                             color = 'rgba(128,128,128,64)'
                         ))
                     )

        if(!is.null(sel1)) {
            p <- p %>%
                plotly::add_annotations(
                    x = df[sel1,1],
                    y = df[sel1,2],
                    text = as.character(label.text),
                    ##text = rep("mylabel",length(sel1)),
                    ##xanchor = 'left',
                    xanchor = 'center',
                    yanchor = 'top',
                    font = list(size=14),
                    xref = "x",
                    yref = "y",
                    showarrow = FALSE,
                    ax = 20,
                    ay = -40
                )
        }

        p <- p %>%
            plotly::layout(
                ## title= 'Scatterplot',
                annotations = annot.rho,
                hovermode = 'closest',
                dragmode= 'select',
                ##plot_bgcolor='rgba(240,240,240, 0.95)',
                ## template = "plotly_dark",
                xaxis = c(title = paste(colnames(df)[1],"   (logFC)"), axis),
                yaxis = c(title = paste(colnames(df)[2],"   (logFC)"), axis)
            )    

        p <- p %>%
            plotly::layout(margin = list(80,80,80,80) )  ## l,r,b,t
        
        p <- p %>%
            ## plotly::config(displayModeBar = FALSE) %>% ## disable buttons
            plotly::config( toImageButtonOptions = list(format='svg', height=800, width=800, scale=1.1)) %>%
            plotly::event_register('plotly_selected') 
        p    
    })

    pairsPlot.opts = shiny::tagList(
        withTooltip( shiny::checkboxInput(ns("pairsplot_labelgenes"),"Label genes",TRUE),
               "Label genes on the plots."),
        withTooltip( shiny::checkboxInput(ns("pairsplot_showselected"),"Show selected",TRUE),
               "Show selected genes.")
        ##tipify( shiny::selectInput(ns("pairs_ntop"),"Number of top genes",c(100,500,1000,2500,999999),selected=1000),
        ##        "Number of top genes ")
    )

    pairsPlot_info = "For the selected contrasts, the <strong>Pairs</strong> panel provides pairwise scatterplots for the differential expression profiles corresponding to multiple contrasts. The main purpose of this panel is to identify similarity or dissimilarity between selected contrasts. When K >= 3 contrasts are selected, the figure shows a KxK scatterplot matrix. When K <= 2, the Pairs panel provides an interactive pairwise scatterplots for the differential expression profiles of the two selected contrasts. The pairs plot is interactive and shows information of each gene with a mouse hover-over. Users can also select a number points by selecting points with the mouse, using the box selection or the lasso selection tool. Note that the selected genes will appear in input panel on the left sidebar as '<custom>' selection."

    shiny::callModule(
        plotModule,
        id = "pairsPlot", 
        func = pairsPlot.PLOT,
        plotlib="plotly",
        title = "Scatterplot (pairs)", label="a",
        options = pairsPlot.opts,
        ##  download.fmt = c("pdf","html"),  ## scatterGL does not work for PDF
        download.fmt = c("html"),
        pdf.width=8, pdf.height=8,
        height = c(fullH-80,700), res=95,
        info.text = pairsPlot_info,
        ##caption = pairsPlot_caption,
        add.watermark = WATERMARK
    )
        
    ##=============================================================================
    ## FOLDCHANGE HEATMAP
    ##=============================================================================
    
    FoldchangeHeatmap.RENDER <- shiny::reactive({
        ##
        ##
        ngs <- inputData()
        if(input$FoldchangeHeatmap_allfc) {
            F <- getFoldChangeMatrix()$fc
        } else {
            F <- getActiveFoldChangeMatrix()$fc
        }
        F <- F[order(-rowMeans(F**2)),]
        F <- F[order(-abs(rowMeans(F))),]

        F1 <- head(F,80)
        F1 <- F1[order(rowMeans(F1)),]
        bh=5;mh=6;bm=4;cclust = TRUE
        mh = min( max(ncol(F1)*0.35, 0.8), 8)
        cclust = input$FoldchangeHeatmap_cluster
        if(input$level=='geneset') {
            bh = 3
            ## mh = min(max(ncol(F1)*0.35, 0.8), 8)
        }
        bm = 10 - mh  ## bottom margin
        at <- input$FoldchangeHeatmap_annotype
        
        par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,3,0))
        plt <- grid::grid.grabExpr({
            frame()
            heatmapWithAnnot(
                F1, anno.type=at,
                bar.height=bh, map.height=mh,
                mar = c(bm,0.5,0.5,1),
                cluster_columns = cclust,
                ##column_dend_height = ggplot2::unit(10,"mm"),
                inset = c(0.01,0.01))            
        })
        ##grid::grid.draw(plt)
        plt        
    })
    
    FoldchangeHeatmap.opts = shiny::tagList(
        withTooltip( shiny::checkboxInput(ns('FoldchangeHeatmap_allfc'), "show all contrasts", TRUE),
               "Show all contrasts or just the selected ones."),
        withTooltip( shiny::checkboxInput(ns('FoldchangeHeatmap_cluster'), "cluster genes", FALSE),
               "Cluster genes (columns)."),
        withTooltip( shiny::radioButtons(ns('FoldchangeHeatmap_annotype'), "annotation type",
                             c("boxplot","barplot"), inline=TRUE),
               "Plot type of column annotation.")
        )
    
    FoldchangeHeatmap_info = "<b>The Connectivity Heatmap</b> shows the most similar profiles as a heatmap. Contrasts that are similar will be clustered close together."

    FoldchangeHeatmap_caption = "<b>Connectivity Heatmap.</b> Similarity of the contrasts profiles as a heatmap. Contrasts that are similar will be clustered close together."
    
    shiny::callModule(
        plotModule,
        "FoldchangeHeatmap", label = "a",
        func = FoldchangeHeatmap.RENDER,
        func2 = FoldchangeHeatmap.RENDER, 
        plotlib="grid",
        options = FoldchangeHeatmap.opts,
        title = "FOLDCHANGE HEATMAP",
        info.text = FoldchangeHeatmap_info,
        ##caption = FoldchangeHeatmap_caption,
        pdf.width=14, pdf.height=6.5,
        height = c(750,750), width = c('auto',1600),
        res = c(90,110),
        add.watermark = WATERMARK
    )

    ##================================================================================
    ## Contrast corrplot 
    ##================================================================================

    ctcorrplot.PLOT <- shiny::reactive({
                        
        ngs <- inputData()
        shiny::req(ngs)
        shiny::req(input$comparisons)
        
        ## res <- pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        res <- getFoldChangeMatrix()

        if(is.null(res)) return(NULL)
        ##validate(shiny::need(NCOL(res$fc)<2, "warning. need multiple comparisons."))
        if(NCOL(res$fc)<2) return(NULL)        
        
        fc0 = res$fc
        qv0 = res$qv

        ntop = 9999
        ntop <- input$ctcorrplot_ntop
        if(ntop=="all") ntop <- 999999
        ntop <- as.integer(ntop)

        allfc <- input$ctcorrplot_allfc
        if(!allfc) {
            comp = input_comparisons()
            if(length(comp)<2) return(NULL)
            kk = match(comp, colnames(fc0))
            fc0 <- fc0[,kk,drop=FALSE]
        }

        ##R.full <- cor(fc0[,], use="pairwise", method="spearman")
        R.full <- cor(apply(fc0,2,rank), use="pairwise")
        jj <- head(order(-rowMeans(fc0**2)),ntop)
        ##R <- cor(fc0[jj,], use="pairwise", method="spearman")
        R <- cor(apply(fc0[jj,],2,rank), use="pairwise")
        R <- round(R,digits=2)
        diag(R) <- 0
        
        notecex=0.001
        notecex=1.1; cex=1.3
        if( nrow(R) > 8)  {notecex=0.95; cex=1.2}
        if( nrow(R) > 20) {notecex=0.75; cex=1.05}
        if( nrow(R) > 50) {notecex=0.65; cex=0.9}
        if( nrow(R) > 80) {notecex=0.0001; cex=0.6}
        
        mar1 <- c(16,18)*1.2
        if(nrow(R) <= 8) { mar1=c(16,18)*2 }
        if(nrow(R) > 30) { mar1=c(16,18)*0.9 }
        if(nrow(R) > 80) { mar1=c(16,18)*0.6 }    
        
        col <- BLUERED(16)
        col <- gplots::colorpanel(64,"royalblue3","grey90","indianred3")
        ##col <- tail(BLUERED(16),8)
        if(min(R, na.rm=TRUE)>=0) col <- tail(col,32)
        if(max(R, na.rm=TRUE)<=0) col <- head(col,32)
        cellnote <- NULL

        col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                   "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                   "#4393C3", "#2166AC", "#053061"))
        

        corrplot::corrplot(R, method = "circle", order="hclust",
                 is.corr = FALSE,
                 cl.lim = c(-1.02,1.02)*max(abs(R),na.rm=TRUE),
                 col = rev(col2(50)),
                 ## mar = c(1,0.2,0.2,1)*0.2*mean(mar1),
                 mar = c(0,0,0,0),                 
                 tl.cex = 0.8*cex,
                 tl.col = "black",
                 tl.srt = 90)
        ##corrplot(R, method = "circle", order="AOE")
        ##return(R)
    })

    ctcorrplot.PLOTLY <- shiny::reactive({

        ## install.packages("heatmaply")

        
        ngs <- inputData()
        shiny::req(ngs)
        ##req(input$comparisons)
        ##res <- pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        res = getFoldChangeMatrix()
        if(is.null(res)) return(NULL)
        if(NCOL(res$fc)<2) return(NULL)

        fc0 = res$fc
        qv0 = res$qv

        ntop=2000
        ntop <- input$ctcorrplot_ntop
        if(ntop=="all") ntop <- 999999
        ntop <- as.integer(ntop)
        
        allfc <- input$ctcorrplot_allfc
        if(!allfc) {
            comp = input_comparisons()
            if(length(comp)<2) return(NULL)
            kk = match(comp, colnames(fc0))
            fc0 <- fc0[,kk,drop=FALSE]
        }

        ##R.full <- cor(fc0[,], use="pairwise", method="spearman")
        R.full <- cor(apply(fc0,2,rank), use="pairwise")
        jj <- head(order(-rowMeans(fc0**2)),ntop)
        ##R <- cor(fc0[jj,], use="pairwise", method="spearman")
        R <- cor(apply(fc0[jj,],2,rank), use="pairwise")
        R <- round(R,digits=2)
        
        col <- BLUERED(16)
        col <- gplots::colorpanel(64,"royalblue3","grey90","indianred3")
        ##col <- tail(BLUERED(16),8)
        if(min(R,na.rm=TRUE)>=0) col <- tail(col,32)
        if(max(R,na.rm=TRUE)<=0) col <- head(col,32)

        bluered.pal <- colorRampPalette(colors = c("royalblue3","grey90","indianred3"))
        cellnote <- NULL
        ##if(input$ctcorrplot_showrho) cellnote <- R

        plt <- heatmaply(
            R, margins = c(250, 200, NA, 0),
            ## k_col = 5, k_row = 5,
            cellnote = cellnote, cellnote_size = 11,
            cellnote_textposition = "middle center", 
            colors = bluered.pal,
            limits = c(-1,1))        

        plt
    })
    
    ctcorrplot.opts = shiny::tagList(
        ##tipify( shiny::checkboxInput(ns('ctcorrplot_showrho'), "show correlation values", FALSE),
        ##"Show correlation values in cells."),
        withTooltip( shiny::checkboxInput(ns('ctcorrplot_allfc'), "show all contrasts", TRUE),
               "Show all contrasts or just the selected ones."),
        ##tipify( shiny::checkboxInput('ctcorrplot_fixed', "fix heatmap", FALSE),
        ##       "Fix heatmap layout when changing number of top genes"),
        withTooltip( shiny::radioButtons(ns('ctcorrplot_ntop'), "number of top genes",
                             c("100","1000","all"),
                             selected="1000", inline=TRUE),
               "Number of top genes to compute correlation values.") 
    )

    ctcorrplot_info = "<strong>Constrast heatmap.</strong> Similarity of the contrasts visualized as a clustered heatmap. Contrasts that are similar will be clustered close together. The numeric value in the cells correspond to the Pearson correlation coefficient between contrast signatures. Red corresponds to positive correlation and blue to negative correlation."
        
    shiny::callModule(
        plotModule,
        id = "ctcorrplot",  label="b",
        func = ctcorrplot.PLOT, plotlib="base",
        func2 = ctcorrplot.PLOT,
        title = "CONTRAST CORRELATION",
        info.text = ctcorrplot_info,
        ##caption = ctcorrplot_caption,
        options = ctcorrplot.opts,
        download.fmt = c("pdf","png"),
        pdf.width = 11, pdf.height = 10,
        height = c(550,720), width = c("auto",1100),
        res=c(80,85),
        add.watermark = WATERMARK
    )
    
    ##================================================================================
    ## Gene/geneset UMAP plots
    ##================================================================================

    ctGeneUMAP.RENDER <- shiny::reactive({
        
        dbg("[ctGeneUMAP.RENDER] reacted!")
        out <- getCurrentSig()
        pos <- out$pos
        
        W <- out$sig
        sel0 <- input_comparisons()
        shiny::req(sel0)
        if(!all(sel0 %in% colnames(W))) return(NULL)
        W <- W[,sel0,drop=FALSE]        
        dim(W)

        sel1=sel2=1
        sel1 <- ctGeneTable_module$rows_selected()
        sel2 <- ctGseaTable_module$rows_selected()

        dbg("[ctGeneUMAP.RENDER] 1: ")
        
        is.single = (NCOL(W)==1)
        hilight = NULL
        hilight2 = NULL
        cexlab = 1.2
        cex = ifelse(ncol(W)>9, 0.5, 0.8)
        gset = NULL
        
        if(length(sel1)>0) {
            df <- getGeneTable()
            sel.gene = rownames(df)[1]            
            sel.gene = rownames(df)[sel1]
            hilight2 = sel.gene
            hilight = c(hilight,sel.gene)
            cexlab = 1.8
        }
        if(length(sel2)>0) {
            gse <- out$gsea
            gset = rownames(gse)[sel2]
            gset.genes <- unlist(getGSETS(gset))
            hilight = c(hilight, gset.genes)
            hilight2 = c(hilight2,hilight)
        }

        dbg("[ctGeneUMAP.RENDER] 2: length.hilight = ",length(hilight))
        dbg("[ctGeneUMAP.RENDER] 2: length.hilight2 = ",length(hilight2))
        
        no.hilight = (length(hilight)==0)
        is.nmf <- all(grepl("NMF",colnames(W)))
        ntop <- as.integer(input$ntop_genes)
        lab.pos <- NULL
        
        ## Spectral NMF vectors
        nc = ceiling(sqrt(ncol(W)))
        nr = ceiling(ncol(W)/nc)
        par(mfrow=c(nr,nc))
        par(mar=c(2.7,2.8,0.7,0.2), mgp=c(1.4,0.5,0), cex.axis=0.9, cex.lab=0.9)
        i=1
        for(i in 1:ncol(W)) {

            dbg("[ctGeneUMAP.RENDER] 3: i = ",i)
            dbg("[ctGeneUMAP.RENDER] 3: w.name = ",colnames(W)[i])

            jj <- match(rownames(pos), rownames(W))
            fc = W[jj,i]
            fc[is.na(fc)] <- 0
            
            if(no.hilight) {
                ## no selected gene or geneset
                sorted.fc <- sort(fc)
                hilight = c(head(names(sorted.fc),ntop/2),tail(names(sorted.fc),ntop/2))
                hilight2 = hilight
            }
            sorted.h2 <- hilight2[order(-fc[hilight2]**2)] 
            ## hilight2 = c(head(sorted.h2,ntop),tail(sorted.h2,ntop))
            hilight2 = head(sorted.h2,ntop)
            opacity = ifelse(length(hilight)>0 || length(hilight2)>0, 0.15, 1)
            zsym <- ifelse(is.nmf,FALSE,TRUE)
            if(min(fc,na.rm=TRUE)>=0) zsym <- FALSE
            
            dbg("[ctGeneUMAP.RENDER] 3: hilight = ",hilight)
            dbg("[ctGeneUMAP.RENDER] 3: hilight2 = ",hilight2)
            dbg("[ctGeneUMAP.RENDER] 3: length.hilight = ",length(hilight))
            dbg("[ctGeneUMAP.RENDER] 3: length.hilight2 = ",length(hilight2))
            dbg("[ctGeneUMAP.RENDER] 3: length.fc = ",length(fc))
            
            plt <- pgx.scatterPlotXY.BASE(
                pos, var=fc,
                lab.pos = lab.pos,
                softmax = 1,
                zsym = zsym,
                cex.lab = cexlab,
                opacity = opacity,
                cex = cex,
                cex.legend = 0.9,
                hilight.cex = 1.3*cex,
                hilight.col = NULL,
                hilight.lwd = 0.5,
                hilight = hilight,
                hilight2 = hilight2,
                set.par = FALSE,
                lab.xpd = FALSE,  ## clip repel
                dlim = c(0.05,0.05),
                bty = 'n',
                xlab = "UMAP-y  (genes)",
                ylab = "UMAP-y  (genes)",                
                legend.pos = 'bottomleft'
            )
            title(colnames(W)[i], cex.main=1.1, line=-0.5)
            if(!is.null(gset)) {
                tt <- substring(tolower(gset),1,80)
                title(tt, cex.main=0.8, line=-1.3, font.main=1)                
                lab.pos <- plt$lab.pos
            }
        }

        dbg("[ctGeneUMAP.RENDER] done!")
        
    })
        
    ctGeneUMAP_info = "<b>CORSA module analysis.</b> Functional analysis of NMF modules."

    ctGeneUMAP_opts = shiny::tagList(
        shiny::radioButtons(ns("ntop_genes"),"Gene labels:",
                     choices = c(0,10,30,100), selected = 10, inline=TRUE)
    )
    
    shiny::callModule(
        plotModule, 
        id = "ctGeneUMAP", ##ns=ns,
        title = "GENE SIGNATURE", label="a",
        func = ctGeneUMAP.RENDER,
        func2 = ctGeneUMAP.RENDER, 
        download.fmt = c("png","pdf"),
        options = ctGeneUMAP_opts,
        info.text = ctGeneUMAP_info,        
        height = c(480, 780), width = c('auto',1200),
        pdf.height = 8, pdf.width =12, 
        res = c(75,95),
        add.watermark = WATERMARK
    )

    ##-----------------------------------------------------------------
    ##----------------------  NMF gsea --------------------------------
    ##-----------------------------------------------------------------
    
    ctGseaUMAP.RENDER <- shiny::reactive({

        dbg("[ctGseaUMAP.RENDER] reacted")
        out <- getCurrentSig()
        gse <- out$gsea

        sel0=1:4
        sel0 <- input_comparisons()
        shiny::req(sel0)
        if(!all(sel0 %in% colnames(gse))) return(NULL)        
        gse <- gse[,sel0,drop=FALSE]        
        gse0 <- gse

        sel2 <- ctGseaTable_module$rows_selected()
        gset <- NULL
        if(length(sel2)>0) {
            gset <- rownames(gse)[sel2]
        }
        
        ii <- grep("HALLMARK",rownames(gse))
        ii <- ctGseaTable_module$rows_all()
        shiny::req(ii)
        gse <- gse[ii,,drop=FALSE]

        ## multiple group layout 
        gse.scores <- lapply(1:ncol(gse), function(i) gse[,i])
        names(gse.scores) <- colnames(gse)            

        ## layout
        ngse <- length(gse.scores)
        nc = ceiling(sqrt(ngse))
        nr = ceiling(ngse/nc)
        cex1 = 1.15
        if(nc==2 && nr>1) cex1 = 0.9
        if(nr==1) cex1 = 0.8
        par(mfrow=c(nr,nc))

        gstype <- input$gstype
        if(gstype=='bar') {
            par(oma=c(0,0,0,0))
            par(mar=c(2.8,2,1.8,3), mgp=c(1.6,0.6,0))            
            ntop = round(28/nr)
            ## maximum 24 terms
            gse.scores <- lapply(gse.scores, function(x)
                head(sort(x,decreasing=TRUE),ntop))
            xlim = c(0,max(abs(unlist(gse.scores))))        
            
            i=1
            for(i in 1:ngse) {
                gsea.barplot(
                    gse.scores[[i]],
                    names=names(gse.scores[[i]]),
                    n = ntop, xlim = xlim,
                    main='', cex.text=cex1)
                title(names(gse.scores)[i], cex.main=1.1, line=+0.3)            
            }
        } 

        if(gstype=='umap') {
            par(mar=c(2.7,2.8,0.7,0.2), mgp=c(1.4,0.5,0), cex.axis=0.9, cex.lab=0.9)
            pos  <- out$gsea.pos
            ntop <- as.integer(input$ntop_gsets)
            cex = ifelse(ngse>9, 0.5, 0.8)
            
            i=1
            for(i in 1:length(gse.scores)) {                
                k <- names(gse.scores)[i]
                if(ncol(gse)==1) {
                    var <- gse0[,1]
                } else {
                    var <- gse0[,k]
                }                
                if(!is.null(gset)) {
                    hmarks <- gset
                } else {
                    var1 <- gse.scores[[i]]
                    ss <- names(sort(var))
                    ss <- intersect(ss, names(var1))
                    hmarks <- c(head(ss,ntop/2),tail(ss,ntop/2))
                }
                var <- var[match(rownames(pos),names(var))]                
                opacity = ifelse(length(hmarks)>0, 0.15, 1)
                zsym <- ifelse(min(var,na.rm=TRUE)>=0, FALSE, TRUE)

                pgx.scatterPlotXY(
                    pos, var=var, 
                    zsym=zsym, set.par=FALSE, softmax=1,
                    cex=cex, cex.legend = 0.9, cex.lab=1.2, bty='n',
                    plotlib='base', col='grey70', dlim=c(0.2,0.08),
                    hilight=hmarks, hilight.col=NULL, opacity=opacity,
                    xlab = "UMAP-y  (genesets)", ylab="UMAP-y  (genesets)",
                    hilight.lwd=0.5, hilight.cex=1.3)
                title(k, cex.main=1.1, line=-0.5)                                
            }
        } 

        
    })
    
    ctGseaUMAP.info = "<b>Functional signature.</b> Functional enrichment signatures for selected comparisons are shown as barplots or UMAP plots."

    ctGseaUMAP.opts = shiny::tagList(
        shiny::radioButtons(ns("gstype"),"Geneset plot type:",
                     choices = c('bar','umap'),
                     selected = 'bar', inline=TRUE),
        shiny::radioButtons(ns("ntop_gsets"),"Geneset labels:",
                     choices = c(0,6,10,30,100), selected = 6, inline=TRUE)
    )
    
    shiny::callModule(
        plotModule, 
        id = "ctGseaUMAP", ##ns=ns,
        title = "FUNCTIONAL SIGNATURE", label="b",
        func  = ctGseaUMAP.RENDER,
        func2 = ctGseaUMAP.RENDER, 
        download.fmt = c("png","pdf"),
        options = ctGseaUMAP.opts,
        info.text = ctGseaUMAP.info,        
        height = c(480, 780), width = c('auto',1200),
        pdf.height = 8, pdf.width =12,
        res = c(75,95),
        add.watermark = WATERMARK
    )    

    ##-------------------------------------------
    ##--------------gene table ------------------
    ##-------------------------------------------

    getGeneTable <- shiny::reactive({

        ngs <- inputData()
        out <- getCurrentSig()        

        W <- out$sig
        sel0 <- 1:ncol(W)
        sel0 <- input_comparisons()
        shiny::req(sel0)
        if(length(sel0)==0) return(NULL)
        if(!all(sel0 %in% colnames(W))) return(NULL)

        ## only genes
        W <- W[rownames(W) %in% rownames(ngs$X),,drop=FALSE]
        W <- W[,sel0,drop=FALSE]

        tt <- NA
        tt <- GENE.TITLE[rownames(W)]
        tt <- substring(tt,1,80)
        df <- data.frame( gene=rownames(W), title = tt, W, check.names=FALSE)
        sel1 <- ctGseaTable_module$rows_selected()
        if(length(sel1)>0) {
            gset <- rownames(out$gsea)[sel1]
            gset.genes <- unlist(getGSETS(gset))
            gg <- intersect(rownames(df),gset.genes)
            df <- df[gg,,drop=FALSE]
        }        
        df
    })
    
    ctGeneTable.RENDER <- shiny::reactive({

        df <- getGeneTable()
        shiny::req(df)
        numeric.cols <- colnames(df)[3:ncol(df)]
        
        DT::datatable(
                df, rownames=FALSE, ## escape = c(-1,-2),
                ## filter = 'top',
                extensions = c('Buttons','Scroller'),
                selection = list(mode='single', target='row', selected=NULL),
                class = 'compact cell-border stripe hover',
                fillContainer = TRUE,
                options=list(
                    dom = 'lfrtip', ##buttons = c('copy','csv','pdf'),
                    ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                    scrollX = TRUE, ##scrollY = TRUE,
                    ##scrollY = 170,
                    scrollY = '70vh',
                    scroller=TRUE, deferRender=TRUE
                )  ## end of options.list 
            ) %>%
            DT::formatSignif(numeric.cols,3) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') 
    })

    ctGeneTable.info = "Genes."
    ctGeneTable.caption = "<b>Module gene table.</b> <b>(a)</b> ..."    
    ctGeneTable.opts <- shiny::tagList(
        ##selectInput(ns('wgcna_ctGeneTable_selmodule'),'module:', choices=NULL)
    )
    
    ctGeneTable_module <- shiny::callModule(
        tableModule, id = "ctGeneTable",
        title = "GENE TABLE", label="I",
        ## caption = wgcna_ctGeneTable_caption,
        func  = ctGeneTable.RENDER, ## ns=ns,
        options = ctGeneTable.opts,
        info.text = ctGeneTable.info,
        height = c(225,750), width=c('100%',1500)
    )
    
    ##-------------------------------------------
    ##-------------- GSEA table --------------
    ##-------------------------------------------
    
    ctGseaTable.RENDER <- shiny::reactive({

        ngs <- inputData()
        out <- getCurrentSig()        
        df <- out$gsea
        sel0 <- input_comparisons()
        shiny::req(sel0)
        if(!all(sel0 %in% colnames(df))) return(NULL)
        df <- df[,sel0,drop=FALSE]        
        numeric.cols <- colnames(df)
        gset <- substring(rownames(df),1,80)
        df <- data.frame(geneset=gset, df, check.names=FALSE)
        
        DT::datatable(
                df, rownames=FALSE, ## escape = c(-1,-2),
                ## filter = 'top',
                extensions = c('Buttons','Scroller'),
                selection = list(mode='single', target='row', selected=NULL),
                class = 'compact cell-border stripe hover',
                fillContainer = TRUE,
                options=list(
                    dom = 'lfrtip', ##buttons = c('copy','csv','pdf'),
                    ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                    scrollX = TRUE, ##scrollY = TRUE,
                    ##scrollY = 170,
                    scrollY = '70vh',
                    scroller=TRUE, deferRender=TRUE
                )  ## end of options.list 
            ) %>%
            DT::formatSignif(numeric.cols,3) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') 
    })

    ctGseaTable.info = "NMF GSEA."
    ctGseaTable.caption = "<b>Module enrichment table.</b> <b>(a)</b> ..."    
    ctGseaTable.opts <- shiny::tagList(
        ##selectInput(ns('wgcna_ctGseaTable_selmodule'),'module:', choices=NULL)
    )
    
    ctGseaTable_module <- shiny::callModule(
        tableModule, id = "ctGseaTable",
        title = "ENRICHMENT TABLE", label="II",
        ## caption = ctGseaTable_caption,
        func = ctGseaTable.RENDER, ## ns=ns,
        options = ctGseaTable.opts,
        info.text = ctGseaTable.info,
        height = c(225,750), width=c('100%',1500)
    )
  })    
} ## end-of-Board 
