IntersectionInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
    )
}

IntersectionUI <- function(id) {
    ns <- NS(id)  ## namespace
    fillCol(
        height = 750,
        if(DEV.VERSION) {
            tabsetPanel(
                id = ns("tabs1"),
                tabPanel("Pairs",uiOutput(ns("cmp_scatterPlotMatrix_UI"))),
                tabPanel("Pairs2",uiOutput(ns("cmp_pairsPlot_UI"))),
                tabPanel("Contrast heatmap",uiOutput(ns("cmp_ctClustering_UI")))
                ## tabPanel("Connectivity map",uiOutput(ns("cmp_connectivitymap_UI")))
            )
        } else {
            tabsetPanel(
                id = ns("tabs1"),
                tabPanel("Pairs",uiOutput(ns("cmp_scatterPlotMatrix_UI"))),
                tabPanel("Contrast heatmap",uiOutput(ns("cmp_ctClustering_UI")))
            )

        }
    )
    
}

IntersectionModule <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE
    fullH = 750       # row height of panel 
    
    ## reactive functions from shared environment
    inputData <- env[["load"]][["inputData"]]
    selected_gxmethods <- env[["expr"]][["selected_gxmethods"]]
    selected_gsetmethods <- env[["enrich"]][["selected_gsetmethods"]]
    
    description = "<b>Intersection analysis</b>. Compare experiments by intersecting
their signature genes. The main goal is to identify contrasts showing
similar profiles. Find genes that are commonly up/down regulated
between two contrasts."
    output$description <- renderUI(HTML(description))

    cmp_infotext =
        "The <strong>Intersection analysis module</strong> enables users to compare multiple contrasts by intersecting the genes of profiles. The main goal is to identify contrasts showing similar profiles.

<br><br>For the selected contrasts, the platform provides volcano plots and pairwise correlation plots between the profiles in the <strong>Pairs</strong> panel. Simultaneously, a Venn diagram with the number of intersecting genes between the profiles is plotted in <strong>Venn diagram</strong> panel. Details of intersecting genes are also reported in an interactive table. A more detailed scatter plot of two profiles is possible under the <strong>Two-pairs</strong> panel. Users can check the pairwise correlations of the contrasts under the <b>Contrast heatmap</b> panel. Alternatively, the <strong>Connectivity Map (CMap)</strong> shows the similarity of the contrasts profiles as a t-SNE plot.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=5' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>

"

    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================
    
    output$inputsUI <- renderUI({
        ui <- tagList(
            tipify( actionLink(ns("cmp_info"), "Tutorial", icon = icon("youtube")),
                   "Show more information about this module"),
            hr(), br(),             
            tipify( selectInput(ns('cmp_comparisons'),'Contrasts:', choices=NULL, multiple=TRUE),
                   "Select the contrasts that you want to compare. If you select N=2 contrast a single scatterplot will be drawn. For N>=3 a scatterplot matrix will be drawn.",
                   placement="top"),

            tipify( actionLink(ns("cmp_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Toggle advanced options.", placement="top"),
            br(),br(),
            conditionalPanel(
                "input.cmp_options % 2 == 1", ns=ns,
                tipify( selectInput(ns("cmp_filter"),"Filter:", choices=NULL, multiple=FALSE),
                       "Filter features", placement="top"),
                conditionalPanel(
                    "input.cmp_filter == '<custom>'", ns=ns,
                    tipify( textAreaInput(ns("cmp_customlist"), NULL, value = NULL,
                                          rows=5, placeholder="Paste your custom gene list"),
                           "Paste a custom list of genes to highlight.", placement="bottom")
                ),
                tipify( radioButtons(ns("cmp_level"),"Level:",
                                     choices=c("gene","geneset"), inline=TRUE),
                       "Select feature level: gene or geneset", placement="top")
            )
        )
        if(DEV.VERSION) {
            uix <- tagList(
                hr(),
                h5("Developer options:"),
                radioButtons(ns('cmp_featuretype'),'Feature type:',
                             choices=c("logFC","logCPM"), inline=TRUE)
            )
            ui <- c(ui, uix)
        }
        ui
    })
    outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!
    
    ## delayed input
    input_cmp_comparisons <- reactive({
        input$cmp_comparisons
    }) %>% debounce(500)
    
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    observeEvent( input$cmp_info, {
        showModal(modalDialog(
            title = HTML("<strong>Intersection Analysis Module</strong>"),
            HTML(cmp_infotext),
            easyClose = TRUE, size="l" ))
    })

    ## update choices upon change of data set 
    observe({
        ngs <- inputData()
        ##req(ngs)
        if(is.null(ngs)) return(NULL)
        comparisons <- colnames(ngs$model.parameters$contr.matrix)
        comparisons <- sort(comparisons)
        updateSelectInput(session, "cmp_comparisons", choices=comparisons,
                          selected=head(comparisons,3))
    })

    ## update choices upon change of feature level
    ##observeEvent( input$cmp_level, {
    observe({
        ngs <- inputData()
        ## req(ngs,input$cmp_level)
        if(is.null(ngs)) return(NULL)
        req(input$cmp_level)
        ##flt.choices = names(ngs$families)
        if(input$cmp_level=="geneset") {
            ft <- names(COLLECTIONS)
            nn <- sapply(COLLECTIONS, function(x) sum(x %in% rownames(ngs$gsetX)))
            ft <- ft[nn >= 10]
        } else {
            ## gene level
            ft <- pgx.getFamilies(ngs,nmin=10,extended=FALSE)
        }
        ft <- sort(ft)
        ## if(input$cmp_level=="gene") ft = sort(c("<custom>",ft))
        ## ft = sort(c("<custom>",ft))        
        updateSelectInput(session, "cmp_filter", choices=ft, selected="<all>")
    })
    
    observe({
        splom.sel <- plotly::event_data("plotly_selected", source="splom")
        sel.keys <- as.character(splom.sel$key)
        if(0 && length(sel.keys)>0) {
            updateSelectInput(session, "cmp_filter", selected="<custom>")
            sel.keys = paste(sel.keys, collapse=" ")
            updateTextAreaInput(session, "cmp_customlist", value=sel.keys)
        }
    })
    
    observeEvent(input$cmp_comparisons, {
        cmp <- input$cmp_comparisons
        if(is.null(cmp)) return(NULL)
        dt.labels = LETTERS[1:length(cmp)]
        updateCheckboxGroupInput(
            session, "cmp_intersection", choices=dt.labels,
            ## selected=dt.labels,
            inline=TRUE )    
    })
    
    ##================================================================================
    ##========================= REACTIVE FUNCTIONS ===================================
    ##================================================================================

    ## selected_gxmethods <- reactive({
    ##     ##sel <- SEL.GXMETHODS()
    ##     req(sel)
    ##     sel
    ## })
    
    getFoldChangeMatrix <- reactive({
        ## 
        ## Get full foldchange matrix from ngs object.
        ##
        ##
        ##
        ##dbg("<intersectionModule:getFoldChangeMatrix> reacted\n")
        fc0 = NULL
        qv0 = NULL
        ngs <- inputData()
        alertDataLoaded(session,ngs)        
        req(ngs)
        
        sel = names(ngs$gset.meta$meta)
        ##sel = input_cmp_comparisons()
        ##sel = intersect(sel, names(ngs$gset.meta$meta))
        ##if(length(sel)==0) return(NULL)
        
        if(input$cmp_level=="geneset") {
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
            if(input$cmp_filter == "<custom>") {
                gsets = strsplit( input$cmp_customlist, split="[, ;]")[[1]]
                if(length(gsets)>0) {
                    gsets = intersect(rownames(ngs$gsetX), gsets)
                }
            } else if(input$cmp_filter != "<all>") {
                gsets = unique(unlist(COLLECTIONS[input$cmp_filter]))
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
            if(input$cmp_filter == "<custom>") {
                genes = strsplit( input$cmp_customlist, split="[, ;]")[[1]]
                if(length(genes)>0) {
                    sel.probes = filterProbes(ngs$genes, genes)
                }
            } else  if(input$cmp_filter != "<all>") {
                sel.probes = filterProbes(ngs$genes, GSETS[[input$cmp_filter]])
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

    
    getActiveFoldChangeMatrix <- reactive({
        
        res = getFoldChangeMatrix()
        ##if(is.null(res)) return(NULL)
        req(res)            

        ## match with selected/active contrasts
        ## comp = head(colnames(res$fc),3)
        comp = input_cmp_comparisons()
        kk = match(comp, colnames(res$fc))
        if(length(kk)==0) return(NULL)
        if(length(kk)==1) kk = c(kk,kk)
        res$fc = res$fc[,kk,drop=FALSE]
        res$qv = res$qv[,kk,drop=FALSE]        
        res$fc.full = res$fc.full[,kk,drop=FALSE]
        res$qv.full = res$qv.full[,kk,drop=FALSE]        

        return(res)
    })
    
    getCPMMatrix <- reactive({
        ## 
        ## Get full foldchange matrix from ngs object.
        ##
        ##
        ##
        ##dbg("[IntersectionModule::getCPMMatrix] reacted\n")

        ngs <- inputData()
        req(ngs)
        
        if(input$cmp_level=="geneset") {
            grp <- ngs$samples$group
            mx = tapply(colnames(ngs$gsetX), grp, function(k) rowMeans(ngs$gsetX[,k,drop=FALSE]))
            mx = do.call(cbind, mx)
            mx0 = mx
            gsets = rownames(mx)
            if(input$cmp_filter != "<all>") {
                gsets = unique(unlist(COLLECTIONS[input$cmp_filter]))
            }
            gsets = intersect(gsets, rownames(mx))
            mx = mx[gsets,,drop=FALSE]
        } else {
            grp <- ngs$samples$group
            mx = tapply(colnames(ngs$X), grp, function(k) rowMeans(ngs$X[,k,drop=FALSE]))
            mx = do.call(cbind, mx)
            mx0 = mx
            sel.probes = rownames(mx) ## default to all probes
            if(input$cmp_filter == "<custom>") {
                genes = strsplit( input$cmp_customlist, split="[, ;]")[[1]]
                if(length(genes)>0) {
                    sel.probes = filterProbes(ngs$genes, genes)
                }
            } else if(input$cmp_filter %in% names(GSETS)) {
                sel.probes = filterProbes(ngs$genes, GSETS[[input$cmp_filter]])
            }
            sel.probes = intersect(sel.probes, rownames(mx))
            mx = mx[sel.probes,,drop=FALSE]
        }    
        mx <- mx[,!duplicated(colnames(mx)),drop=FALSE]
        res = list(mx=mx, mx.full=mx0)
        return(res) 
    })
    
    getSignificanceCalls <- reactive({
        ## Gets the matrix of significance calls.
        ##
        ngs <- inputData()

        sel = names(ngs$gset.meta$meta)
        sel = input_cmp_comparisons()
        sel = intersect(sel, names(ngs$gset.meta$meta))
        if(length(sel)==0) return(NULL)        
        res <- getFoldChangeMatrix()    
        fc <- res$fc[,sel,drop=FALSE]
        qv <- res$qv[,sel,drop=FALSE]

        fdr=0.05;lfc=0.2
        fdr = as.numeric(input$cmp_fdr)
        lfc = as.numeric(input$cmp_lfc)
        dt = sign(fc) * (qv <= fdr & abs(fc) >= lfc)
        dt[is.na(dt)] = 0
        ## add label of venn intersection region
        dt.labels = LETTERS[1:ncol(dt)]

        venn.intersection = apply( 1*(dt!=0), 1, function(x)
            paste(dt.labels[which(x==1)],collapse=""))
        dt = data.frame( intersection=venn.intersection, dt, check.names=FALSE )
        return(dt)
    })


    getSignificantFoldChangeMatrix <- reactive({
        ##
        ## Filters FC matrix with significance and user-defined
        ## intersection region.
        dt <- getSignificanceCalls()
        req(dt)

        isect = input$cmp_intersection
        fc0 = getFoldChangeMatrix()$fc  
        if( length(isect) == 0) {
            fc1 = fc0
        } else {
            ## only genes at least significant in one group
            jj = which(rowSums(dt[,2:ncol(dt),drop=FALSE]!=0)>0)
            if(length(jj)==0) return(NULL)
            dt = dt[jj,,drop=FALSE]

            ## check same sign
            if(input$cmp_include=="up/down") {
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
            intersection = paste0(input$cmp_intersection,collapse="")
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
        
        return(fc1)
    })

    ##======================================================================    
    ## Pairs (SPLOM)
    ##======================================================================

    ##======================================================================
    ## Scatterplot matrix in plotly
    ##
    ## From: https://plot.ly/r/splom/
    ##======================================================================
    
    cmp_scatterPlotMatrix.PLOT <- reactive({

        dbg("[IntersectionModule::cmp_scatterPlotMatrix.PLOT]  reacted\n")

        require(ggplot2)
        require(plotly)
        ##require(GGally)

        featuretype = "logFC"
        if(DEV.VERSION) {
            req(input$cmp_featuretype)
            featuretype <- input$cmp_featuretype
        }

        res=NULL
        if(featuretype=="logCPM") {
            ##res = pgx.getMetaFoldChangeMatrix(ngs, what="meta")
            res = getCPMMatrix()
            ##if(is.null(res)) return(NULL)
            req(res)    
            fc0 = res$mx.full
            fc1 = res$mx
        } else if(featuretype=="logFC") {
            ##res = pgx.getMetaFoldChangeMatrix(ngs, what="meta")
            res = getActiveFoldChangeMatrix()
            fc0 = res$fc.full
            fc1 = res$fc
        }
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
            ##ntop <- input$cmp_splom_ntop        
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
        if(input$cmp_splom_highlight) {
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
        if(input$cmp_level == "gene") {
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
        
        if(ncol(df)>=3) {
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
            
            p <- plot_ly(df, source="splom", key=rownames(df) ) %>%
                add_trace(
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
                add_annotations(
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
                layout(
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
            ## %>% style(diagonal = list(visible = F))

        } else {

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

            p <- plot_ly( data=df[,1:2], x = df[,1], y= df[,2],
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
                add_annotations(
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
                layout(
                    ## title= 'Scatterplot',
                    annotations = annot.rho,
                    hovermode = 'closest',
                    dragmode= 'select',
                    ##plot_bgcolor='rgba(240,240,240, 0.95)',
                    ## template = "plotly_dark",
                    xaxis = c(title = paste(colnames(df)[1],"   (logFC)"), axis),
                    yaxis = c(title = paste(colnames(df)[2],"   (logFC)"), axis)
                ) 
        }

        p <- p %>%
            layout(margin = list(80,80,80,80) )  ## l,r,b,t

        p <- p %>%
            ##config(displayModeBar = FALSE) %>% ## disable buttons
            config(modeBarButtonsToRemove = setdiff(all.plotly.buttons,"toImage") ) %>%
            config(toImageButtonOptions = list(format='svg',
                                               height=800, width=800, scale=1.1)) %>%
            config(displaylogo = FALSE) %>% 
            event_register('plotly_selected') 

        dbg("cmp_scatterPlotMatrix:: done\n")
        p    
    })

    cmp_scatterPlotMatrix.opts = tagList(
        tipify( checkboxInput(ns("cmp_splom_highlight"),"Highlight genes",TRUE),
               "Enable highlighting genes on the plots. Users can highlight points by selecting them with the mouse, using the box selection or the lasso selection tool.")
        ##tipify( selectInput(ns("cmp_splom_ntop"),"Number of top genes",c(100,500,1000,2500,999999),selected=1000),
        ##        "Number of top genes ")
    )

    cmp_scatterPlotMatrix_info = "For the selected contrasts, the <strong>Pairs</strong> panel provides pairwise scatterplots for the differential expression profiles corresponding to multiple contrasts. The main purpose of this panel is to identify similarity or dissimilarity between selected contrasts. When K >= 3 contrasts are selected, the figure shows a KxK scatterplot matrix. When K <= 2, the Pairs panel provides an interactive pairwise scatterplots for the differential expression profiles of the two selected contrasts. The pairs plot is interactive and shows information of each gene with a mouse hover-over. Users can also select a number points by selecting points with the mouse, using the box selection or the lasso selection tool. Note that the selected genes will appear in input panel on the left sidebar as '<custom>' selection."

    callModule(
        plotModule,
        id = "cmp_scatterPlotMatrix", 
        func = cmp_scatterPlotMatrix.PLOT,
        plotlib = "plotly",
        title = "Scatterplot matrix (pairs)", label="a",
        options = cmp_scatterPlotMatrix.opts,
        ## download.fmt = c("pdf","html"),  ## scatterGL does not work for PDF
        ## download.fmt = c("html"),
        pdf.width=8, pdf.height=8,
        height = c(fullH-80,700), res=95,
        info.text = cmp_scatterPlotMatrix_info
        ##caption = cmp_scatterPlotMatrix_caption
    )
    ##output <- attachModule(output, cmp_scatterPlotMatrix_module)

    
    ##======================================================================
    ## Venn diagram
    ##======================================================================

    cmp_venndiagram.RENDER %<a-% reactive({
        
        dt = getSignificanceCalls()
        if(is.null(dt) || nrow(dt)==0) return(NULL)
        
        dt1 = dt[,2:ncol(dt),drop=FALSE]
        label = LETTERS[1:ncol(dt1)]
        colnames(dt1) = label
        if(ncol(dt1)==1) {
            dt1 <- cbind(dt1, dt1)
        }
        include = "both"
        if(input$cmp_include=="up/down") {
            include = c("up","down")
        }    
        dt1 = dt1[,1:min(5,ncol(dt1))]
        
        par(mfrow=c(1,1), mar=c(1,1,3,1)*0, bty="n")
        par(oma=c(0.0,0,0,0))
        limma::vennDiagram(
                   dt1,  main=NULL, cex.main=0.2, cex=1.2, mar=c(0,0,2,0),
                   include=include, bty="n", fg=grey(0.7),
                   circle.col=c("turquoise", "salmon","lightgreen","orange") )
        tt = paste(label,"=",colnames(dt)[-1])
        legend("topleft", legend=tt, bty='n', cex=0.9, y.intersp=0.95,
               inset=c(0.04,-0.01), xpd=TRUE)
        
    })

    ##output$intersection_table <- DT::renderDataTable({
    cmp_venntable.RENDER <- reactive({
        ngs <- inputData()
        req(ngs)
        
        ## get foldchanges
        fc0 = getSignificantFoldChangeMatrix()  ## isolate??
        
        if(is.null(fc0) || nrow(fc0)==0) return(NULL)
        
        fc0 <- fc0[order(-rowMeans(fc0)),,drop=FALSE]        
        fc0 = round(fc0, digits=3)
        colnames(fc0) = paste0("fc.",LETTERS[1:ncol(fc0)])
        ##fc0 = data.frame(fc0)
        
        ## add gene name/title
        if(input$cmp_level == "gene") {
            gene = as.character(ngs$genes[rownames(fc0),"gene_name"])
            gene.tt = substring( GENE.TITLE[gene],1,50)
            gene.tt = as.character(gene.tt)
            ##fc0 = data.frame( name=name, title=gene.tt, fc0)
            fc0 = data.frame( name=gene, fc0, check.names=FALSE)
        } else {
            name = substring(rownames(fc0),1,50)
            name[is.na(name)] = "NA"
            fc0 = data.frame(name=name, fc0, check.names=FALSE)
        }
        
        df = data.frame( fc0, check.names=FALSE)
        ##dt <- dt[rownames(fc0),]    
        ##D <- cbind(intersection=dt$intersection, D)
        DT::datatable(df, class='compact cell-border stripe',
                      rownames=FALSE,
                      extensions = c('Scroller'), selection='none',
                      options=list(
                          ## pageLength = 40, ##lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE, ## scrollY = TRUE,
                          scrollY = 180, scroller=TRUE, deferRender=TRUE,
                          dom = 'lfrtip'                      
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  
    })

    FDR.VALUES2 <- c(1e-9,1e-6,1e-3,0.01,0.05,0.1,0.2,0.5,1)
    cmp_venndiagram.opts = tagList(
        ##    checkboxGroupInput(ns('cmp_intersection'),NULL, choices=c("A","B","C"), inline=TRUE )    
        fillRow(
            flex=c(1,1), ## height=80,       
            tipify( selectInput(ns("cmp_fdr"),"FDR", choices=FDR.VALUES2, selected=0.20),
                   "Threshold for false discovery rate",
                   placement="top", options = list(container = "body")),
            tipify( selectInput(ns("cmp_lfc"),"logFC threshold",
                                choices=c(0,0.2,0.5,1,2,5),
                                selected=0.5),
                   "Threshold for fold-change (log2 scale)",
                   placement="top", options = list(container = "body"))
        ),
        br(),br(),br(),br(),
        radioButtons(ns('cmp_include'),'Counting:', choices=c("both","up/down"), inline=TRUE)
    )
   
    cmp_venntable_buttons <- inputPanel(
        div(checkboxGroupInput(
            ns('cmp_intersection'), NULL,
            choices=c("A","B","C"), ## selected=c("A","B","C"),
            inline=TRUE ),
            style="font-size:0.85em; margin-top:-10px;")
    )

    callModule(
        plotModule,
        id = "cmp_venndiagram", 
        func = cmp_venndiagram.RENDER,
        func2 = cmp_venndiagram.RENDER,
        title = "Venn diagram", label="b",
        ##caption = cmp_venntable_buttons,
        info.text = "The Venn diagram visualizes the number of intersecting genes between the profiles. The list of intersecting genes with further details is also reported in an interactive table below, where users can select and remove a particular contrasts from the intersection analysis.",
        options = cmp_venndiagram.opts,
        pdf.width=8, pdf.height=8,
        height = 0.40*fullH, res=72
    )

    callModule(
        tableModule,
        id = "cmp_venntable", 
        func = cmp_venntable.RENDER,
        ##caption = cmp_venntable_buttons,
        title = "Intersecting genes", label="d",
        info.text = "Table of intersecting genes", 
        info.width = "500px",
        height = 0.4*fullH
    )

    ##================================================================================
    ## Cumuative FC
    ##================================================================================

    cmp_cumFCplot.RENDER %<a-% reactive({

        ngs <- inputData()
        sel = names(ngs$gx.meta$meta)
        req(input$cmp_comparisons)
        sel = input_cmp_comparisons()
        if(is.null(sel) || length(sel)==0 || sel[1]=="") return(NULL)
        sel = intersect(sel, names(ngs$gx.meta$meta))
        
        ##fc = sapply(ngs$gx.meta$meta[1:3], function(x) x$meta.fx)
        ##rownames(fc) <- rownames(ngs$gx.meta$meta[[1]])    
        fc = getSignificantFoldChangeMatrix()  ## isolate??xs
        fc <- fc[,sel,drop=FALSE]
        if(input$cmp_cumFCplot_abs) {
            fc <- abs(fc)  
        }
        fc[is.na(fc)] <- 0
        fc <- fc[order(-rowMeans(fc**2,na.rm=TRUE)),,drop=FALSE]
        
        ## add some empty rows (keeps barplot bar-widths equal)
        fc.na <- matrix(0,nrow=100,ncol=ncol(fc))
        fc <- rbind(fc,fc.na)
        col1 = grey.colors(ncol(fc),start=0.15)

        if(input$cmp_level=="geneset") {
            fc.top <- head(fc,24)
            fc.top <- fc.top[order(rowMeans(fc.top,na.rm=TRUE)),,drop=FALSE]            
            rownames(fc.top) <- tolower(rownames(fc.top))
            rownames(fc.top) <- gsub(".*[:]","",rownames(fc.top))
            ## rownames(fc.top) <- gsub("[ ]*\\(.*","",rownames(fc.top))
            par(mar=c(4,4,2,2))
            par(mfrow=c(1,3), mar=c(7.8,0,0,1), mgp=c(2.1,0.8,0) )
            par(mfrow=c(1,3), mar=c(5,0,0,1), mgp=c(2.1,0.8,0) )
            ##barplot(t(fc.top), las=3, cex.names=0.81, col=col1,
            ##        ylim=ylim, ylab="cumulative logFC")
            frame()
            frame()
            pgx.stackedBarplot(
                fc.top, cex.names=1.3, col=col1,
                hz=TRUE,  ## ylim=ylim,
                xlab="cumulative logFC")
            legend("topleft", legend=colnames(fc.top),
                   fill = col1, cex=1.0, y.intersp=0.8)
        } else {
            NTOP = 40
            fc.top <- head(fc,NTOP)
            fc.top <- fc.top[order(rowMeans(fc.top,na.rm=TRUE)),,drop=FALSE]            
            par(mar=c(8,4,2,2))
            par(mfrow=c(1,1), mar=c(9,4,1,1), mgp=c(2.4,1,0) )
            ##barplot(t(fc.top), las=3, cex.names=0.81, col=col1,
            ##        ylim=ylim, ylab="cumulative logFC")
            pgx.stackedBarplot(fc.top, cex.names=0.77, col=col1,
                               ## ylim=ylim,
                               ylab="cumulative logFC")        
            legend("topleft", legend=colnames(fc.top),
                   fill = col1, cex=0.75, y.intersp=0.78)
        }
        
    })

    cmp_cumFCplot.opts = tagList(
        checkboxInput(ns('cmp_cumFCplot_abs'),'Absolute foldchange')
    )

    callModule(
        plotModule,
        id = "cmp_cumFCplot", label="c",
        func = cmp_cumFCplot.RENDER,
        func2 = cmp_cumFCplot.RENDER,
        csvFunc = getSignificantFoldChangeMatrix,
        download.fmt = c("pdf","png"),
        ##caption = cmp_venntable_buttons,
        title = "Cumulative fold-change", 
        info.text = "</b>Cumulative fold-change.</b> This plot visualizes the cumulative fold-change of genes shared between the profiles.",
        options = cmp_cumFCplot.opts,
        pdf.width = 8, pdf.height = 6,
        height = c(0.45*fullH,550), width=c('auto',1000),
        res=c(75,110)
    )
    
    ##================================================================================
    ## Single-pair scatter plot
    ##================================================================================

    cmp_pairsPlot.PLOT <- reactive({
        require(ggplot2)
        require(plotly)
        ##require(GGally)
        
        ##res = pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        res = getActiveFoldChangeMatrix()
        ##if(is.null(res)) return(NULL)
        req(res)    
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
        
        ## highlighted genes (color blue)
        sel.row = pairsEnrichmentTable$rows_selected()
        hi.genes <- rownames(fc0)
        if(length(sel.row)>0) {
            df <- getPairsEnrichmentTable()
            gs <- rownames(df)[sel.row]
            hi.genes <- GSETS[[gs]]
        }
        
        ## subsample for speed: take top1000 + 1000
        df <- data.frame(fc0)
        if(0) {
            ntop = 99999
            ##ntop <- input$cmp_splom_ntop                    
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
        if(input$cmp_pairsplot_showselected) {
            width1 <- 0.3 + 0.5*is.sel
        }

        ## reorder so the selected genes don't get overlapped
        jj <- order(is.sel|is.hi + is.sel&is.hi)
        df <- df[jj,]
        df.color <- df.color[jj]
        width1 <- width1[jj]

        sel1 <- NULL
        if(input$cmp_pairsplot_labelgenes) {
            sel1 <- match(label.text0, rownames(df))  ## index for labeled
        }
        
        ## Tooltip text for all 
        tt <- rownames(df)  ## strip prefix
        ## tt <- sub("","",tt)  ## strip prefix??
        if(input$cmp_level == "gene") {
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
        
        p <- plot_ly( data=df[,1:2], x = df[,1], y= df[,2],
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
                add_annotations(
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
            layout(
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
            layout(margin = list(80,80,80,80) )  ## l,r,b,t
        
        p <- p %>%
            ## config(displayModeBar = FALSE) %>% ## disable buttons
            config( toImageButtonOptions = list(format='svg', height=800, width=800, scale=1.1)) %>%
            event_register('plotly_selected') 
        p    
    })

    cmp_pairsPlot.opts = tagList(
        tipify( checkboxInput(ns("cmp_pairsplot_labelgenes"),"Label genes",TRUE),
               "Label genes on the plots."),
        tipify( checkboxInput(ns("cmp_pairsplot_showselected"),"Show selected",TRUE),
               "Show selected genes.")
        ##tipify( selectInput(ns("cmp_pairs_ntop"),"Number of top genes",c(100,500,1000,2500,999999),selected=1000),
        ##        "Number of top genes ")
    )

    cmp_pairsPlot_info = "For the selected contrasts, the <strong>Pairs</strong> panel provides pairwise scatterplots for the differential expression profiles corresponding to multiple contrasts. The main purpose of this panel is to identify similarity or dissimilarity between selected contrasts. When K >= 3 contrasts are selected, the figure shows a KxK scatterplot matrix. When K <= 2, the Pairs panel provides an interactive pairwise scatterplots for the differential expression profiles of the two selected contrasts. The pairs plot is interactive and shows information of each gene with a mouse hover-over. Users can also select a number points by selecting points with the mouse, using the box selection or the lasso selection tool. Note that the selected genes will appear in input panel on the left sidebar as '<custom>' selection."

    callModule(
        plotModule,
        id = "cmp_pairsPlot", 
        func = cmp_pairsPlot.PLOT,
        plotlib="plotly",
        title = "Scatterplot (pairs)", label="a",
        options = cmp_pairsPlot.opts,
        ##  download.fmt = c("pdf","html"),  ## scatterGL does not work for PDF
        download.fmt = c("html"),
        pdf.width=8, pdf.height=8,
        height = c(fullH-80,700), res=95,
        info.text = cmp_pairsPlot_info
        ##caption = cmp_pairsPlot_caption
    )

    
    cmp_pairsSelectedTable.RENDER <- reactive({

        ngs <- inputData()
        req(ngs)
        
        ## get foldchanges
        ##fc0 = getSignificantFoldChangeMatrix()  ## isolate??
        fc0 = getActiveFoldChangeMatrix()$fc  ## isolate??        
        if(is.null(fc0) || nrow(fc0)==0) return(NULL)
        
        fc0 <- fc0[order(-rowMeans(fc0)),,drop=FALSE]        
        fc0 = round(fc0, digits=3)
        colnames(fc0) = paste0("fc.",LETTERS[1:ncol(fc0)])
        ##fc0 = data.frame(fc0)
        
        pairs.sel <- plotly::event_data("plotly_selected", source="pairs")
        sel.keys <- as.character(pairs.sel$key)
        if(0 && length(sel.keys)>0) {
            updateSelectInput(session, "cmp_filter", selected="<custom>")
            sel.keys = paste(sel.keys, collapse=" ")
            updateTextAreaInput(session, "cmp_customlist", value=sel.keys)
        }
        
        df = data.frame( fc0, check.names=FALSE)
        DT::datatable(df, class='compact cell-border stripe',
                      rownames = TRUE,
                      extensions = c('Scroller'),
                      selection='none',
                      options=list(
                          ## pageLength = 40, ##lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE, ## scrollY = TRUE,
                          scrollY = 225, scroller=TRUE, deferRender=TRUE,
                          dom = 'lfrtip'                      
                      )  ## end of options.list 
                      ) %>%
            ##formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  
    })

      
    getPairsEnrichmentTable <- reactive({

        ngs <- inputData()
        req(ngs)
        
        ## get foldchanges
        ##fc0 = getSignificantFoldChangeMatrix()  ## isolate??
        fc0 = getActiveFoldChangeMatrix()$fc  ## isolate??        
        if(is.null(fc0) || nrow(fc0)==0) return(NULL)

        ##sel.genes <- names(which(GSETxGENE[42174,]!=0))        
        pairs.sel <- plotly::event_data("plotly_selected", source="pairs")
        sel.genes <- as.character(pairs.sel$key)
        if(length(sel.genes) == 0) return(NULL)
                    
        ## fisher test
        cat("[cmp_pairsEnrichmentTable.RENDER] head(sel.genes)=",head(sel.genes),"\n")
        ii <- setdiff(match(toupper(sel.genes), colnames(GSETxGENE)),NA)
        N <- cbind(k1=Matrix::rowSums(GSETxGENE!=0), n1=ncol(GSETxGENE),
                   k2=Matrix::rowSums(GSETxGENE[,ii]!=0), n2=length(ii) )
        rownames(N) = rownames(GSETxGENE)
        ##N <- N[which(!(N[,1]==0 & N[,3]==0)), ]
        N <- N[which(N[,1]>0 | N[,3]>0), ]
        odds.ratio = ( N[,3]/ N[,4]) / ( N[,1]/ N[,2]) 
        dim(N)
        
        require(corpora) 
        pv <- corpora::fisher.pval( N[,1], N[,2], N[,3], N[,4], log.p=FALSE)
        head(pv)
        names(pv) <- rownames(N)
        pv = pv[match(names(odds.ratio),names(pv))]
        ##qv = p.adjust(pv, method="fdr")
        qv = p.adjust(pv, method="bonferroni")
        gset <- substring(names(pv),1,60)
        df = data.frame( geneset=gset, odds.ratio=odds.ratio, p.fisher=pv, q.fisher=qv)
        dim(df)

        ##df <- round(df, digits=3)
        df <- df[order(df$p.fisher),]
        df
    })
    
    cmp_pairsEnrichmentTable.RENDER <- reactive({

        df <- getPairsEnrichmentTable()
        req(df)
        numeric.cols <- 2:ncol(df)

        DT::datatable(
                df, class='compact cell-border stripe',
                rownames = FALSE,
                fillContainer = TRUE,
                extensions = c('Scroller'), ##selection='none',
                selection = list(mode='single', target='row', selected=1),                
                options=list(
                    ## pageLength = 40, ##lengthMenu = c(20, 30, 40, 60, 100, 250),
                    scrollX = TRUE, ## scrollY = TRUE,
                    scrollY = 225, scroller=TRUE, deferRender=TRUE,
                    dom = 'lfrtip'                      
                )  ## end of options.list 
            ) %>%
            formatSignif(numeric.cols,3) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  
    })
    
    pairsSelectedTable <- callModule(
        tableModule,
        id = "cmp_pairsSelectedTable", 
        func = cmp_pairsSelectedTable.RENDER,
        ##caption = cmp_pairsSelected_buttons,
        title = "Selected genes", label="b",
        info.text = "Table of selected genes", 
        info.width = "500px",
        height = 0.45*fullH
    )
    
    pairsEnrichmentTable <- callModule(
        tableModule,
        id = "cmp_pairsEnrichmentTable", 
        func = cmp_pairsEnrichmentTable.RENDER,
        ##caption = cmp_pairsSelected_buttons,
        title = "Functional enrichment", label="c",
        info.text = "Table", 
        info.width = "500px",
        height = 0.45*fullH
    )

    
    ##================================================================================
    ## Contrast heatmap 
    ##================================================================================

    cmp_ctheatmap.PLOT %<a-% reactive({
                        
        ngs <- inputData()
        req(ngs)
        req(input$cmp_comparisons)
        
        res <- pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        res = getFoldChangeMatrix()

        if(is.null(res)) return(NULL)
        ##validate(need(NCOL(res$fc)<2, "warning. need multiple comparisons."))
        if(NCOL(res$fc)<2) return(NULL)        
        
        fc0 = res$fc
        qv0 = res$qv

        ntop = 9999
        ntop <- input$cmp_ctheatmap_ntop
        if(ntop=="all") ntop <- 999999
        ntop <- as.integer(ntop)

        allfc <- input$cmp_ctheatmap_allfc
        if(!allfc) {
            comp = input_cmp_comparisons()
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
        col <- colorpanel(64,"royalblue3","grey90","indianred3")
        ##col <- tail(BLUERED(16),8)
        if(min(R,na.rm=TRUE)>=0) col <- tail(col,32)
        if(max(R,na.rm=TRUE)<=0) col <- head(col,32)
        cellnote <- NULL


        col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                   "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                   "#4393C3", "#2166AC", "#053061"))        
        require(corrplot)
        corrplot(R, method = "circle", order="hclust",
                 col = rev(col2(50)), mar=c(1,0.2,0.2,1) * 0.2*mean(mar1),
                 tl.cex = 0.65*cex, tl.col="black", tl.srt = 90)
        ##corrplot(R, method = "circle", order="AOE")
        ##return(R)
    })

    cmp_ctheatmap.PLOTLY <- reactive({

        ## install.packages("heatmaply")
        require(heatmaply)
        
        ngs <- inputData()
        req(ngs)
        ##req(input$cmp_comparisons)
        ##res <- pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        res = getFoldChangeMatrix()
        if(is.null(res)) return(NULL)
        if(NCOL(res$fc)<2) return(NULL)

        fc0 = res$fc
        qv0 = res$qv

        ntop=2000
        ntop <- input$cmp_ctheatmap_ntop
        if(ntop=="all") ntop <- 999999
        ntop <- as.integer(ntop)
        
        allfc <- input$cmp_ctheatmap_allfc
        if(!allfc) {
            comp = input_cmp_comparisons()
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
        col <- colorpanel(64,"royalblue3","grey90","indianred3")
        ##col <- tail(BLUERED(16),8)
        if(min(R,na.rm=TRUE)>=0) col <- tail(col,32)
        if(max(R,na.rm=TRUE)<=0) col <- head(col,32)

        bluered.pal <- colorRampPalette(colors = c("royalblue3","grey90","indianred3"))
        cellnote <- NULL
        ##if(input$cmp_ctheatmap_showrho) cellnote <- R

        plt <- heatmaply(
            R, margins = c(250, 200, NA, 0),
            ## k_col = 5, k_row = 5,
            cellnote = cellnote, cellnote_size = 11,
            cellnote_textposition = "middle center", 
            colors = bluered.pal,
            limits = c(-1,1))        

        plt
    })
    
    cmp_ctheatmap.opts = tagList(
        ##tipify( checkboxInput(ns('cmp_ctheatmap_showrho'), "show correlation values", FALSE),
        ##"Show correlation values in cells."),
        tipify( checkboxInput(ns('cmp_ctheatmap_allfc'), "show all contrasts", TRUE),
               "Show all contrasts or just the selected ones."),
        ##tipify( checkboxInput('cmp_ctheatmap_fixed', "fix heatmap", FALSE),
        ##       "Fix heatmap layout when changing number of top genes"),
        tipify( radioButtons(ns('cmp_ctheatmap_ntop'), "number of top genes",
                             c("100","1000","all"),
                             selected="1000", inline=TRUE),
               "Number of top genes to compute correlation values.") 
    )

    cmp_ctheatmap_info = "<strong>Constrast heatmap.</strong> Similarity of the contrasts visualized as a clustered heatmap. Contrasts that are similar will be clustered close together. The numeric value in the cells correspond to the Pearson correlation coefficient between contrast signatures. Red corresponds to positive correlation and blue to negative correlation."
        
    callModule(
        plotModule,
        id = "cmp_ctheatmap",  label="a",
        func = cmp_ctheatmap.PLOT, plotlib="base",
        func2 = cmp_ctheatmap.PLOT,
        ##func = cmp_ctheatmap.PLOTLY, plotlib="plotly",
        ##cmp_ctheatmap.PLOTLY, plotlib="generic", renderFunc="renderPlotly",
        info.text = cmp_ctheatmap_info,
        ##caption = cmp_ctheatmap_caption,
        options = cmp_ctheatmap.opts,
        download.fmt = c("pdf","html"),
        pdf.width = 11, pdf.height = 10,
        height = c(fullH-50,720), width = c("auto",1100), res=c(80,85)
    )
    
    ##================================================================================
    ## CONNECTIVITY MAP
    ##================================================================================
    
    ##ntop=1000;nb=100
    getNeighbourhoodFoldChangeMatrix <- reactive({

        ngs <- inputData()
        res = pgx.getMetaFoldChangeMatrix(ngs, what="trend.limma")
        ##fc0 = sapply(ngs$gx.meta$meta, function(x) unclass(x$fc)[,"trend.limma"])
        fc0 <- res$fc
        rownames(fc0) <- toupper(gsub(".*:","",rownames(fc0)))

        ## normalize???
        fc0 <- scale(fc0, center=TRUE)

        ## set missing values are set to nearly zero
        fc0[is.na(fc0)] <- 1e-4 * rnorm(sum(is.na(fc0)))

        ngenes = 200
        ngenes = as.integer(input$cmp_topgenes)
        fc0 = head(fc0[order(-apply(fc0,1,sd,na.rm=TRUE)),,drop=FALSE],ngenes)
        
        ## --------- compute correlation distance 
        fc0 <- apply(fc0,2,rank,na.last="keep")  ## rank correlation??
        rho <- cor(fc0, use="pairwise")

        res <- list(fc=fc0, rho=rho)
        return(res)
    })

    cmp_connectivitymap.RENDER <- reactive({

        ngs <- inputData()
        req(ngs)
        
        ## get the fold-changes of selected comparison and neighbourhood
        res <- getNeighbourhoodFoldChangeMatrix()
        fc  <- res$fc
        fc[is.na(fc)] <- 0

        ##validate(need(NCOL(fc)>=2, "warning. need multiple comparisons."))
        if(NCOL(fc)<2) return(NULL)        
        
        add.negative=FALSE
        ##add.negative <- ("add negative" %in% input$fc_cmap_options)
        if(add.negative) {
            comparisons=1
            comparisons = input_cmp_comparisons()
            if(is.null(comparisons)) return(NULL)
            ## add NEGATIVE phenotype node
            kk <- which(colnames(fc) %in% comparisons)
            neg.fc <- -fc[,kk]
            colnames(neg.fc) <- paste0("NEG:",colnames(neg.fc))
            fc <- cbind(fc, neg.fc)
            ##k <- which(rownames(pos)==comparison)
            ##rho <- cor(fc[,], fc[,k])[,1]  ## correlation to QUERY node
        }
        
        ## -------- compute t-SNE
        if(input$cmp_cmapclust=="tsne") {
            require(Rtsne)
            perplexity <- pmax(min(ncol(fc)/5,30),2)
            perplexity
            sfc <- scale(fc)
            if(ncol(sfc)<=6) sfc <- cbind(sfc,sfc,sfc,sfc)
            sfc <- sfc + 1e-2*matrix(rnorm(length(sfc)),nrow(sfc),ncol(sfc))
            pos <- Rtsne( t(sfc), is_distance=FALSE, check_duplicates=FALSE,
                         ## pca = TRUE, partial_pca = TRUE,
                         perplexity=perplexity, num_threads=4)$Y
            pos <- pos[1:ncol(fc),]
            xlab="tSNE-x"
            ylab="tSNE-y"
        } else {
            require(irlba)
            sfc <- scale(fc)
            if(ncol(sfc)<=3) sfc <- cbind(sfc,sfc,sfc)
            pos <- irlba(sfc, nv=2)$v
            pos <- pos[1:ncol(fc),]
            xlab="PC1"
            ylab="PC2"
        }
        rownames(pos) = colnames(fc)
        dim(pos)
        pos = scale(pos) ## scale 
        
        ## prepare plotting
        lab <- rownames(pos)
        ##if(input$fc_cmapshownames) lab <- rownames(pos)
        ##if(!("show names" %in% input$fc_cmap_options)) lab <- rep("",nrow(pos))
        jj <- grep("\\]",rownames(pos))
        grp <- rep("[this_data]",nrow(pos))
        grp[jj] <- sub("].*","]",rownames(pos)[jj])
        table(grp)
        ##rr <- abs(rho)**2
        ##sign<- c("-","+")[1 + 1*(sign(rho)>0)]
        ##sign[k] <- "*"
        ##tt <- paste0(rownames(pos)," (",sign,")")
        tt <- rownames(pos)
        stype <- (1 + 1*(grp=="[this_data]"))
        labels_size <- ifelse( nrow(pos) < 20, 20, 15)
        labels_size <- ifelse( nrow(pos) > 400, 12, labels_size)
        labels_size <- 0.7 * labels_size
        
        require(scatterD3)
        scatterD3( pos[,1], pos[,2], transitions=TRUE,
                  xlab = xlab, ylab = ylab,
                  legend_width=0, tooltip_text=tt,
                  lab=lab, labels_size = labels_size, point_opacity = 0.8,
                  size_var=stype, size_range = c(30,150), size_lab="type",
                  ## size_var=rr, size_range = c(20,200)*2.2,
                  ## size_lab="connectivity score",
                  ## symbol_var=sign, symbol_lab="type",
                  col_var=grp, col_lab="data set")

    })

    cmp_connectivitymap.opts = tagList(
        ##tipify( selectInput(ns('cmp_cmapsets'),"Dataset:", choices=NULL, multiple=TRUE),
        ##"Select datasets to compare with external contrast profiles."),
        tipify(radioButtons(ns('cmp_cmapclust'),"Layout:",c("tsne","pca"),inline=TRUE),
               "Choose the plot layout: t-SNE or PCA"),
        tipify(radioButtons(
            ns('cmp_topgenes'),'Top genes:',c(50,200,1000), inline=TRUE,selected=200),
            "Specify the number of top genes for calculating the distances.")
        ##checkboxGroupInput('fc_cmap_options',NULL,c('show names','add negative'),
        ##                   selected=c('show names'),inline=TRUE),
    )

    cmp_connectivitymap_info = "<b>The Connectivity Map (CMap)</b> shows the similarity of the contrasts profiles as a t-SNE plot. Contrasts that are similar will be clustered close together, contrasts that are different are placed farther away. For comparison with external signatures, users can select multiple public datasets in the settings under 'Dataset'."

    cmp_connectivitymap_caption = "<b>Connectivity Map (CMap).</b> The CMap shows the similarity of the contrasts as a t-SNE plot. Contrasts that are similar will be clustered close together, contrasts that are different are placed farther away."
    
    callModule(
        plotModule,
        "cmp_connectivitymap", label="b",
        func = cmp_connectivitymap.RENDER,
        plotlib="scatterD3",
        options = cmp_connectivitymap.opts,
        pdf.width=8, pdf.height=8,
        height = c(fullH-110), res=90,
        title = "Connectivity map",
        info.text = cmp_connectivitymap_info
        ##caption = cmp_connectivitymap_caption,
    )
    
    ##-------------------------------------------------------
    ##---------- UI LAYOUT ----------------------------------
    ##-------------------------------------------------------
    
    cmp_venn_caption = "<b>(a)</b> <b>Venn diagram</b> showing the number of overlapping genes for multiple contrasts. <b>(b)</b> <b>Cumulative fold-change plot</b> of genes in the selected overlap region."
    
    cmp_scatterPlotMatrix_caption = "<b>(a)</b> <b>Pairs plot.</b> Pairwise scatterplots for two or more differential expression profiles for multiple selected contrasts. Similar profiles will show high correlation with points close to the diagonal.  <b>(b)</b> <b>Venn diagram</b> showing the number of overlapping genes for multiple contrasts. <b>(c)</b> <b>Cumulative fold-change plot</b> of genes in the selected overlap region."

    output$cmp_scatterPlotMatrix_UI <- renderUI({
        fillCol(
            ## id = ns("expr_topgenes"),
            height = fullH,
            flex=c(NA,0.02,1), ##height = 370,
            div(HTML(cmp_scatterPlotMatrix_caption),class="caption"),
            br(),
            fillRow(
                flex=c(1.3,0.12,1), ##height = 370,
                plotWidget(ns("cmp_scatterPlotMatrix")),
                br(),
                fillCol(
                    flex = c(0.9,NA,1),
                    plotWidget(ns("cmp_venndiagram")),
                    cmp_venntable_buttons,
                    plotWidget(ns("cmp_cumFCplot"))
                )
            )
        )
    })
    ##outputOptions(output, "cmp_scatterPlotMatrix_UI", suspendWhenHidden=FALSE) ## important!!!    

    output$cmp_pairsPlot_UI <- renderUI({
        fillCol(
            ## id = ns("expr_topgenes"),
            height = fullH,
            flex=c(NA,0.02,1), ##height = 370,
            div(HTML("caption"),class="caption"),
            br(),
            fillRow(
                flex=c(1.1,0.12,1), ##height = 370,
                plotWidget(ns("cmp_pairsPlot")),
                br(),
                fillCol(
                    flex = c(1,0.1,1),
                    tableWidget(ns("cmp_pairsSelectedTable")),
                    br(),
                    tableWidget(ns("cmp_pairsEnrichmentTable"))
                )
            )
        )
    })
    ##outputOptions(output, "cmp_scatterPlotMatrix_UI", suspendWhenHidden=FALSE) ## important!!!

    cmp_ctClusteringUI_caption = "<b>(a)</b> <b>Constrast heatmap.</b> Similarity of the contrasts visualized as a clustered heatmap. The numeric value in the cells correspond to the Pearson correlation coefficient. Red corresponds to positive correlation and blue to negative correlation. <b>(b)</b> <b>Connectivity Map.</b> The CMap shows the similarity of the contrasts as a t-SNE plot. Contrasts that are similar will be clustered close together, contrasts that are different are placed farther away."
    
    output$cmp_ctClustering_UI <- renderUI({
        fillCol(
            ## id = ns("expr_topgenes"),
            height = fullH,
            flex = c(NA,0.035,1), ##height = 370,
            div(HTML(cmp_ctClusteringUI_caption), class="caption"),
            br(),
            fillRow(
                flex = c(1.2,0.1,1),
                height = fullH - 120,
                plotWidget(ns("cmp_ctheatmap")),
                br(),
                plotWidget(ns("cmp_connectivitymap"))
            )
        )
    })
    ##outputOptions(output, "cmp_ctClustering_UI", suspendWhenHidden=FALSE) ## important!!!
    
} ## end-of-Module 
