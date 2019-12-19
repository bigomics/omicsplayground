ProfilingInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
    )
}

ProfilingUI <- function(id) {
    ns <- NS(id)  ## namespace
    fillRow(
        flex = c(1.6,0.10,1),
        height = 780,
        tabsetPanel(
            tabPanel("Cell type",uiOutput(ns("pr_icp_UI"))),
            tabPanel("Markers",uiOutput(ns("pr_markersplot_UI"))),
            tabPanel("CNV",uiOutput(ns("pr_cnaModule_UI"))),
            tabPanel("iTALK",uiOutput(ns("italk_panel_UI"))),
            tabPanel("Monocle",uiOutput(ns("monocle_panel_UI")))            
        ),
        br(),
        tabsetPanel(
            tabPanel("Phenotypes",uiOutput(ns("pr_phenoModule_UI"))),
            tabPanel("Proportions",uiOutput(ns("pr_crosstabModule_UI"))),
            tabPanel("CytoPlot",uiOutput(ns("pr_cytoModule_UI")))      
        )
    )
}

ProfilingModule <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE
    inputData <- env[["load"]][["inputData"]]
    fullH = 780  ## row height of panel
    tabH = 200  ## row height of panel    
    
    description = "<b>Single-Cell Profiling</b>. Visualize the distribution of (inferred)
immune cell types, expressed genes and pathway activation."
    output$description <- renderUI(HTML(description))

    pr_infotext =
        "The <strong>Cell Profiling Module</strong> infers the type of cells using computational deconvolution methods and reference datasets from the literature. Currently, we have implemented a total of 8 methods and 9 reference datasets to predict immune cell types (4 datasets), tissue types (2 datasets), cell lines (2 datasets) and cancer types (1 dataset). However, we plan to expand the collection of methods and databases and to infer other cell types.

<br><br>The <strong>Proportions tab</strong> visualizes the interrelationships between two categorical variables (so-called cross tabulation). Although this feature is very suitable for a single-cell sequencing data, it provides useful information about the proportion of different cell types in samples obtained by the bulk sequencing method.

<br><br>For each combination of gene pairs, the platform can generate a cytometry-like plot of samples under the <strong>Cytoplot</strong> tab. The aim of this feature is to observe the distribution of samples in relation to the selected gene pairs. For instance, when applied to single-cell sequencing data from immunological cells, it can mimic flow cytometry analysis and distinguish T helper cells from the other T cells by selecting the CD4 and CD8 gene combination.

<br><br>The <strong>Markers</strong> section provides potential marker genes, which are the top N=36 genes with the highest standard deviation within the expression data across the samples. For every gene, it produces a t-SNE plot of samples, with samples colored in red when the gene is overexpressed in corresponding samples. Users can also restrict the marker analysis by selecting a particular functional group in which genes are divided into 89 groups, such as chemokines, transcription factors, genes involved in immune checkpoint inhibition, and so on.

<br><br>It is also possible to perform a copy number variation analysis under the <strong>CNV tab</strong>. The copy number is estimated from gene expression data by computing a moving average of the relative expression along the chromosomes. CNV generates a heatmap of samples versus chromosomes, where samples can be annotated further with a phenotype class provided in the data."

    output$inputsUI <- renderUI({
        ui <- tagList(
            ##actionLink("pr_info", "More details ...")
            tipify(actionLink(ns("pr_info"), "Info", icon=icon("info-circle")),
                   "Show more information about this module."),
            hr(),br(),
            tipify(selectInput(ns("pr_samplefilter"),"Filter samples:", choices=NULL, multiple=TRUE),
                   "Filter relevant samples (cells).", placement="top", options = list(container = "body")),
            tipify( actionLink(ns("pr_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Toggle options", placement="top"),
            br(),br(),
            conditionalPanel(
                "input.pr_options % 2 == 1", ns=ns, 
                tagList(
                    tipify(radioButtons(ns('pr_clustmethod'),NULL,c("tsne","pca"), inline=TRUE, selected="tsne"),
                           "Specify a layout for the figures: t-SNE or PCA-based layout.",
                           options = list(container = "body")),
                    tipify(checkboxInput(ns("pr_group"), "group", FALSE),
                           "Group/ungroup the samples (cells).", options = list(container = "body"))
                )
            )
        )
        ui
    })
    outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================

    observeEvent(input$pr_info, {
        showModal(modalDialog(
            title = HTML("<strong>Cell Profiling Module</strong>"),
            HTML(pr_infotext),
            easyClose = TRUE, size="l" ))
    })

    ## update filter choices upon change of data set 
    observe({
        ngs <- inputData()
        req(ngs)
        ## levels for sample filter
        levels <- getLevels(ngs$Y)
        updateSelectInput(session, "pr_samplefilter", choices=levels)
    })

    observe({
        ngs <- inputData()
        req(ngs)
        refsets = "LM22"
        refsets <- sort(names(ngs$deconv))
        refsel <- grep("LM22",refsets,value=TRUE)
        updateSelectInput(session,"pr_refset",choices=refsets, selected=refsel)
        ## updateSelectInput(session,"pr_refset2",choices=refsets, selected=refsel)
        
        dcmethods <- names(ngs$deconv[[1]])
        dcsel <- intersect(c("meta.prod","meta"),dcmethods)[1]
        updateSelectInput(session, "pr_dcmethod", choices=dcmethods, selected=dcsel)

    })

    ##================================================================================
    ##========================= REACTIVE FUNCTIONS ===================================
    ##================================================================================

    pfGetClusterPositions <- reactive({
        ngs <- inputData()
        req(ngs)
        
        ##zx <- filtered_matrix1()
        zx = ngs$X
        kk <- selectSamplesFromSelectedLevels(ngs$Y, input$pr_samplefilter)
        if(length(kk)==0) return(NULL)
        
        zx <- zx[,kk,drop=FALSE]
        zx = head(zx[order(-apply(zx,1,sd)),],1000)
        zx = t(scale(t(zx)))  ## scale??

        if(FALSE && input$pr_group) {
            grp = ngs$Y[colnames(zx),"group"]
            zx = t(apply( zx, 1, function(x) tapply(x, grp, mean)))
        }

        pos = NULL
        if(input$pr_clustmethod=="tsne") {
            require(Rtsne)
            tsne.dim = 3
            ##do3d <- grepl("3D", as.character(input$pr_clustmethod))
            do3d <- ("3D" %in% input$pca.options)
            tsne.dim = c(2,3)[ 1 + 1*do3d]
            force.compute = FALSE
            ## force.compute = TRUE        
            if(!force.compute && tsne.dim==2 && !is.null(ngs$tsne2d) ) {
                pos <- ngs$tsne2d[colnames(zx),]
            } else if(!force.compute && tsne.dim==3 && !is.null(ngs$tsne3d) ) {
                pos <- ngs$tsne3d[colnames(zx),]
            } else {
                perplexity = max(min(30,ncol(zx)/4),1)
                perplexity = min((ncol(zx)-1)/3, 30)
                perplexity
                pos <- Rtsne( t(zx), dim=tsne.dim, check_duplicated=FALSE,
                             num_threads=100, ##Y_init=Y_init, 
                             perplexity=perplexity )$Y
            }
        } else {
            ##cat("pfGetClusterPositions:: computing PCA/SVD...\n")
            require(irlba)
            ##pos <- cmdscale(dist(t(zx)), k=3)
            ##pos = svd(zx,nv=3)$v
            pos = irlba(zx,nv=3)$v        
        }
        
        pos = scale(pos) ## scale 
        colnames(pos) = paste0("dim",1:ncol(pos))
        rownames(pos) = colnames(zx)    

        ##if(input$pca.gx=="<cluster>")
        require(igraph)
        dist = as.dist(dist(pos))
        gr = graph_from_adjacency_matrix(
            1.0/dist, diag=FALSE, mode="undirected")
        clust <- cluster_louvain(gr)$membership

        ##cat("pfGetClusterPositions:: done!\n")
        return( list(pos=pos, clust=clust) )
    })

    ##================================================================================
    ## Cell type
    ##================================================================================

    getDeconvResults <- reactive({
        ngs <- inputData()
        req(ngs)

        method = "meta"
        method <- input$pr_dcmethod
        if(is.null(method)) return(NULL)
        
        refset = "LM22"
        refset <- input$pr_refset
        if(!("deconv" %in% names(ngs))) return(NULL)
        results <- ngs$deconv[[refset]][[method]]
        ## threshold everything (because DCQ can be negative!!!)
        results <- pmax(results,0)
        
        return(results)
    })

    pr_icpplotFUNC <- reactive({
        require(RColorBrewer)
        
        ngs <- inputData()
        req(ngs)
        
        clust <- pfGetClusterPositions()
        if(is.null(clust)) return(NULL)
        pos <- ngs$tsne2d
        pos <- clust$pos
        
        score <- ngs$deconv[["LM22"]][["meta"]]
        score = getDeconvResults()
        if(is.null(score) || length(score)==0  ) return(NULL)
        
        ## normalize
        score <- score[rownames(pos),,drop=FALSE]
        score[is.na(score)] <- 0
        score <- pmax(score,0)
        ##score <- score - min(score,na.rm=TRUE) + 0.01 ## subtract background??
        ##score <- score / (1e-20 + sqrt(rowMeans(score**2,na.rm=TRUE)))
        score <- score / (1e-20 + rowSums(score))
        score <- tanh(score/mean(abs(score)))
        score <- score / max(score,na.rm=TRUE)
        summary(as.vector(score))

        ## take top10 features
        jj.top <- unique(as.vector(apply(score,1,function(x) head(order(-x),10))))
        score <- score[,jj.top]
        score <- score[,order(-colMeans(score**2))]    
        score <- score[,1:min(50,ncol(score))]
        ii <- hclust(dist(score))$order
        jj <- hclust(dist(t(score)))$order
        score <- score[ii,jj]
        score0 <- score
        pos <- pos[rownames(score),]
        
        if(input$pr_group && input$pr_view!="distribution") {
            grp <- ngs$samples[rownames(score),"group"]
            pos <- apply(pos,2,function(x) tapply(x,grp,median))
            score <- apply(score,2,function(x) tapply(x,grp,mean))
            ii <- hclust(dist(score))$order
            jj <- hclust(dist(t(score)))$order
            score <- score[ii,jj]        
        }    
        b0 <- 1 + 0.85*pmax(30 - ncol(score), 0)

        if(input$pr_view=="dotmap") {
            require(corrplot)
            ##gx.heatmap(score)
            par(mfrow=c(1,1), mar=c(0,0,8,1), oma=c(1,1,1,1)*0.5 )
            score3 <- score**1.5
            rownames(score3) <- paste("",rownames(score3),"  ")
            tl.srt=90
            tl.cex=ifelse(nrow(score)>60,0.7,0.85)
            if(max(sapply(rownames(score3),nchar))>30) tl.srt=45
            corrplot( t(score3), mar=c(b0,1,4,1),
                     cl.lim = c(0,max(score3)), cl.pos = "n",
                     tl.cex = tl.cex, tl.col = "grey20",
                     tl.srt = tl.srt )
        }

        if(input$pr_view=="heatmap") {

            ##if(PRO.VERSION) {
            usermode = USERMODE()
            if(!is.null(usermode) && usermode >= 'PRO') {            
                kk <- head(colnames(score)[order(-colMeans(score**2))],18)
                kk <- intersect(colnames(score),kk)
                all.scores <- ngs$deconv[["LM22"]]
                all.scores <- ngs$deconv[[input$pr_refset]]
                if(input$pr_group && input$pr_view!="distribution") {
                    grp <- ngs$samples[rownames(all.scores[[1]]),"group"]
                    for(i in 1:length(all.scores)) {
                        all.scores[[i]] <- apply(all.scores[[i]],2,
                                                 function(x) tapply(x,grp,mean))
                        ii <- rownames(score)
                        all.scores[[i]] <- all.scores[[i]][ii,kk]
                    }
                }    

                nm <- length(all.scores)
                m=3;n=2
                if(nm>6) {m=3;n=3}
                if(nm>9) {m=4;n=3}
                rr <- 2+max(nchar(colnames(score)))/2
                par(mfrow=c(m,n), mar=c(0,0.3,2,0.3), oma=c(10,0,0,rr), xpd=TRUE)
                k=1
                for(k in 1:length(all.scores)) {
                    ii <- rownames(score)
                    score1 <- all.scores[[k]][ii,kk]
                    ##score1 <- score1[rownames(score0),kk]
                    if(k%%n!=0) colnames(score1) <- rep("",ncol(score1))
                    if((k-1)%/%n!=(nm-1)%/%n) rownames(score1) <- rep("",nrow(score1))
                    score1 <- score1 / (1e-8+rowSums(score1))
                    if(nrow(score1) > 100)  rownames(score1) <- rep("",nrow(score1))
                    gx.imagemap( t(score1**1), cex=0.85, main="", clust=FALSE)
                    title(main=names(all.scores)[k], cex.main=1.1, line=0.4, font.main=1)
                }
                
            } else {
                score1 <- score
                score1 <- score1 / (1e-8+rowSums(score1))            
                if(nrow(score1) > 100)  rownames(score1) <- rep("",nrow(score))
                gx.heatmap( t(score1**2), scale="none", 
                           cexRow=1, cexCol=0.6, col=heat.colors(16),
                           mar=c(b0,15), key=FALSE, keysize=0.5)
            }
            
        }    
        
        if(input$pr_view=="distribution") {
            
            cex1 = 1.2
            cex1 <- 0.9*c(2.2,1.1,0.6,0.3)[cut(nrow(pos),breaks=c(-1,40,200,1000,1e10))]
            klrpal = colorRampPalette(c("grey90", "grey50", "red3"))(16)
            ##klrpal = paste0(col2hex(klrpal),"AA")    
            klrpal = paste0(col2hex(klrpal),"66")
            
            par(mfrow=c(5,5), mar=c(0.2,0.2,1.8,0.2), oma=c(1,1,1,1)*0.8 )
            par(mfrow=c(5,5), mar=c(0,0.2,0.5,0.2), oma=c(1,1,6,1)*0.5)
            if(ncol(score)>25) par(mfrow=c(6,6), mar=c(0,0.2,0.5,0.2)*0.6)
            i=1    
            
            jj <- head(order(-colMeans(score**2)),36)
            jj <- jj[order(colnames(score)[jj])]
            colnames(score)[jj]
            for(j in jj) {
                gx = pmax(score[,j],0)
                gx = 1+round(15*gx/(1e-8+max(score)))
                klr0 = klrpal[gx]
                ii <- order(gx)
                ## ii <- sample(nrow(pos))
                plot( pos[ii,], pch=19, cex=1*cex1, col=klr0[ii],
                     xlim=1.2*range(pos[,1]), ylim=1.2*range(pos[,2]),
                     fg = gray(0.8), bty = "o", xaxt='n', yaxt='n',
                     xlab="", ylab="")
                legend( "topleft", legend=colnames(score)[j], bg="#AAAAAA88",
                       cex=1.2, text.font=1, y.intersp=0.8, bty="n",
                       inset=c(-0.05,-0.0) )
            }

            refset <- input$pr_refset
            mtext(refset, outer=TRUE, line=0.5, cex=1.0)
            
        }
    })

    pr_icp.opts = tagList(
        ##radioButtons('pr_view',NULL,c("distribution","dotmap","heatmap"), inline=TRUE),
        tipify(selectInput(ns("pr_view"),"plot type:", choices=c("distribution","dotmap","heatmap")),
               "Specify the plot type: distribution, dotmap, or heatmap.",
               placement="top", options = list(container = "body")),
        tipify(selectInput(ns("pr_refset"), "reference:", choices=NULL),
               "Select a reference dataset for the cell type prediction.",
               placement="top", options = list(container = "body")),
        tipify(selectInput(ns("pr_dcmethod"),"method:", choices=NULL),
               "Choose a method for the cell type prediction.",
               placement="top", options = list(container = "body"))
    )

    pr_icp_info = "<strong>Cell type profiling</strong> infers the type of cells using computational deconvolution methods and reference datasets from the literature. Currently, we have implemented a total of 8 methods and 9 reference datasets to predict immune cell types (4 datasets), tissue types (2 datasets), cell lines (2 datasets) and cancer types (1 dataset). However, we plan to expand the collection of methods and databases and to infer other cell types."

    pr_icp_captionSAVE = "> **Cell type profiling (deconvolution).** The cell type can be computationally inferred using so-called deconvolution methods by matching the (absolute) expression profile to a known reference set. This reference set can be a cell type reference database but also cancer types, tissue types or cell lines. Reference set: `r renderText(input$pr_refset)`. Deconvolution method: `r renderText(input$pr_dcmethod)`"

    pr_icp_caption = "<b>Cell type profiling (deconvolution).</b> The cell type can be computationally inferred using so-called deconvolution methods by matching the (absolute) expression profile to a known reference set. This reference set can be a cell type reference database but also cancer types, tissue types or cell lines."
    
    pr_icp_module <- plotModule(
        "pr_icpplot", ns=ns, func=pr_icpplotFUNC,
        title = "Cell type profiling (deconvolution)",
        options = pr_icp.opts,
        pdf.width=8, pdf.height=8, res=85,
        info.text = pr_icp_info,
        caption = pr_icp_caption
    )
    output <- attachModule(output, pr_icp_module)

    output$pr_icp_UI <- renderUI({
        fillCol(
            height = fullH,
            moduleWidget(pr_icp_module, ns=ns)
        )
    })
    outputOptions(output, "pr_icp_UI", suspendWhenHidden=FALSE) ## important!!!
    
    ##================================================================================
    ## Markers
    ##================================================================================
    
    ##output$pr_markersplot <- renderPlot({
    pr_markersplotFUNC <- reactive({
        ##if(!input$tsne.all) return(NULL)
        require(RColorBrewer)
        
        ngs <- inputData()
        req(ngs)
        
        clust <- pfGetClusterPositions()
        if(is.null(clust)) return(NULL)
        pos <- ngs$tsne2d
        pos <- clust$pos
        
        ##markers <- ngs$families[["CD family"]]
        if(is.null(input$pr_features)) return(NULL)
        if(input$pr_features=="") return(NULL)

        term = ""
        if(input$pr_level=="gene") {
            markers <- ngs$families[["Transcription factors (ChEA)"]]
            if(input$pr_search!="") {
                term = input$pr_search
                jj <- grep(term, ngs$genes$gene_name, ignore.case=TRUE )
                markers <- ngs$genes$gene_name[jj]
                term = paste("filter:",term)
            } else if(input$pr_features %in% names(ngs$families)) {
                markers <- ngs$families[[input$pr_features]]
                term = input$pr_features
            } else {
                markers <- ngs$genes$gene_name
            }
            ##markers <- intersect(markers, rownames(ngs$X))       
            markers <- intersect(markers,ngs$genes$gene_name)
            jj <- match(markers,ngs$genes$gene_name)
            pmarkers <- intersect(rownames(ngs$genes)[jj],rownames(ngs$X))
            gx <- ngs$X[pmarkers,rownames(pos),drop=FALSE]
        } else if(input$pr_level=="geneset") {
            ##markers <- ngs$families[["Immune checkpoint (custom)"]]
            markers <- COLLECTIONS[[1]]
            if(is.null(input$pr_features)) return(NULL)
            ft <- input$pr_features
            if(input$pr_search=="" && ft %in% names(COLLECTIONS)) {
                markers <- COLLECTIONS[[input$pr_features]]
                markers <- intersect(markers, rownames(ngs$gsetX))
                term = input$pr_features
            } else if(input$pr_search!="") {
                term = input$pr_search
                jj <- grep(term, rownames(ngs$gsetX), ignore.case=TRUE )
                markers <- rownames(ngs$gsetX)[jj]
                term = paste("filter:",term)
            } else {
                markers <- rownames(ngs$gsetX)
            }
            gx <- ngs$gsetX[markers,rownames(pos),drop=FALSE]
        } else {
            cat("fatal error")
            return(NULL)
        }
        
        ## prioritize gene with large variance (groupwise)
        grp <- as.character(ngs$samples[rownames(pos),"group"])
        zx <- t(apply(gx,1,function(x) tapply(x,grp,mean)))
        gx <- gx[order(-apply(zx,1,sd)),,drop=FALSE]
        gx <- gx - min(gx,na.rm=TRUE) + 0.01 ## subtract background??    
        rownames(gx) = sub(".*:","",rownames(gx))
        
        ##gx <- tanh(gx/sd(gx) ) ## softmax
        cex1 = 1.0
        cex1 <- 0.8*c(2.2,1.1,0.6,0.3)[cut(nrow(pos),breaks=c(-1,40,200,1000,1e10))]
        klrpal <- colorRampPalette(c("grey90", "grey80", "grey70", "grey60","red4", "red3"))(16)
        klrpal = colorRampPalette(c("grey90", "grey60", "red3"))(16)
        klrpal = paste0(col2hex(klrpal),"66")

        NP=25
        if(input$pr_level=="gene") NP=36
        top.gx = head(gx,NP)  ## match number of plot below!
        top.gx = top.gx[order(rownames(top.gx)),,drop=FALSE]
        top.gx = pmax(top.gx,0)
        ##top.gx <- tanh(top.gx/mean(top.gx))

        plevel="gene"
        plevel <- input$pr_level
        
        par(mfrow=c(1,1)*sqrt(NP), mar=c(0,0.2,0.5,0.2)*0.6, oma=c(1,1,1,1)*0.5)
        par(mfrow=c(1,1)*sqrt(NP), mar=c(0,0.2,0.5,0.2)*0.6, oma=c(1,1,6,1)*0.5)
        ##par(mfrow=c(6,6), mar=c(0,0.2,0.5,0.2), oma=c(1,1,1,1)*0.5)
        i=1    
        for(i in 1:min(NP,nrow(top.gx))) {

            colvar = pmax(top.gx[i,],0) 
            colvar = 1+round(15*(colvar/(0.7*max(colvar)+0.3*max(top.gx))))
            klr0 = klrpal[colvar]

            ii <- order(colvar)
            ##ii <- sample(nrow(pos))
            plot( pos[ii,], pch=19, cex=cex1, col=klr0[ii],
                 xlim=1.1*range(pos[,1]), ylim=1.1*range(pos[,2]),
                 fg = gray(0.8), bty = "o",
                 xaxt='n', yaxt='n', xlab="tSNE1", ylab="tSNE2")

            if(plevel=="gene") {
                gene <- sub(".*:","",rownames(top.gx)[i])
                ##title( gene, cex.main=1.0, line=0.3, col="grey40", font.main=1)
                legend( "topleft", legend=gene, bg="#AAAAAA88",
                       cex=1.2, text.font=1, y.intersp=0.8,
                       bty="n", inset=c(-0.05,-0.0) )
            } else {
                gset <- sub(".*:","",rownames(top.gx)[i])
                gset1 <- breakstring(substring(gset,1,80),24,force=TRUE)
                gset1 <- tolower(gset1) 
                ##title( gset1, cex.main=0.9, line=0.4, col="grey40", font.main=1)
                legend( "topleft", legend=gset1, cex=0.95, bg="#AAAAAA88",
                       text.font=2, y.intersp=0.8, bty="n",
                       inset=c(-0.05,-0.0) )
            }
        }
        mtext(term, outer=TRUE, cex=1.0, line=0.6)
        
    })

    pr_markersplot.opts = tagList(
        tipify(selectInput(ns("pr_level"),"Level:", choices=c("gene","geneset")),
               "Specify the level of the marker analysis: gene or gene set level.",
               placement="top", options = list(container = "body")),
        tipify(selectInput(ns("pr_features"),"Feature set:", choices=NULL, multiple=FALSE),
               "Select a particular functional group for the analysis.",
               placement="top", options = list(container = "body")),
        tipify(textInput(ns("pr_search"),"Filter:"),
               "Filter markers by a specific keywords.",
               placement="top", options = list(container = "body"))
    )

    pr_markersplot_info = "The Markers section produces for the top marker genes, a t-SNE with samples colored in red when the gene is overexpressed in corresponding samples. The top genes (N=36) with the highest standard deviation are plotted. <p>In the plotting options, users can also restrict the marker analysis by selecting a particular functional group in which genes are divided into 89 groups, such as chemokines, transcription factors, genes involved in immune checkpoint inhibition, and so on."

    pr_markersplot_caption = "<b>T-SNE distribution of expression of marker genes.</b> Good biomarkers will show a distribution pattern strongly correlated with some phenotype. The top genes with the highest standard deviation are shown. The red color shading is proportional to the (absolute) expression of the gene in corresponding samples." 

    pr_markersplot_module <- plotModule(
        "pr_markersplot", ns=ns, func=pr_markersplotFUNC,
        title = "Distribution of marker genes",
        pdf.height = 10, pdf.width=10, res=85,
        options = pr_markersplot.opts,
        info.text = pr_markersplot_info,
        caption = pr_markersplot_caption
    )
    output <- attachModule(output, pr_markersplot_module)

    observe({
        ngs <- inputData()
        req(ngs,input$pr_level)    
        choices <- names(ngs$families)
        selected = grep("^CD",choices,ignore.case=TRUE,value=TRUE)[1]
        if(input$pr_level=="geneset") {
            choices <- names(COLLECTIONS)
            selected = grep("HALLMARK",names(COLLECTIONS),ignore.case=TRUE,value=TRUE)
        }
        updateSelectInput(session, "pr_features", choices=choices, selected=selected)

    })
    
    output$pr_markersplot_UI <- renderUI({
        fillCol(
            height = fullH,
            moduleWidget(pr_markersplot_module, ns=ns)
        )
    })

    ##================================================================================
    ## CNV
    ##================================================================================

    getCNAfromExpression <- reactive({
        ngs <- inputData()
        req(ngs)
        
        ##source("../R/pgx-cna.R");source("../R/gx-heatmap.r")
        cat("<profiling:getCNAfromExpression> calculating CNV using heuristic smoothing...\n")
        withProgress( message='calculating CNV...', value=0.66, {
            res <- pgx.CNAfromExpression(ngs, nsmooth=40, downsample=1)
        })
        return(res)
    })

    getCNAfromExpression.inferCNV <- reactive({
        ngs <- inputData()
        req(ngs)
        
        ##source("../R/pgx-cna.R");source("../R/gx-heatmap.r")
        cat("<profiling:getCNAfromExpression> calculating CNV using inferCNV...\n")
        withProgress( message='calculating inferCNV...', value=0.66, {
            res <- pgx.inferCNV(ngs, refgroup=NULL)
        })
        return(res)
    })

    ##output$pr_cnaplot <- renderPlot({
    pr_cnaplotFUNC <- reactive({
        require(RColorBrewer)
        ##return(NULL)    
        ngs <- inputData()
        req(ngs,input$pr_cna_method,input$pr_cna_annotvar,input$pr_cna_orderby)
        
        if(input$pr_cna_method=="inferCNV") {
            res <- getCNAfromExpression.inferCNV()
            par(mfrow=c(1,1))
            grid::grid.raster(res$png)
        } else {
            res <- getCNAfromExpression()
            annotvar=NA
            annotvar <- input$pr_cna_annotvar
            if(annotvar=="<none>") annotvar <- NULL
            order.by <- input$pr_cna_orderby
            pgx.plotCNAHeatmap(ngs, res, annot=annotvar, order.by=order.by)
        }
        
    })

    pr_cna.opts = tagList(
        tipify(radioButtons(ns("pr_cna_method"),label="Method:", choices=c("sma40","inferCNV"),inline=TRUE),
               "Select the computational method for CNV inference. The <tt>sma40</tt> method uses a fast moving average of the relative expression values with a window of 40 genes. <tt>inferCNV</tt> uses the method inferCNV of the Trinity CTAT Project (warning: this method is very slow!)."),
        tipify(selectInput(ns("pr_cna_annotvar"),label="Annotate with:", choices=NULL, multiple=FALSE),
               "Select what annotation variable to show together with the heatmap", placement = "top"),
        ##checkboxGroupInput('pr_cnaoptions','',c("bin20"), inline=TRUE),
        ##radioButtons('pr_cnaplottype',NULL,c("image","heatmap","splitmap"), inline=TRUE),
        tipify(radioButtons(ns("pr_cna_orderby"),"Order samples by:",c("clust","pc1"), inline=TRUE),
               "Select how to order the vertical (sample) axis: clustering or according the loading of the first principal component.")
    )

    pr_cnaModule_info = "<strong>Copy number variation (CNV)</strong> analysis. The copy number is estimated from gene expression data by computing a moving average of the relative expression along the chromosomes. CNV generates a heatmap of samples versus chromosomes, where samples can be annotated further with a phenotype class provided in the data."

    pr_cnaModule_caption = "<b>Inferred copy number variation (CNV)</b>. Genomic copy number is estimated from gene expression data by computing a moving average of the relative expression along the chromosomes. The heatmap shows the inferred copy number of the samples versus chromosomes. The samples are annotated further with phenotype information on the right side of the figure."

    pr_cnaModule <- plotModule(
        "pr_cnaplot", ns=ns, func=pr_cnaplotFUNC,
        options = pr_cna.opts,
        title = "Inferred copy number variation (CNV)",
        pdf.width=10, pdf.height=8, res=110,
        info.text = pr_cnaModule_info,
        caption = pr_cnaModule_caption
    )
    output <- attachModule(output, pr_cnaModule)

    observe({
        ngs <- inputData()
        ##if(is.null(ngs)) return(NULL)
        req(ngs)        
        ## levels for sample filter
        annotvar <- c(colnames(ngs$Y),"<none>")
        updateSelectInput(session, "pr_cna_annotvar", choices=annotvar)
        
    })

    output$pr_cnaModule_UI <- renderUI({
        fillCol(
            height = fullH,
            moduleWidget(pr_cnaModule, ns=ns)
        )
    })

    ##================================================================================
    ## iTALK
    ##================================================================================
    
    italk_getResults <- reactive({
        ngs <- inputData()
        ## if(is.null(ngs)) return(NULL)
        req(ngs)
        
        require(iTALK)
        db <- iTALK::database
        db.genes <- unique(c(db$Ligand.ApprovedSymbol,db$Receptor.ApprovedSymbol))
        length(db.genes)
        ##genes <- intersect(genes, rownames(ngs$X))
        pp <- rownames(ngs$X)
        db.genes <- intersect(db.genes, toupper(ngs$genes[pp,"gene_name"]))

        ## make groups
        ph <- "cell.type"
        ph <- input$italk_groups
        ct <- as.character(ngs$samples[,ph])    
        table(ct)
        
        ##data <- data.frame(cell_type=ct, t(log2(1 + ngs$counts[genes,])))
        pp1 <- rownames(ngs$genes)[match(db.genes, ngs$genes$gene_name)]
        gx <- t(ngs$X[pp1,])
        colnames(gx) <- db.genes
        gx0 <- apply(gx,2,function(x) tapply(x,ct,mean))

        top_genes <- 50
        top_genes <- input$italk_netview_topgenes
        
        data1 <- data.frame(cell_type=ct, gx)
        dim(data1)
        
        ## find the ligand-receptor pairs from highly expressed genes
        cell_col <- rep(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a'),99)
        ct.names = unique(as.character(data1$cell_type))
        cell_col <- cell_col[1:length(ct.names)]
        names(cell_col) <- ct.names
        
        comm_type='cytokine'
        comm_type <- input$italk_category
        
        mode = "absolute"
        ##mode <- input$italk_mode
        if(mode=="absolute") {
            highly_exprs_genes <- rawParse(data1, top_genes=50, stats='mean')    
            res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
            res_cat <- res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
        } else {
            ## contrast <- input$fa_contrast
            ## group <- ngs$model.parameters$exp.matrix[,contrast]
            ## data1$compare_groups <- group
            ## data1 <- data1[which(data1$compare_groups!=0),]        
            ## ## find DEGenes of regulatory T cells and NK cells between these 2 groups
            ## deg_t  <- DEG(data %>% filter(cell_type=='CD4Tcells'),method='Wilcox',contrast=c(2,1))
            ## deg_nk <- DEG(data %>% filter(cell_type=='NKcells'),method='Wilcox',contrast=c(2,1))
            ## ## find significant ligand-receptor pairs and do the plotting
            ## res_cat <- FindLR(deg_t,deg_nk,datatype='DEG',comm_type=comm_type)
            ## res_cat <- res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC,decreasing=T),]        
        }
        
        is.ligand = (colnames(gx0) %in% db[,"Ligand.ApprovedSymbol"])
        is.receptor = (colnames(gx0) %in% db[,"Receptor.ApprovedSymbol"])
        lr.type = c("L","R","LR")[ 1*is.ligand + 2*is.receptor]
        names(lr.type) <- colnames(gx0)
        table(lr.type)
        res <- list( table=res_cat, exprs=gx0, cell_col=cell_col, lr.type=lr.type)
        return(res)
    })

    italk_netview.RENDER <- reactive({
        res <- italk_getResults()
        req(res)
        ##if(is.null(res)) return(NULL)
        res_cat <- res$table
        ## Communication graph
        NetView(res_cat, col=res$cell_col, vertex.label.cex=1, arrow.width=1, edge.max.width=5)
    })

    italk_LRPlot.RENDER <- reactive({
        ## Circos plot
        res <- italk_getResults()
        req(res)
        ##if(is.null(res)) return(NULL)
        res_cat <- res$table
        if(nrow(res_cat)<1) return(NULL) 

        ntop=25
        ntop = as.integer(input$italk_LRPlot_ntop)
        res_top <- head(res_cat,ntop)
        LRPlot(res_top, datatype='mean count', cell_col=res$cell_col,
               link.arr.lwd = head(res_cat$cell_from_mean_exprs,ntop),
               link.arr.width = head(res_cat$cell_to_mean_exprs,ntop))
        comm_type <- isolate(input$italk_category)
        title((paste(comm_type,"genes     ")), line=0.5)
    })

    italk_heatmap.RENDER <- reactive({
        ## Expression heatmap
        ngs <- inputData()
        res <- italk_getResults()
        req(ngs,res)
        ##if(is.null(res)) return(NULL)
        
        res_cat <- res$table
        ntop=50
        ntop = as.integer(input$italk_LRPlot_ntop)
        res_top <- head(res_cat,ntop)
        genes_top <- sort(unique(c(res_top$ligand,res_top$receptor)))
        if(length(genes_top)==0) return(NULL)
        gx0  <- t(res$exprs[,genes_top,drop=FALSE])
        rownames(gx0) <- paste0(rownames(gx0)," (",res$lr.type[rownames(gx0)],")")

        par(oma=c(3,2,3,0))
        gx.heatmap(gx0, scale="none", mar=c(15,8),
                   cexRow=1, cexCol=1.3,
                   key=FALSE, keysize=0.6)
    })


    require(shinyBS)

    italk_LRPlot_info = "The Ligand-Receptor plot visualizes the communication structure of ligand-receptor genes as
a circle plot. The width of the arrow represents the expression level/log fold change of the ligand; while the width of arrow head represents the expression level/log fold change of the receptor. Different color and the type of the arrow stands for whether the ligand and/or receptor are upregulated or downregulated. For further information, see iTALK R package (Wang et al., BioRxiv 2019)."

    italk_LRPlot_module <- plotModule(
        "italk_LRPlot", ns=ns, func=italk_LRPlot.RENDER,
        title = "Ligand-Receptor plot", label="a",
        info.text = italk_LRPlot_info,
        options = tagList(
            tipify( selectInput(ns("italk_LRPlot_ntop"),"ntop pairs",
                                choices=c(10,15,25,50,75,100),selected=25),
                   "Select the maximum number of LR pairs to include in the LR plot.",
                   placement="top")
        ),
        pdf.width=6, pdf.height=8, res=80
    )
    output <- attachModule(output, italk_LRPlot_module)


    italk_heatmap_info = "The heatmap visualizes the expression level/log fold change of the ligand/receptor genes. For further information, see iTALK R package (Wang et al., BioRxiv 2019)."
    
    italk_heatmap_module <- plotModule(
        "italk_heatmap", ns=ns, func=italk_heatmap.RENDER,
        title = "Expression heatmap", label="b",
        info.text = italk_heatmap_info,
        options = tagList(),
        pdf.width=6, pdf.height=8, res=80
    )
    output <- attachModule(output, italk_heatmap_module)

    
    italk_netview_info = "The NetView plot visualizes the communication structure of ligand-receptor genes as a graph. The colors represent different types of cells as a structure and the width of edges represent the strength of the communication. Labels on the edges show exactly how many interactions exist between two types of cells. For further information, see iTALK R package (Wang et al., BioRxiv 2019)."

    italk_netview_module <- plotModule(
        "italk_netview", ns=ns, func=italk_netview.RENDER,
        title = "NetView", label="c",
        info.text = italk_netview_info,
        options = tagList(
            tipify( selectInput(ns("italk_netview_topgenes"),"top genes",
                                choices=c(10,25,50,75,100),selected=50),
                   "Select the number of topgenes to search for ligand-receptor pairs.",
                   placement="top" )
        ),
        pdf.width=6, pdf.height=8, res=80
    )
    output <- attachModule(output, italk_netview_module)

    observe({
        ngs <- inputData()
        req(ngs)
        ph <- colnames(ngs$samples)
        updateSelectInput(session, "italk_groups", choices=ph)
    })


    italk_panel_caption = "<b>Visualization of ligand-receptor interaction.</b> The iTALK R package is designed to profile and visualize the ligand-receptor mediated intercellular cross-talk signals from single-cell RNA sequencing data (scRNA-seq). iTALK uses a manually curated list of ligand-receptor gene pairs further classified into 4 categories based on the primary function of the ligand: cytokines/chemokines, immune checkpoint genes, growth factors, and others. <b>(a)</b> The Ligand-Receptor plot visualizes the communication structure of ligand-receptor genes as a circle plot. <b>(b)</b> The heatmap visualizes the expression level/log fold change of the ligand/receptor genes." 

    output$italk_panel_UI <- renderUI({
        fillCol(
            height = fullH,
            flex = c(NA,1,0.05,NA), 
            inputPanel(
                tipify( selectInput(ns("italk_groups"),"Group by:", choices=""),
                       "Select the phenotype parameter to divide samples into groups.",
                       placement="right"),
                tipify( selectInput(ns("italk_category"),"Gene category:",
                                    choices=c('cytokine','growth factor','checkpoint','other')),
                       "Select the gene category for finding L/R pairs.",
                       placement="right")
                ##radioButtons("italk_mode","Mode", c("absolute","differential"))
            ),
            fillRow(
                ##flex = c(1,1,1),
                flex = c(1,1),
                moduleWidget(italk_LRPlot_module, ns=ns),
                moduleWidget(italk_heatmap_module, ns=ns)
                ##moduleWidget(italk_netview_module, ns=ns)
            ),
            br(),
            div(HTML(italk_panel_caption),class="caption")
        )
    })
    outputOptions(output, "italk_panel_UI", suspendWhenHidden=FALSE) ## important!!!
    
    
    ##================================================================================
    ## Trajectory (dev)
    ##================================================================================

    monocle_getResults <- reactive({

        require(monocle3)
        
        ngs <- inputData()
        req(ngs)
        ##if(is.null(ngs)) return(NULL)

        dbg("monocle_getResults: reacted")
        
        ## Create a Progress object
        progress <- shiny::Progress$new()
        on.exit(progress$close())    
        progress$set(message = "Calculating trajectories", value = 0)
        
        ##------------------------------------------------------------
        ## Step 0: Make Monocle object from ngs
        ##------------------------------------------------------------
        NGENES=2000
        jj <- head(order(-apply(log(1+ngs$counts),1,sd)),NGENES)
        ngs$counts <- ngs$counts[jj,]
        ngs$genes  <- ngs$genes[rownames(ngs$counts),]
        ngs$genes$gene_short_name <- ngs$genes$gene_name
        ngs$samples$.cluster <- ngs$samples$cluster  ## save to avoid name clash
        ngs$samples$cluster <- NULL
        cds <- monocle3::new_cell_data_set(
                             expression_data = ngs$counts,
                             cell_metadata = ngs$samples,
                             gene_metadata = ngs$genes)

        ##------------------------------------------------------------
        ## Step 1: Normalize and pre-process the data
        ##------------------------------------------------------------
        progress$inc(0.1, detail = "preprocessing")
        ##cds <- monocle3::preprocess_cds(cds, num_dim = 100, residual_model_formula_str = "~ batch")
        cds <- monocle3::preprocess_cds(cds, num_dim = 40)
        ##plot_pc_variance_explained(cds)
        
        ##------------------------------------------------------------
        ## Step 2: Reduce the dimensions using UMAP
        ##------------------------------------------------------------
        progress$inc(0.1, detail = "reduce dimensions")
        cds <- monocle3::reduce_dimension(cds)  ## default is UMAP
        ##cds <- monocle3::reduce_dimension(cds, reduction_method="tSNE")

        ##------------------------------------------------------------
        ## Step 3: Cluster the cells
        ##------------------------------------------------------------
        progress$inc(0.2, detail = "clustering cells")
        cds <- cluster_cells(cds, reduction_method="UMAP")
        
        ##------------------------------------------------------------
        ## Step 4: Learn trajectory graph
        ##------------------------------------------------------------
        progress$inc(0.2, detail = "learning graph")
        cds <- learn_graph(cds)

        ##------------------------------------------------------------
        ## Step 5: Order cells
        ##------------------------------------------------------------
        ##cds <- order_cells(cds)
        ##plot_cells(cds)

        ##------------------------------------------------------------
        ## After: update selectors
        ##------------------------------------------------------------
        GENES = rownames(cds)
        updateSelectInput(session, "monocle_plotgene", choices=GENES)
        grps = setdiff(colnames(cds@colData),c("Size_Factor"))
        updateSelectInput(session, "monocle_groupby", choices=grps,
                          selected=".cluster" )

        progress$inc(0.2, detail = "done")

        dbg("monocle_getResults: DONE!")
        
        return(cds)
    })

    monocle_topMarkers <- reactive({
    })

    monocle_plotTopMarkers.RENDER <- reactive({

        dbg("monocle_plotTopMarkers.RENDER: reacted")
        
        cds <- monocle_getResults()
        require(monocle3)
        req(cds,input$monocle_groupby)

        dbg("monocle_plotTopMarkers.RENDER: 1")
        
        ## Find marker genes expressed by each cluster
        pheno1 = "cluster"
        pheno1 = input$monocle_groupby
        if(pheno1=="cluster") pheno1 <- ".cluster"
        pheno1
        marker_test_res = top_markers(cds, group_cells_by=pheno1, cores=4)

        dbg("monocle_plotTopMarkers.RENDER: 2")
        
        NTOP = 1
        NTOP = 3
        NTOP = as.integer(input$monocle_ntop)
        top_specific_markers = marker_test_res %>%
            filter(fraction_expressing >= 0.10) %>%
            group_by(cell_group) %>%
            top_n(NTOP, specificity)
        ##top_n(1, marker_test_q_value)
        ##top_n(1, pseudo_R2)

        dbg("monocle_plotTopMarkers.RENDER: 3")
        
        top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))    
        scale_max1 = 0.8 * log(max(assay(cds))+0.1)
        ##par(mar=c(0,0,0,0))
        g <- monocle3::plot_genes_by_group(cds,
                                           top_specific_marker_ids,
                                           group_cells_by = pheno1,
                                           ##ordering_type="maximal_on_diag",
                                           scale_max = scale_max1, max.size = 5)
        ## g <- qplot(seq(0,4*pi,0.1), sin(seq(0,4*pi,0.1)))
        g <- g + theme(plot.margin=unit(c(1,1,1,1)*0.5,"cm"))

        dbg("monocle_plotTopMarkers.RENDER: DONE!")

        return(g)
    })

    monocle_plotTrajectory.RENDER <- reactive({

        dbg("monocle_plotTrajectory.RENDER: reacted")

        require(monocle3)
        cds <- monocle_getResults()
        req(cds,input$monocle_groupby)    

        pheno1 = ".cluster"
        pheno1 = input$monocle_groupby
        if(pheno1=="cluster") pheno1 <- ".cluster"    

        if(pheno1=="pseudotime") {
            root_cells <- input$monocle_rootcells
            cds = order_cells(cds, root_cells=root_cells)
        }
        
        size1 <- ifelse( ncol(cds) < 100, 2, 1)

        g <- monocle3::plot_cells(cds,
                                  cell_size = size1,
                                  group_label_size = 5,
                                  color_cells_by = pheno1,
                                  label_groups_by_cluster=FALSE,
                                  label_leaves=FALSE,
                                  label_branch_points=FALSE)
        g <- g + theme(plot.margin=unit(c(1,1,1,1)*0.5,"cm"))
        return(g)
    })

    monocle_plotGene.RENDER <- reactive({

        dbg("monocle_plotGene.RENDER: reacted")
        
        cds <- monocle_getResults()
        req(cds, input$monocle_plotgene)   
        require(monocle3)
        genes = "CD14"
        genes = input$monocle_plotgene
        
        size1 <- ifelse( ncol(cds) < 100, 2, 1)
        g <- monocle3::plot_cells(cds, genes = genes[1],
                                  cell_size = size1, 
                                  show_trajectory_graph=FALSE,
                                  label_cell_groups=FALSE)
        g <- g + theme(plot.margin=unit(c(1,1,1,1)*0.5,"cm"))
        return(g)
    })


    ##------- plotTopMarkers module -------
    monocle_plotTopMarkers.MODULE <- plotModule(
        "monocle_plotTopMarkers", ns=ns,
        func = monocle_plotTopMarkers.RENDER, 
        title = "Top markers by group", label="a", 
        info.text = "The heatmap visualizes the expression of group-specific markers.",
        options = tagList(
            selectInput(ns("monocle_groupby"),"Group by:",choices=NULL),
            selectInput(ns("monocle_ntop"),"ntop:",choices=c(1,2,3,4,5,10,25),selected=5)
        ),
        pdf.width=5, pdf.height=8, res=72, plotlib="ggplot"
    )
    output <- attachModule(output, monocle_plotTopMarkers.MODULE)


    ##------- plotTrajectory module -------
    monocle_plotTrajectory.opts = tagList()
    monocle_plotTrajectory.MODULE <- plotModule(
        "monocle_plotTrajectory", ns=ns,
        func = monocle_plotTrajectory.RENDER,
        title = "Single-cell trajectory", label="b",
        info.text = "Single-cell trajectory analysis how cells choose between one of several possible end states. Reconstruction algorithms can robustly reveal branching trajectories, along with the genes that cells use to navigate these decisions.",
        options = monocle_plotTrajectory.opts,
        pdf.width=8, pdf.height=8, res=80, plotlib="ggplot"
    )
    output <- attachModule(output, monocle_plotTrajectory.MODULE)

    ##------- plotGene module -------
    monocle_plotGene.opts = tagList(    
        selectInput(ns("monocle_plotgene"),"Gene:",choices=NULL)
    )
    monocle_plotGene.MODULE <- plotModule(
        "monocle_plotGene", ns=ns, func=monocle_plotGene.RENDER,
        title = "Gene expression", label="c",
        info.text = ".",
        options = monocle_plotGene.opts,
        pdf.width=8, pdf.height=8, res=80, plotlib="ggplot"
    )
    output <- attachModule(output, monocle_plotGene.MODULE)

    monocle_panel_caption = "<b>Single-cell trajectory analysis</b>.  <b>(a)</b> Heatmap visualizing the expression of the group-specific top markers. <b>(b)</b> Single-cell trajectory analysis how cells choose between one of several possible end states. Reconstruction algorithms can robustly reveal branching trajectories, along with the genes that cells use to navigate these decisions. <b>(c)</b> Gene expression distribution for selected marker gene."

    output$monocle_panel_UI <- renderUI({
        fillCol(
            height = fullH,
            flex = c(1,0.05,NA), 
            ##------- Page layout -------
            fillRow(
                flex = c(1.2,1),
                moduleWidget(monocle_plotTopMarkers.MODULE, outputFunc="plotOutput", ns=ns),
                fillCol(
                    flex = c(1,1),
                    moduleWidget(monocle_plotTrajectory.MODULE, outputFunc="plotOutput", ns=ns),
                    moduleWidget(monocle_plotGene.MODULE, outputFunc="plotOutput", ns=ns)
                )
            ),
            br(),
            div(HTML(monocle_panel_caption), class="caption")
        )
    })
    outputOptions(output, "monocle_panel_UI", suspendWhenHidden=FALSE) ## important!!!
    
    ##================================================================================
    ## Phenotypes
    ##================================================================================

    ##output$pr_phenoplot <- renderPlot({
    pr_phenoplotFUNC <- reactive({
        ##if(!input$tsne.all) return(NULL)
        require(RColorBrewer)
        ngs <- inputData()
        ##if(is.null(ngs)) return(NULL)
        req(ngs)
        
        clust <- pfGetClusterPositions()
        if(is.null(clust)) return(NULL)
        
        pos <- ngs$tsne2d
        pos <- clust$pos
        sel <- rownames(pos)
        pheno = colnames(ngs$Y)

        ## layout
        par(mfrow = c(3,2), mar=c(0.3,0.7,2.8,0.7))
        if(length(pheno)>6) par(mfrow = c(4,3), mar=c(0.3,0.4,2.8,0.4)*0.8)
        if(length(pheno)>12) par(mfrow = c(5,4), mar=c(0.2,0.2,2.5,0.2)*0.8)
        
        cex1 <- 1.6*c(1.8,1.3,0.8,0.5)[cut(nrow(pos),breaks=c(-1,40,200,1000,1e10))]    
        cex1 = cex1 * ifelse(length(pheno)>6, 0.8, 1)
        cex1 = cex1 * ifelse(length(pheno)>12, 0.8, 1)

        ## is it a float/number???
        is.num <- function(y, fmin=0.1) {
            suppressWarnings(numy <- as.numeric(as.character(y)))
            t1 <- !all(is.na(numy)) && is.numeric(numy)
            t2 <- (length(unique(y))/length(y)) > fmin
            (t1 && t2)
        }
        
        i=6    
        for(i in 1:min(20,length(pheno))) {

            px=4
            px=pheno[i]
            y = ngs$Y[sel,px]
            y[which(y %in% c(NA,""," ","NA","na"))] <- NA
            if(sum(!is.na(y))==0) next

            if(is.num(y)) {
                klrpal = colorRampPalette(c("grey90", "grey50", "red3"))(16)
                y = rank(as.numeric(y))
                ny <- round(1 + 15*(y - min(y))/(max(y)-min(y)))
                klr0 = klrpal[ny]
            } else {
                y = factor(as.character(y))
                klrpal = COLORS
                klrpal <- paste0(col2hex(klrpal),"99")
                klr0 = klrpal[y]
            }

            jj = which(is.na(klr0))
            if(length(jj)) klr0[jj] <- "#AAAAAA22"
            plot( pos, pch=19, cex=cex1, col=klr0,fg = gray(0.5), bty = "o",
                 xaxt='n', yaxt='n', xlab="tSNE1", ylab="tSNE2")

            if(!is.num(y)) {
                if(input$pr_labelmode=="legend") {
                    legend("bottomright", legend=levels(y), fill=klrpal,
                           cex=0.9, y.intersp=0.8, bg="white")
                } else {
                    grp.pos <- apply(pos,2,function(x) tapply(x,y,mean,na.rm=TRUE))
                    grp.pos <- apply(pos,2,function(x) tapply(x,y,median,na.rm=TRUE))
                    nvar <- length(setdiff(y,NA))
                    if(nvar==1) {
                        grp.pos <- matrix(grp.pos,nrow=1)
                        rownames(grp.pos) <- setdiff(y,NA)[1]
                    }
                    
                    labels = rownames(grp.pos)
                    ## title("\u2591\u2592\u2593") 
                    boxes = sapply(nchar(labels),function(n) paste(rep("\u2588",n),collapse=""))
                    boxes = sapply(nchar(labels),function(n) paste(rep("█",n),collapse=""))
                    ##boxes = sapply(nchar(labels),function(n) paste(rep("#",n),collapse=""))
                    cex2 <- c(1.4,1.2,1,0.8)[cut(length(labels),breaks=c(-1,5,10,20,999))]    
                    text( grp.pos, labels=boxes, cex=cex2, col="#CCCCCC99")
                    text( grp.pos, labels=labels, font=2, cex=1.1*cex2, col="white")
                    text( grp.pos, labels=labels, font=2, cex=cex2)
                    ##text( grp.pos[,], labels=rownames(grp.pos), font=2, cex=cex1**0.5)
                }
            }
            title(tolower(pheno[i]), cex.main=1.3, line=0.5, col="grey40")
        }
    })


    pr_phenoplot.opts <- tagList(
        tipify( radioButtons(ns('pr_labelmode'),'Label:',c("groups","legend"), inline=TRUE),
               "Select whether you want the group labels to be plotted inside the plots or in a seperate legend.")
    )

    pr_phenoModule_info = "<strong>Phenotypes.</strong> This figure visualizes the distribution of the phenotype information. In the plot options, you can select whether you want the group labels to be plotted inside the plots or in a seperate legend."

    pr_phenoModule_caption = "<b>Phenotype plots.</b> The plots show the distribution of the phenotypes superposed on the t-SNE clustering. Often, we can expect the t-SNE distribution to be driven by the particular phenotype that is controlled by the experimental condition or unwanted batch effects."
   
    pr_phenoModule <- plotModule(
        "pr_phenoplot", ns=ns, func=pr_phenoplotFUNC,
        options = pr_phenoplot.opts,
        pdf.width=5, pdf.height=8, res=90,
        info.text = pr_phenoModule_info,
        caption = pr_phenoModule_caption
    )
    output <- attachModule(output, pr_phenoModule)


    output$pr_phenoModule_UI <- renderUI({    
        fillCol(
            height = fullH,
            flex = c(1), 
            moduleWidget(pr_phenoModule, ns=ns)
        )
    })
    outputOptions(output, "pr_phenoModule_UI", suspendWhenHidden=FALSE) ## important!!!
    
    ##================================================================================
    ## Proportions
    ##================================================================================

    ##output$statsplot <- renderPlot({
    pr_crosstabPlotFUNC <- reactive({
        ##if(!input$tsne.all) return(NULL)
        require(RColorBrewer)
        ngs <- inputData()

        req(ngs)
        req(input$pr_crosstabpheno, input$pr_crosstabvar, input$pr_crosstabgene)
        
        scores = ngs$deconv[[1]][[1]]  ## just an example...
        if(input$pr_crosstabvar == "<cell type>") {
            scores <- getDeconvResults()
            if(is.null(scores)) return(NULL)
            scores <- pmax(scores,0) ## ??
        } else {
            x <- as.character(ngs$Y[,1])
            x <- as.character(ngs$Y[,input$pr_crosstabvar])
            x[is.na(x)] <- "_"
            scores <- model.matrix( ~ 0 + x )
            rownames(scores) <- rownames(ngs$Y)
            colnames(scores) <- sub("^x","",colnames(scores))
        }
        dim(scores)

        ## restrict to selected sample set
        kk <- head(1:nrow(scores),1000)
        kk <- 1:nrow(scores)
        kk <- selectSamplesFromSelectedLevels(ngs$Y, input$pr_samplefilter)
        scores <- scores[kk,,drop=FALSE]
        scores <- scores[,which(colSums(scores)>0),drop=FALSE]
        scores[which(is.na(scores))] <- 0    
        dim(scores)

        ## expected counts per stat level
        ##kk.counts <- colSums(ngs$counts[,kk,drop=FALSE])  ## total count of selected samples
        kk.counts <- colSums(2**ngs$X[,kk,drop=FALSE])  ## approximate counts from log2X
        grp.counts <- ( t(scores / rowSums(scores)) %*% matrix(kk.counts,ncol=1))[,1]  
        
        getProportionsTable <- function(pheno, is.gene=FALSE) {    
            y <- NULL
            ##if("gene" %in% input$pr_crosstaboptions) {
            if( is.gene ) {
                xgene <- ngs$genes[rownames(ngs$X),]$gene_name
                gx <- ngs$X[which(xgene == pheno),kk,drop=FALSE]
                gx.highTH <- mean(gx, na.rm=TRUE)
                y <- paste(pheno,c("low","high"))[ 1 + 1*(gx >= gx.highTH)]
                table(y)
            } else if(pheno %in% colnames(ngs$samples)) {
                y <- ngs$samples[kk,1]
                y <- ngs$samples[kk,pheno]
                pheno <- tolower(pheno)
            } else if(pheno == "<cell type>") {
                res1 <- getDeconvResults()
                res1 <- pmax(res1,0) ## ??
                res1 <- res1[kk,,drop=FALSE]
                ##res1 <- res1[,which(colSums(res1)>0),drop=FALSE]
                y <- colnames(res1)[max.col(res1)]  ## take maximum col??
                remove(res1)
                pheno <- "<cell type>"
            } else {
                return(NULL)
            }   
            
            ## calculate proportions by group
            grp <- factor(as.character(y))
            ngrp <- length(levels(grp))
            grp.score <- apply(scores, 2, function(x) tapply(x,grp,mean,na.rm=TRUE))
            ngrp
            if(ngrp==1) {
                grp.score <- matrix(grp.score,nrow=1)
                rownames(grp.score) <- y[1]
                colnames(grp.score) <- colnames(scores)
            }
            
            ## weighted counts
            grp.score[is.na(grp.score)] <- 0
            grps <- levels(grp)
            grp.score
            fy <-  (table(y) / sum(!is.na(y)))
            jj <- match(rownames(grp.score),names(fy))
            grp.score <-  grp.score * as.vector(fy[jj])
            ## normalize to total 100% 
            grp.score <- grp.score / (1e-20+sum(grp.score))
            
            ## reduce to maximum number of items (x-axis)
            if(0) {
                ##jj <- which(colSums(grp.score) > 0.001)
                jj <- order(-colSums(grp.score))
                j1 <- head(jj, 25)  ## define maximum number of items
                j0 <- setdiff(jj, j1)
                grp.score0 <- grp.score[,j1,drop=FALSE]
                grp.counts0 <- grp.counts[j1]
                if(length(j0)>0) {
                    grp.score0 <- cbind( grp.score0, "other"=rowSums(grp.score[,j0,drop=FALSE]))
                    grp.counts0 <- c(grp.counts0, "other"=sum(grp.counts[j0]))
                }
                grp.score <- grp.score0
                grp.counts <- grp.counts0
                grp.score <- t( t(grp.score) / (1e-20+colSums(grp.score)))  ##
            }
            
            ## normalize to total 100% and reduce to maximum number of items (y-axis)
            jj <- order(-rowSums(grp.score))
            j1 <- head(jj, 10)  ## define maximum number of items
            j0 <- setdiff(jj, j1)
            grp.score0 <- grp.score[j1,,drop=FALSE]
            if(length(j0)>0) {
                grp.score0 <- rbind( grp.score0, "other"=colSums(grp.score[j0,,drop=FALSE]))
            }
            grp.score <- grp.score0
            grp.score <- t( t(grp.score) / (1e-20+colSums(grp.score)))  ##

            ## cluster columns??
            ##dist1 <- dist(t(scale(grp.score)))
            dist1 <- dist(t(grp.score))
            dist1[is.na(dist1)] <- mean(dist1,na.rm=TRUE)
            jj <- hclust(dist1)$order
            grp.score <- grp.score[,jj,drop=FALSE]
            return(grp.score)
        }

        ## select phenotype variable
        head(ngs$samples)
        pheno=1
        pheno="cluster"
        pheno="activated"
        pheno="cell.type"
        pheno="<cell type>"
        pheno <- input$pr_crosstabpheno
        if(is.null(pheno)) return(NULL)

        ##pheno="cluster"
        grp.score1 <- getProportionsTable(pheno, is.gene=FALSE)    
        grp.score2 <- NULL
        gene = ngs$genes$gene_name[1]
        gene = input$pr_crosstabgene
        if(gene != "<none>") {
            grp.score2 <- getProportionsTable(pheno=gene, is.gene=TRUE)
            kk <- colnames(grp.score2)[order(grp.score2[1,])]
            grp.score2  <- grp.score2[,match(kk,colnames(grp.score2))]
            grp.score1  <- grp.score1[,match(kk,colnames(grp.score1))]
        }
        
        jj <- match(colnames(grp.score1),names(grp.counts))
        grp.counts <- grp.counts[jj] / 1e6  ## per million
        names(grp.counts) = colnames(grp.score1)

        ##-------------- plot by estimated cell.type ----------------------
        
        ##par(mar = c(4,6,2,3))
        layout(matrix(c(1,2,3), 3,1), heights=c(2,4,3))
        if(!is.null(grp.score2)) {
            layout(matrix(c(1,2,3,4), 4,1), heights=c(2.2,1,4,2))
        }
        
        ## top bar with counts
        par(mar = c(0,5,5.3,3), mgp=c(2.0,0.8,0) )
        xlim <- c(0,1.2*length(grp.counts))  ## reserves space for legend
        barplot( grp.counts, col="grey50", ylab="counts (M)", cex.axis=0.8,
                cex.lab=0.8, names.arg=NA, xpd=NA, xlim=1.3*xlim, ## log="y", 
                ylim=c(0.01, max(grp.counts)), yaxs="i" )
        ## title(pheno, cex.main=1.2, line=2, col="grey40")
        
        ## middle plot (gene)
        if(!is.null(grp.score2)) {
            klrpal2 = COLORS[1:nrow(grp.score2)]
            klrpal2 = rev(grey.colors(nrow(grp.score2),start=0.45))
            par(mar = c(0,5,0.3,3), mgp=c(2.4,0.9,0) )
            barplot( 100*grp.score2, col=klrpal2, las=3, srt=45,
                    xlim=1.3*xlim, ylim=c(0,99.99),
                    names.arg = rep(NA, ncol(grp.score2)),
                    ylab="(%)", cex.axis=0.90)
            legend(1.02*xlim[2], 100, legend=rownames(grp.score2),
                   fill=klrpal2, xpd=TRUE, cex=0.8, y.intersp=0.8,
                   bg="white", bty="n")
        }

        if(1) {
            ## main proportion graph
            klrpal1 = COLORS[1:nrow(grp.score1)]
            par(mar = c(4,5,0.3,3), mgp=c(2.4,0.9,0) )
            barplot( 100*grp.score1, col=klrpal1, las=3, srt=45, xlim=1.3*xlim,
                    ylim=c(0,99.99), ylab="proportion (%)", cex.axis=0.90)
            legend(1.02*xlim[2], 100, legend=rev(rownames(grp.score1)),
                   fill=rev(klrpal1), xpd=TRUE, cex=0.8, y.intersp=0.8,
                   bg="white", bty="n")
        }
        
    })


    pr_crosstab.opts <- tagList(
        tipify(selectInput(ns("pr_crosstabvar"),label="x-axis:", choices=NULL, multiple=FALSE),
                           "Choose a predefined phenotype group on the x-axis.",
                           placement="top", options = list(container = "body")),
        tipify(selectInput(ns("pr_crosstabpheno"),label="y-axis:", choices=NULL, multiple=FALSE),
               "Choose a predefined phenotype group on the y-axis.",
               placement="top", options = list(container = "body")),
        tipify(selectInput(ns("pr_crosstabgene"),label="gene:", choices=NULL, multiple=FALSE),
               "Visualize the expression barplot of a gene by specifying the gene name.",
               placement="top", options = list(container = "body"))
        ##checkboxGroupInput('pr_crosstaboptions','',c("gene"), inline=TRUE, width='50px')
        ##selectInput("pr_crosstabgene",label=NULL, choices=NULL, multiple=FALSE),
        ##br(), cellArgs=list(width='80px')
    )


    pr_crosstabModule_info = "The <strong>Proportions tab</strong> visualizes the interrelationships between two categorical variables (so-called cross tabulation). Although this feature is very suitable for a single-cell sequencing data, it provides useful information about the proportion of different cell types in samples obtained by the bulk sequencing method."
    
    pr_crosstabModule_caption = "<b>Proportion plot.</b> Plot visualizing the overlap between two categorical variables (so-called cross tabulation). Although this feature is very suitable for a single-cell sequencing data, it provides useful information about the proportion of different cell types in samples obtained by the bulk sequencing method."
   
    pr_crosstabModule <- plotModule(
        "pr_crosstabPlot", ns=ns, func=pr_crosstabPlotFUNC,
        options = pr_crosstab.opts,
        pdf.width=6, pdf.height=8, res=110,
        info.text = pr_crosstabModule_info,
        caption = pr_crosstabModule_caption
    )
    output$pr_crosstabPlot <- pr_crosstabModule$render
    output$pr_crosstabPlot_pdf <- pr_crosstabModule$pdf

    observe({
        ngs <- inputData()
        ##if(is.null(ngs)) return(NULL)
        req(ngs)

        ##if(is.null(input$pr_crosstaboptions)) return(NULL)
        pheno0 <- grep("group|sample|donor|id|batch",colnames(ngs$samples),invert=TRUE,value=TRUE)
        pheno0 <- grep("sample|donor|id|batch",colnames(ngs$samples),invert=TRUE,value=TRUE)
        kk <- selectSamplesFromSelectedLevels(ngs$Y, input$pr_samplefilter)
        nphenolevel <- apply(ngs$samples[kk,pheno0,drop=FALSE],2,function(v) length(unique(v)))
        pheno0 = pheno0[which(nphenolevel>1)]
        genes <- sort(as.character(ngs$genes$gene_name))
        pheno1 <- c("<cell type>",pheno0)
        genes1 <- c("<none>",genes)
        updateSelectInput(session, "pr_crosstabvar", choices=pheno1)
        updateSelectInput(session, "pr_crosstabpheno", choices=pheno1)
        updateSelectInput(session, "pr_crosstabgene", choices=genes1)

    })


    output$pr_crosstabModule_UI <- renderUI({    
        fillCol(
            height = fullH,
            flex = c(1), 
            moduleWidget(pr_crosstabModule, ns=ns)
        )
    })
    outputOptions(output, "pr_crosstabModule_UI", suspendWhenHidden=FALSE) ## important!!!

    ##================================================================================
    ## CytoPlot
    ##================================================================================

    ##output$pr_cytoplot <- renderPlot({
    pr_cytoplotFUNC <- reactive({
        ##if(!input$tsne.all) return(NULL)
        require(RColorBrewer)
        ngs <- inputData()
        ##if(is.null(ngs)) return(NULL)
        req(ngs)
        
        if(is.null(input$pr_cytovar1)) return(NULL)
        if(is.null(input$pr_cytovar2)) return(NULL)
        if(input$pr_cytovar1=="") return(NULL)
        if(input$pr_cytovar2=="") return(NULL)

        kk <- selectSamplesFromSelectedLevels(ngs$Y, input$pr_samplefilter)    
        gene1 <- input$pr_cytovar1
        gene2 <- input$pr_cytovar2
        ##if(gene1 == gene2) return(NULL)
        par(mfrow=c(1,1), mar=c(12,4,6,2))
        pgx.cytoPlot( ngs, gene1, gene2, samples=kk, cex=0.8,
                     col="grey60", cex.names=1)
        
    })

    pr_cyto.opts = tagList(
        tipify(selectInput(ns("pr_cytovar1"),label="x-axis:", choices=NULL, multiple=FALSE),
               "Select your prefered gene on the x-axis.",
               placement="top", options = list(container = "body")),
        tipify(selectInput(ns("pr_cytovar2"),label="y-axis:", choices=NULL, multiple=FALSE),
               "Choose your prefered gene on the y-axis.",
               placement="top", options = list(container = "body"))
    )

    pr_cytoModule_info = "For each combination of gene pairs, the platform can generate a cytometry-like plot of samples under the Cytoplot tab. The aim of this feature is to observe the distribution of samples in relation to the selected gene pairs. For instance, when applied to single-cell sequencing data from immunological cells, it can mimic flow cytometry analysis and distinguish T helper cells from the other T cells by selecting the CD4 and CD8 gene combination."
    
    pr_cytoModule_caption = "<b>Cyto plot.</b> This plot shows the distribution of samples in relation to the expression of selected gene pairs. It mimics the scatter plots used for gating in flow cytometry analysis."    
    
    pr_cytoModule <- plotModule(
        "pr_cytoplot", ns=ns, func=pr_cytoplotFUNC,
        options = pr_cyto.opts,
        pdf.width=6, pdf.height=8, res=80,
        info.text = pr_cytoModule_info,
        caption = pr_cytoModule_caption
    )
    output$pr_cytoplot <- pr_cytoModule$render
    output$pr_cytoplot_pdf <- pr_cytoModule$pdf

    observe({
        ngs <- inputData()
        ##if(is.null(ngs)) return(NULL)
        req(ngs)

        xgenes <- ngs$genes[rownames(ngs$X),]$gene_name
        genes <- sort(as.character(xgenes))
        g1 <- grep("^CD4|^CD8",genes,value=TRUE)[1]
        g2 <- grep("^CD79|^CD3[DEG]|^CD37",genes,value=TRUE)[1]
        if(length(g1)==0) g1 <- genes[1]
        if(length(g2)==0) g2 <- genes[2]
        updateSelectInput(session, "pr_cytovar1", choices=genes, selected=g1)
        updateSelectInput(session, "pr_cytovar2", choices=genes, selected=g2)
    })
    

    output$pr_cytoModule_UI <- renderUI({    
        fillCol(
            height = fullH,
            flex = c(1), 
            moduleWidget(pr_cytoModule, ns=ns)
        )
    })
    outputOptions(output, "pr_cytoModule_UI", suspendWhenHidden=FALSE) ## important!!!
    
    return(NULL)
}
