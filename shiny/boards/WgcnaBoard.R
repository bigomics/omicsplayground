##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing WgcnaBoard")

WgcnaInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

WgcnaUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    ui <- shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("WGCNA",uiOutput(ns("wgcna_UI"))),
            shiny::tabPanel("Modules",uiOutput(ns("modules_UI"))),
            shiny::tabPanel("Eigengenes",uiOutput(ns("eigen_UI"))),
            shiny::tabPanel("Intramodular",uiOutput(ns("intra_UI")))
        )
    )
    ui
}

WgcnaBoard <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE

    inputData <- env[["load"]][["inputData"]]

    fullH = 700  ## full height of page
    rowH1 = 250  ## row 1 height
    rowH2 = 440  ## row 2 height
    
    description = "<b>WGCNA Analysis.</b> Weighted correlation network analysis (WGCNA) is a gene-level cluster analysis method based on pairwise correlations between genes. It allows one to define modules (clusters), intramodular hubs, and network nodes with regard to module membership, to study the relationships between co-expression modules."
    output$description <- shiny::renderUI(shiny::HTML(description))

    infotext ="Weighted gene co-expression network analysis (WGCNA) is a systems biology method for describing the correlation patterns among genes across microarray samples. Weighted correlation network analysis can be used for finding clusters (modules) of highly correlated genes, for summarizing such clusters using the module eigengene or an intramodular hub gene, for relating modules to one another and to external sample traits (using eigengene network methodology), and for calculating module membership measures. Correlation networks facilitate network based gene screening methods that can be used to identify candidate biomarkers or therapeutic targets."
    

    ##================================================================================
    ##========================= OUTPUT UI ============================================
    ##================================================================================
    
    wgcna_caption = "<b>WGCNA module detection.</b> <b>(a)</b> Modules are detected as branches of the resulting cluster tree using the dynamic branch cutting approach. Genes inside a given module are summarized with the module eigengene. The module eigengene of a given module is defined as the first principal component of the standardized expression profiles. <b>(b)</b> Scale independence and mean connectivity plots to determine the soft threshold <b>(c)</b> Topological overlap matrix visualized as heatmap <b>(d)</b> Dimensionality reduction maps colored by WGCNA module <b>(e)</b> Graph network of WGCNA modules."

    output$wgcna_UI <- shiny::renderUI({
        shiny::fillCol(
            flex = c(NA,0.02,0.75,0.08,1),
            height = fullH,
            shiny::div(shiny::HTML(wgcna_caption), class="caption"),
            shiny::br(),
            shiny::fillRow(
                flex = c(1,1),                
                plotWidget(ns('wgcna_geneDendro')),
                ##plotWidget(ns('wgcna_sampleDendro')),
                plotWidget(ns('wgcna_topologyPlots'))                
            ),
            shiny::br(),
            shiny::fillRow(
                flex = c(1.3,1,0.8),
                plotWidget(ns('wgcna_TOMplot')),
                plotWidget(ns('wgcna_umap')),
                plotWidget(ns('wgcna_moduleGraph'))
            )

        )
        
    })
    shiny::outputOptions(output, "wgcna_UI", suspendWhenHidden=FALSE) ## important!!!

    modules_caption =
"<b>WGCNA functional analysis.</b> <b>(a)</b> Module–trait analysis identifies modules that are significantly associated with the measured clinical traits by quantifying the association as the correlation of the eigengenes with external traits. <b>(b)</b> Partial correlation network of genes most correlated to the eigengene. <b>(c)</b> Module enrichment plot of top most enriched genesets. <b>(d)</b> Table of genes in the selected module. <b>(e)</b> Functional enrichment of the module calculated using Fisher's exact test."
    
    output$modules_UI <- shiny::renderUI({
        shiny::fillCol(
            flex = c(NA, 0.05, 2, 0.07, 1),
            ##flex = c(2,1),
            height = fullH,
            shiny::div(shiny::HTML(modules_caption), class="caption"),
            shiny::br(),
            shiny::fillRow(
                flex = c(1.2, 1.1, 0.9),
                plotWidget(ns('wgcna_moduleTrait')),                
                plotWidget(ns('wgcna_corGraph')),
                plotWidget(ns('wgcna_enrichPlot'))
            ),
            shiny::br(),
            shiny::fillRow(
                flex = c(1, 0.08, 2.3),
                tableWidget(ns("wgcna_geneTable")),
                shiny::br(),
                tableWidget(ns("wgcna_enrichTable"))
            )
        )
        
    })
    shiny::outputOptions(output, "modules_UI", suspendWhenHidden=FALSE) ## important!!!


    eigen_caption =
"<b>WGCNA eigengene analysis.</b> <b>(a)</b> It is often interesting to visualizing the network of eigengenes and study the relationships among the found modules. One can use the eigengenes as represen- tative profiles and quantify module similarity by eigengene correlation. <b>(b)</b> For each module, we also define a quantitative measure of 'module membership' (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes to every module."
    
    output$eigen_UI <- shiny::renderUI({
        shiny::fillCol(
            flex = c(NA,0.04,2),
            height = fullH,
            shiny::div(shiny::HTML(eigen_caption), class="caption"),
            shiny::br(),
            shiny::fillRow(
                flex=c(1,0.06,2.5),
                ## plotWidget(ns('wgcna_eigenHeatmap')),
                plotWidget(ns('wgcna_eigenClustering')),
                shiny::br(),
                plotWidget(ns('wgcna_eigenCorrelation'))
            )
        )
        
    })
    shiny::outputOptions(output, "eigen_UI", suspendWhenHidden=FALSE) ## important!!!

    intra_caption =
        "<b>WGCNA intramodular analysis.</b> We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules."
    
    output$intra_UI <- shiny::renderUI({
        shiny::fillCol(
            flex = c(NA,0.04,2,1),
            height = fullH,
            shiny::div(shiny::HTML(intra_caption), class="caption"),
            shiny::br(),
            shiny::fillRow(
                flex=c(1,0.06,2.5),
                ## plotWidget(ns('wgcna_eigenHeatmap')),
                plotWidget(ns('wgcna_intraHeatmap')),
                shiny::br(),
                plotWidget(ns('wgcna_intraScatter'))
            )
        )
        
    })
    ##outputOptions(output, "intra_UI", suspendWhenHidden=FALSE) ## important!!!
    
    
    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================

    
    output$inputsUI <- shiny::renderUI({
        ui <- shiny::tagList(
            shiny::actionLink(ns("info"), "Info", icon=icon("info-circle")),
            shiny::hr(), shiny::br(),             
            
            ## data set parameters
            shiny::selectInput(ns('selected_module'),'select module', choices=NULL),
            shiny::actionButton(ns("compute"),"Compute!",icon=icon("running"),
                         class="run-button"),            
            shiny::br(),
            shiny::br(),
            shiny::actionLink(ns("options"), "Options", icon=icon("cog", lib = "glyphicon")),
            shiny::br(),br(),br(),
            shiny::conditionalPanel(
                "input.options % 2 == 1", ns=ns,
                shiny::tagList(
                    shiny::selectInput(ns("ngenes"),"Number genes:",
                                choices = c(500,1000,2000,4000,8000),
                                selected = 1000),
                    shiny::selectInput(ns("minmodsize"),"Min. module size",
                                choices = c(10,30,100,250),
                                selected = 30 ),
                    shiny::selectInput(ns("power"),"Power", c(2,4,6,10), selected=6),
                    shiny::selectInput(ns("deepsplit"),"deepsplit", 0:4, selected=2),
                    shiny::selectInput(ns("cutheight"),"Merge cut height",
                                choices = c(0.05, 0.10, 0.25, 0.5, 0.9, 0.999),
                                selected = 0.25)
                )
            )
        )

    })
    shiny::outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!
    
    ##================================================================================
    ##======================= PRECOMPUTE FUNCTION ====================================
    ##================================================================================
    
    ##wgcna.compute <- shiny::reactive({
    wgcna.compute <- shiny::eventReactive( {
        input$compute
        ngs <- inputData()
        1
    }, {

        ngs <- inputData()

        if("wgcna" %in% names(ngs)) {
            message("[wgcna.compute] >>> using pre-computed WGCNA results...")
            return( ngs$wgcna )
        }
                
        pgx.showSmallModal("Calculating WGCNA...<br>please wait")        
        progress <- shiny::Progress$new()
        on.exit(progress$close())    
        progress$set(message = "Calculating WGCNA...", value = 0)
        message("[wgcna.compute] >>> calculating WGCNA...")
        if(0) {
            shinyalert::shinyalert(
                title = "",
                text = "No WGCNA data found in PGX object. Computing now.. "
            )
        }
        

        WGCNA::enableWGCNAThreads()
        
        if(0) {
            liv <- read.csv("~/Downloads/LiverFemale3600.csv")
            X <- as.matrix(liv[,9:ncol(liv)])
            rownames(X) <- liv$gene_symbol
            X <- X[order(-apply(X,1,sd)),]
            X <- X[!duplicated(rownames(X)),]
            dim(X)
        }
        
        X <- as.matrix(ngs$X)
        dim(X)
        X <- X[order(-apply(X,1,sd,na.rm=TRUE)),]
        X <- X[!duplicated(rownames(X)),]
        
        minmodsize=30;power=6;cutheight=0.25;deepsplit=2;ngenes=1000
        
        ngenes <- input$ngenes
        minmodsize <- as.integer(input$minmodsize)
        power      <- as.numeric(input$power)
        cutheight  <- as.numeric(input$cutheight)
        deepsplit  <- as.integer(input$deepsplit)

        datExpr <- t(head(X,ngenes))
        progress$inc(0.1, "computing WGCNA modules...")            
        net = WGCNA::blockwiseModules(
                         datExpr, power = power,
                         TOMType = "unsigned", minModuleSize = minmodsize,
                         reassignThreshold = 0, mergeCutHeight = cutheight,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         deepSplit = deepsplit,
                         ## saveTOMs = TRUE, saveTOMFileBase = "WCNA.tom",
                         verbose = 3)
        names(net)
        table(net$colors)

        ## clean up traits matrix
        datTraits <- ngs$samples
        ## no dates please...
        isdate <- apply(datTraits, 2, is.Date)
        datTraits <- datTraits[,!isdate,drop=FALSE]

        ## Expand multi-class discrete phenotypes into binary vectors
        ##datTraits1 <- datTraits
        tr.class <- sapply(type.convert(datTraits),class) 
        sel1 <- which(tr.class %in% c("factor"))
        sel2 <- which(tr.class %in% c("integer","numeric"))
        tr1 <- datTraits[,0]
        if(length(sel1)) {
            tr1 <- expandPhenoMatrix(datTraits[,sel1,drop=FALSE], drop.ref=FALSE)
        }
        ## keeping numeric phenotypes
        tr2 <- datTraits[,sel2,drop=FALSE]
        datTraits <- cbind(tr1, tr2)
        
        ## get colors of eigengene modules
        me.genes <- tapply( names(net$colors), net$colors, list)
        names(me.genes) <- paste0("ME",names(me.genes))        
        color1 <- labels2rainbow(net)
        me.colors <- color1[!duplicated(color1)]
        names(me.colors) <- paste0("ME",names(me.colors))
        me.colors <- me.colors[names(me.genes)]        
        progress$inc(0.4,"")
        
        if(1) {            
            message("[wgcna.compute] >>> calculating WGCNA clustering...")
            progress$inc(0.1, "computing dim reductions...")
            
            X1 <- t(datExpr)
            X1 <- t(scale(datExpr))
            ##pos <- Rtsne::Rtsne(X1)$Y

            ##dissTOM <- 1 - abs(cor(datExpr))**6
            dissTOM  <- 1 - WGCNA::TOMsimilarityFromExpr(datExpr, power=power)
            rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)            
            clust <- pgx.clusterBigMatrix(dissTOM, methods=c("umap","tsne","pca"), dims=c(2))
            ##pos <- pgx.clusterBigMatrix(t(X1), methods="tsne", dims=2)[[1]]
            ##pos <- pgx.clusterBigMatrix(dissTOM, methods="pca", dims=2)[[1]]
            names(clust)
            if("cluster.genes" %in% names(ngs)) {
                clust[['umap2d']] <- ngs$cluster.genes$pos[['umap2d']]
            }
            progress$inc(0.2)                    

        }

        if(1) {

            message("[wgcna.compute] >>> calculating WGCNA module enrichments...")
            progress$inc(0,"calculating module enrichment...")

            gmt <- GSETS
            gmt <- GSETS[grep("HALLMARK|GOBP|^C[1-9]",names(GSETS))]
            gse <- NULL
            ##bg <- unlist(me.genes)
            bg <- toupper(rownames(ngs$X))
            i=1
            for(i in 1:length(me.genes)) {
                gg <- toupper(me.genes[[i]])
                rr <- gset.fisher( gg, gmt, background=bg, fdr=1 )
                rr <- cbind( module = names(me.genes)[i],
                            geneset = rownames(rr), rr)
                rr <- rr[order(rr$p.value),,drop=FALSE]
                if(i==1) gse <- rr
                if(i>1) gse <- rbind(gse, rr)
            }
            rownames(gse) <- NULL

            progress$inc(0.3)

        }
        
        ## construct results object
        out <- list(
            datExpr = datExpr,
            datTraits = datTraits,
            net = net,
            gse = gse,
            clust = clust,
            me.genes = me.genes,
            me.colors = me.colors
        )

        shiny::updateSelectInput(session, "selected_module", choices = names(me.genes), sel="ME1" )
        
        message("[wgcna.compute] >>> done!")
        beepr::beep(2)  ## short beep
        shiny::removeModal()
        
        out
    })

    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    shiny::observeEvent( input$info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>WGCNA Analysis Board</strong>"),
            shiny::HTML(infotext),
            easyClose = TRUE ))
    })
    
    
    ##================================================================================
    ##======================= PLOTTING FUNCTIONS =====================================
    ##================================================================================

    ##----------------------------------------
    ##------------ samples dendro ------------
    ##----------------------------------------

    labels2rainbow <- function(net) {
        hc <- net$dendrograms[[1]]
        nc <- length(unique(net$colors))
        n <- length(net$colors)
        ii <- hc$order
        col1 <- labels2colors(net$colors)                
        col.rnk <- rank(tapply(1:n,col1[ii],mean))
        new.col <- rainbow(nc)[col.rnk]
        ## new.col <- heat.colors(nc)[col.rnk]
        names(new.col) <- names(col.rnk)
        new.col["grey"] <- "#AAAAAA"
        new.col
        new.col <- new.col[col1]
        names(new.col) <- net$colors
        new.col
    }

    ##----------------------------------------
    ##--------- topology analysis ------------
    ##----------------------------------------

    wgcna_topologyPlots.RENDER %<a-% shiny::reactive({

        message("[wgcna_topologyPlots.RENDER] reacted")

        out <- wgcna.compute()

        ## Choose a set of soft-thresholding powers
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        ## Call the network topology analysis function
        sft = pickSoftThreshold(out$datExpr, powerVector = powers, verbose = 5)

        ## Plot the results:
        ## sizeGrWindow(9, 5)
        par(mfrow = c(1,2), mar=c(3.3,3,1,1), mgp=c(2,0.8,0))
        cex1 = 0.9;
        ## Scale-free topology fit index as a function of the soft-thresholding power
        KEGGgraph::plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             type="n",
             xlab = "Soft threshold (power)",
             ylab = "SFT model fit (signed R^2)",
             ## main = paste("Scale independence")
             main=NULL
             )

        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             labels = powers, cex=cex1, col="red")
        ## this line corresponds to using an R^2 cut-off of h
        abline(h=0.90,col="red")

        ## Mean connectivity as a function of the soft-thresholding power
        KEGGgraph::plot(sft$fitIndices[,1], sft$fitIndices[,5], type="n",
             xlab="Soft threshold (power)",
             ylab="Mean connectivity", 
             ##main = paste("Mean connectivity")
             main=NULL
             )
        text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    })
    
    wgcna_topologyPlots_opts = shiny::tagList()

    wgcna_topologyPlots_info = "<b>WGCNA topology analysis.</b> Analysis of network topology for various soft-thresholding powers. The left panel shows the scale-free fit index (y-axis) as a function of the soft-thresholding power (x-axis). The right panel displays the mean connectivity (degree, y-axis) as a function of the soft-thresholding power (x-axis)."

    shiny::callModule(
        plotModule, 
        id = "wgcna_topologyPlots", ##ns=ns,
        title = "Scale independence and mean connectivity", label="b",
        func = wgcna_topologyPlots.RENDER,
        func2 = wgcna_topologyPlots.RENDER, 
        download.fmt = c("png","pdf"),
        ## options = wgcna_sampleDendro_opts,
        info.text = wgcna_topologyPlots_info,        
        height = c(rowH1, 600), width = c('auto',1200),
        pdf.width=10, pdf.height=5, res=c(72,100),
        add.watermark = WATERMARK
    )

    ##----------------------------------------
    ##-------------- TOM plot ----------------
    ##----------------------------------------
    
    wgcna_TOMplot.RENDER %<a-% shiny::reactive({
        
        message("[wgcna_TOMplot.RENDER] reacted")
        
        out <- wgcna.compute()
        net <- out$net
        datExpr <- out$datExpr
        geneTree = net$dendrograms[[1]]
        ##moduleColors <- labels2colors(out$net$colors)
        moduleColors <- labels2rainbow(out$net)
        MEs <- out$net$MEs
        
        ## Calculate topological overlap anew: this could be done
        ## more efficiently by saving the TOM calculated during
        ## module detection, but let us do it again here.
        power = 6
        power  <- as.numeric(shiny::isolate(input$power))
        dissTOM = 1 - TOMsimilarityFromExpr(datExpr, power = power)
        rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)
        
        nSelect = 999999
        nSelect = 400
        ## For reproducibility, we set the random seed
        set.seed(10)
        select = head( 1:ncol(dissTOM), nSelect)
        selectTOM = dissTOM[select, select];
        ## There’s no simple way of restricting a clustering tree
        ## to a subset of genes, so we must re-cluster.
        ##selectTree = hclust(as.dist(selectTOM), method = "ward.D2")
        selectTree = hclust(as.dist(selectTOM), method = "average")
        selectColors = moduleColors[select];
        ## Taking the dissimilarity to a power, say 10, makes the plot
        ## more informative by effectively changing the color palette;
        ## setting the diagonal to NA also improves the clarity of the
        ## plot
        plotDiss = selectTOM^7;
        diag(plotDiss) = NA;
        myheatcol = gplots::colorpanel(250,'red','orange','lemonchiffon')
        myheatcol = gplots::colorpanel(250,'lemonchiffon','orange','red')

        par(oma=c(2,0,0,0))
        plotly::layout(matrix(c(0, 0, 5, 0,
                        0, 0, 2, 0,
                        4, 1, 3, 6), nr=3, byrow=T),
               widths  = c(2.3,0.5,10,1.8),
               heights = c(2.3,0.5,10) )

        WGCNA::TOMplot(
                   plotDiss, selectTree, selectColors, col=myheatcol,
                   setLayout = FALSE,
                   main = NULL
                   ##main = paste0("Network heatmap plot (subsampled)")
               )

        if(0) {

            source("~/Playground/omicsplayground/R/gx-heatmap.r")
            D <- selectTOM
            diag(D) <- 0
            ann <- data.frame(col=selectColors)
            rownames(ann) <- colnames(D)

            par(mfrow=c(1,1))
            ii <- selectTree$order
            image(D[ii,ii])

            Heatmap(D[ii,ii],
                    top_annotation = HeatmapAnnotation(module=ann[ii,]),
                    cluster_rows=FALSE,
                    cluster_columns=FALSE)
            
            gx=D;col.annot=ann;symm=TRUE;scale="none"            
            gx.heatmap(D, symm=TRUE, scale="none",
                       dist.method = "euclidean",
                       col.dist.method = "euclidean",
                       clust.method = "average",
                       nmax = 99999,
                       col.annot=ann, verbose=3, annot.ht=2)


            gx.heatmap(plotDiss[ii,ii], clust.method=NULL)
            
        }
        
        ## add color legend
        frame()            
        me.names  <- colnames(MEs)
        me.nr <- as.integer(sub("ME","",me.names))
        ii <- order(me.nr)
        label.colors <- labels2rainbow(net)
        me.colors <- label.colors[!duplicated(names(label.colors))]
        me.colors <- me.colors[as.character(me.nr)]
        
        legend(-0.1,1, legend=me.names[ii], fill=me.colors[ii],
               cex=1.2, bty='n', x.intersp=0.5)
        
    })

    wgcna_TOMplot_opts = shiny::tagList()
    wgcna_TOMplot_info = "<b>WGCNA Topological Overlap Matrix (TOM) heatmap.</b>"    

    shiny::callModule(
        plotModule, 
        id = "wgcna_TOMplot", ##ns=ns,
        title="TOM heatmap", label="c",
        func  = wgcna_TOMplot.RENDER,
        func2 = wgcna_TOMplot.RENDER, 
        download.fmt = c("png","pdf"),
        ## options = wgcna_geneDendro_opts,
        info.text = wgcna_TOMplot_info,        
        height = c(rowH2, 650), width = c('auto',1000),
        pdf.width=10, pdf.height=5, res=c(72,90),
        add.watermark = WATERMARK
    )


    ##----------------------------------------
    ##------------ samples dendro ------------
    ##----------------------------------------

    wgcna_sampleDendro.RENDER %<a-% shiny::reactive({

        message("[wgcna_sampleDendro.RENDER] reacted")

        out <- wgcna.compute()
        datExpr <- out$datExpr
        pheno <- out$datTraits
        
        sampleTree2 = hclust(dist(datExpr), method = "average")
        ipheno <- apply(pheno, 2, function(x) as.numeric(factor(x)))
        colnames(ipheno) <- colnames(pheno)
        rownames(ipheno) <- rownames(pheno)
        traitColors = WGCNA::numbers2colors(ipheno, signed = FALSE)
        ## Plot the sample dendrogram and the colors underneath.
        par(mfrow=c(1,1), mar=c(1,1,1,1)*0)
        WGCNA::plotDendroAndColors(
                   sampleTree2, traitColors[,],
                   groupLabels = colnames(ipheno),
                   cex.colorLabels = 0.8, cex.dendroLabels = 0.9, cex.rowText = 0.8,
                   marAll = c(0.2, 5, 0.2, 0.2),
                   ## main = "Sample dendrogram and trait heatmap"
                   main = NULL
               )

    })
    
    wgcna_sampleDendro_opts = shiny::tagList()
    wgcna_sampleDendro_info = "<b>WGCNA sample dendrogram and trait heatmap.</b>"    

    shiny::callModule(
        plotModule, 
        id = "wgcna_sampleDendro", ##ns=ns,
        title = "Sample dendrogram and trait heatmap", label="b",
        func = wgcna_sampleDendro.RENDER,
        func2 = wgcna_sampleDendro.RENDER, 
        download.fmt = c("png","pdf"),
        ## options = wgcna_sampleDendro_opts,
        info.text = wgcna_sampleDendro_info,        
        height = c(rowH1, 650), width = c('auto',1000),
        pdf.width=10, pdf.height=5, res=c(72,90),
        add.watermark = WATERMARK
    )

    ##----------------------------------------
    ##------------ gene dendro ---------------
    ##----------------------------------------
    
    wgcna_geneDendro.RENDER %<a-% shiny::reactive({

        message("[wgcna_geneDendro.RENDER] reacted")

        out <- wgcna.compute()
        net <- out$net        

        ## Convert labels to colors for plotting
        ##mergedColors = labels2colors(net$colors)
        mergedColors <- labels2rainbow(net)
        ## Plot the dendrogram and the module colors underneath
        plotDendroAndColors(
            dendro = net$dendrograms[[1]],
            colors = mergedColors[net$blockGenes[[1]]],
            ##"Module colors",
            dendroLabels = FALSE, hang = 0.03,
            addGuide = FALSE, guideHang = 0.05,
            marAll = c(0.2, 5, 0.4, 0.2),
            main = NULL
        )
        

    })

    wgcna_geneDendro_opts = shiny::tagList()
    wgcna_geneDendro_info = "<b>WGCNA gene dendrogram and gene modules.</b>"    

    shiny::callModule(
        plotModule, 
        id = "wgcna_geneDendro", ##ns=ns,
        title="Gene dendrogram and gene modules", label="a",
        func  = wgcna_geneDendro.RENDER,
        func2 = wgcna_geneDendro.RENDER, 
        download.fmt = c("png","pdf"),
        ## options = wgcna_geneDendro_opts,
        info.text = wgcna_geneDendro_info,        
        height = c(rowH1, 650), width = c('auto',1000),
        pdf.width=10, pdf.height=5, res=c(72,90),
        add.watermark = WATERMARK
    )

    ##----------------------------------------
    ##------------ module-trait --------------
    ##----------------------------------------
    
    wgcna_moduleTrait.RENDER %<a-% shiny::reactive({
        
        message("[wgcna_moduleTrait.RENDER] reacted")

        out <- wgcna.compute()
        net <- out$net
        datExpr <- out$datExpr
        datTraits <- out$datTraits
        ## moduleLabels <- as.character(out$net$colors)
        ## moduleColors <- labels2colors(out$net$colors)
        moduleColors <- labels2rainbow(out$net)
        MEs <- out$net$MEs
        
        ## Define numbers of genes and samples
        nGenes = ncol(datExpr);
        nSamples = nrow(datExpr);

        if(0) {
            ## Recalculate MEs with color as labels
            MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
            MEs  = orderMEs(MEs0)
        }
        
        moduleTraitCor = cor(MEs, out$datTraits, use = "pairwise");
        moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
        
        textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "")
        textMatrix = signif(moduleTraitCor, 2)

        dim(textMatrix) = dim(moduleTraitCor)
        dim(moduleTraitCor)
        
        sel1 <- 1:nrow(moduleTraitCor)
        sel2 <- 1:ncol(moduleTraitCor)        

        sel2 <- sort(head(order(-colMeans(abs(moduleTraitCor))),40)) ## conditions
        sel1 <- sort(head(order(-rowMeans(abs(moduleTraitCor[,sel2]))),12)) ## eigenvectors
        
        sel2 <- sort(head(order(-colMeans(pmax(moduleTraitCor,0))),40)) ## conditions
        sel1 <- sort(head(order(-rowMeans(pmax(moduleTraitCor[,sel2],0))),12)) ## eigenv        

        message("[wgcna_moduleTrait.RENDER] sel1 = ",paste(sel1,collapse=" "))
        message("[wgcna_moduleTrait.RENDER] sel2 = ",paste(sel2,collapse=" "))

        par(mar = c(3, 12, 1.6, 1.5))
        ## Display the correlation values within a heatmap plot
        labeledHeatmap(Matrix = t(moduleTraitCor[sel1,sel2]),
                       yLabels = colnames(out$datTraits)[sel2],
                       xLabels = colnames(MEs)[sel1],
                       xSymbols = colnames(MEs)[sel1],
                       xLabelsAngle = 90,                               
                       colorLabels = FALSE,
                       colors = greenWhiteRed(50),
                       textMatrix = t(textMatrix[sel1,sel2]),
                       setStdMargins = FALSE,
                       cex.text = 0.7,
                       cex.lab = 0.9,
                       zlim = c(-1,1),
                       ##main = paste("Module-trait relationships")
                       main = NULL
                       )
        
        ## Will display correlations and their p-values
        ##sizeGrWindow(10,6)
        
    })

    wgcna_moduleTrait_opts = shiny::tagList(
        shiny::checkboxInput(ns("traits_binarize"),"binarize continuous vars", FALSE)
    )
    wgcna_moduleTrait_info = "<b>WGCNA module and trait relationship.</b>"    
    
    shiny::callModule(
        plotModule, 
        id = "wgcna_moduleTrait", ##ns=ns,
        title="Module-Trait relationships", label="a",
        func  = wgcna_moduleTrait.RENDER,
        func2 = wgcna_moduleTrait.RENDER, 
        download.fmt = c("png","pdf"),
        options = wgcna_moduleTrait_opts,
        info.text = wgcna_moduleTrait_info, info.width="200px",        
        height = c(420, 650), width = c('auto',1000),
        pdf.width=10, pdf.height=5, res=c(72,90),
        add.watermark = WATERMARK
    )

    ##----------------------------------------
    ##------------ module-graph --------------
    ##----------------------------------------
    
    wgcna_moduleGraph.RENDER %<a-% shiny::reactive({

        message("[wgcna_moduleGraph.RENDER] reacted")

        out <- wgcna.compute()
        net <- out$net
        datExpr <- out$datExpr
        datTraits <- out$datTraits
        ## moduleLabels <- as.character(out$net$colors)
        ## moduleColors <- labels2colors(out$net$colors)
        moduleColors <- labels2rainbow(out$net)
        ii <- which(!duplicated(out$net$colors))
        me.colors <- moduleColors[ii]
        names(me.colors) <- out$net$colors[ii]
        me.colors
        
        ## Recalculate MEs with color as labels
        ##MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
        ##MEs  = orderMEs(MEs0)
        MEs <- out$net$MEs
        dim(MEs)
        
        clust <- hclust(dist(t(MEs)))
        clust
        phylo <- ape::as.phylo(clust)
        gr <- igraph::as.igraph(phylo, directed=FALSE)

        is.tip <- grepl("^ME",V(gr)$name)
        module.nr   <- as.integer(sub("^ME|Node.*","",V(gr)$name))
        module.size <- table(out$net$colors) 
        module.size <- module.size / mean(module.size)
        module.size
        
        igraph::V(gr)$label <- igraph::V(gr)$name
        igraph::V(gr)$label[!is.tip] <- NA
        ##V(gr)$color <- WGCNA::labels2colors(module.nr)
        igraph::V(gr)$color <- me.colors[as.character(module.nr)]
        igraph::V(gr)$size  <- 24 * (module.size[as.character(module.nr)])**0.5
        igraph::V(gr)$size[is.na(igraph::V(gr)$size)] <- 0

        par(mfrow=c(1,1), mar=c(1,1,1,1)*0)
        KEGGgraph::plot(
            gr,
            layout=layout.kamada.kawai,
            ##vertex.color = "lightblue",
            ##vertex.size = 20*is.tip,
            vertex.label.cex = 0.8,
            ## vertex.label.dist = 1.8,
            edge.width = 3
        )

        
    })

    wgcna_moduleGraph_opts = shiny::tagList()
    wgcna_moduleGraph_info = "<b>WGCNA module graph.</b>"    

    shiny::callModule(
        plotModule, 
        id = "wgcna_moduleGraph", ##ns=ns,
        title="Module graph", label="e",
        func  = wgcna_moduleGraph.RENDER,
        func2 = wgcna_moduleGraph.RENDER, 
        download.fmt = c("png","pdf"),
        ## options = wgcna_geneDendro_opts,
        info.text = wgcna_moduleGraph_info,        
        height = c(420, 650), width = c('auto',1000),
        pdf.width=10, pdf.height=5, res=c(72,90),
        add.watermark = WATERMARK
    )

    ##----------------------------------------
    ##--------- Eigengenes heatmap --------------
    ##----------------------------------------
    
    wgcna_eigenHeatmap.RENDER %<a-% shiny::reactive({

        message("[wgcna_eigenHeatmap.RENDER] reacted")

        out <- wgcna.compute()
        net <- out$net
        datExpr <- out$datExpr
        datTraits <- out$datTraits
        
        MEs <- net$MEs
        rho <- cor(MEs, datExpr)
        sdx <- apply(datExpr,2,sd)
        ## rho <- t(t(rho) * sdx**2)
        dim(rho)

        if(input$mask_markers) {
            ME <- as.character(net$colors)
            M <- t(model.matrix( ~ 0 + ME))[rownames(rho),]
            rho <- rho * pmax(M,0.33)
        }
            
        ntop <- floor(200 / ncol(MEs))
        ntop        
        ii <- apply(rho, 1, function(x) head(order(-x),ntop))
        ii <- unique(as.vector(t(ii)))
        ii <- head(ii,70)

        gx.heatmap(t(rho[,ii]), keysize=0.2, mar=c(4,5), key=FALSE,
                   cexRow=0.85, cexCol=1, scale="none" )

        
    })

    wgcna_eigenHeatmap_opts = shiny::tagList(
        shiny::checkboxInput(ns("mask_markers"),"mask_markers", TRUE)
    )
    wgcna_eigenHeatmap_info = "<b>WGCNA Eigengene correlation heatmap.</b> The heatmap shows the correlation of genes to the module eigengenes."    

    shiny::callModule(
        plotModule, 
        id = "wgcna_eigenHeatmap", ##ns=ns,
        title="Eigengene correlation heatmap", label="a",
        func  = wgcna_eigenHeatmap.RENDER,
        func2 = wgcna_eigenHeatmap.RENDER, 
        download.fmt = c("png","pdf"),
        options = wgcna_eigenHeatmap_opts,
        info.text = wgcna_eigenHeatmap_info,        
        height = c(fullH,650), width = c('auto',650),
        pdf.width=6, pdf.height=10, res=c(72,90),
        add.watermark = WATERMARK
    )


    ##----------------------------------------
    ##--------------- gene table -------------
    ##----------------------------------------
    
    wgcna_geneTable.RENDER <- shiny::reactive({

        out <- wgcna.compute()
        
        k <- input$selected_module
        genes <- out$me.genes[[k]]
        tt <- GENE.TITLE[toupper(genes)]
        rho <- cor( out$datExpr[,genes], out$net$MEs[,k])[,1]
        
        df <- data.frame( module=k, gene=genes, me.rho=rho, title = tt )
        numeric.cols <- grep("score|value|ratio|rho",colnames(df))
        
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

    wgcna_geneTable_info = "Genes in the selected WGCNA module."
    wgcna_geneTable_caption = "<b>Module gene table.</b> <b>(a)</b> ..."    
    wgcna_geneTable.opts <- shiny::tagList(
        ##selectInput(ns('wgcna_geneTable_selmodule'),'module:', choices=NULL)
    )
    
    wgcna_geneTable_module <- shiny::callModule(
        tableModule, id = "wgcna_geneTable",
        ## caption = wgcna_geneTable_caption,
        func  = wgcna_geneTable.RENDER, ## ns=ns,
        options = wgcna_geneTable.opts,
        info.text = wgcna_geneTable_info,
        title = "Module genes", label="d",
        height = c(250,650)
    )


    ##----------------------------------------
    ##--------- enrichment table -------------
    ##----------------------------------------
    
    wgcna_enrich_table <- shiny::reactive({
        out <- wgcna.compute()
        df  <- out$gse
        k <- input$selected_module
        message("[wgcna_enrichTable.RENDER] 1 : k = ",k)        
        if(length(k)==0 || k=="") k <- "<all>"
        message("[wgcna_enrichTable.RENDER] 2 : k = ",k)        
        if(k %in% df$module) {
            df <- df[df$module == k,,drop=FALSE ]
        }
        df$score <- df$odd.ratio * -log10(df$p.value)
        df <- df[,c("module","geneset","score","p.value","q.value",
                    "odd.ratio","overlap","genes")]
        df <- df[order(-df$score),]
        df
    })

    wgcna_enrichTable.RENDER <- shiny::reactive({

        df <- wgcna_enrich_table()
        numeric.cols <- grep("score|value|ratio",colnames(df))
        
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

    wgcna_enrichTable_info = "In this table, users can check mean expression values of features across the conditions for the selected genes."
    wgcna_enrichTable_caption = "<b>Module enrichment table.</b> <b>(a)</b> ..."    
    wgcna_enrichTable.opts <- shiny::tagList(
        ##selectInput(ns('wgcna_enrichTable_selmodule'),'module:', choices=NULL)
    )
    
    wgcna_enrichTable_module <- shiny::callModule(
        tableModule, id = "wgcna_enrichTable",
        ## caption = wgcna_enrichTable_caption,
        func  = wgcna_enrichTable.RENDER, ## ns=ns,
        options = wgcna_enrichTable.opts,
        info.text = wgcna_enrichTable_info,
        title = "Module enrichment", label="e",
        height = c(250,650)
    )
    
    ##----------------------------------------
    ##--------- enrichment plot --------------
    ##----------------------------------------
    
    wgcna_enrichPlot.RENDER %<a-% shiny::reactive({

        df <- wgcna_enrich_table()
        if(is.null(df) || nrow(df)==0) return(NULL)        
        ii <- wgcna_enrichTable_module$rows_all()
        shiny::req(ii)
        df <- df[ii,,drop=FALSE]        
        df <- head(df, 20)
        gs.top <- df$geneset
        xlim0 <- c(0, max(df$score))
        col1 <- c("lightskyblue1","lightpink")[1 + 1*(df$q.value<0.05)]
        par(mar=c(4.5,1,1,1))        
        barplot( rev(df$score), horiz=TRUE, width=0.8, space=0.25, xlim=xlim0,
                border=NA, col=rev(col1), xlab="score  (odd.ratio * -log10p)")
        text(0, (nrow(df):1)-0.48, gs.top, adj=0, pos=4, cex=0.8)

    })

    wgcna_enrichPlot_info = "Module enrichment plot."
    wgcna_enrichPlot_caption = "<b>Module enrichment plot.</b> <b>(a)</b> ..."    
    wgcna_enrichPlot.opts <- shiny::tagList(
        ##selectInput(ns('wgcna_enrichPlot_selmodule'),'module:', choices=NULL)
    )
    
    
    shiny::callModule(
        plotModule,
        id = "wgcna_enrichPlot", label = "c",
        func = wgcna_enrichPlot.RENDER,
        func2 = wgcna_enrichPlot.RENDER,        
        info.text = wgcna_enrichPlot_info,
        options = wgcna_enrichPlot.opts,
        title = "Enrichment plot",
        ## caption = topEnriched_caption
        caption2 = wgcna_enrichPlot_info,
        ##pdf.width = 14, pdf.height = 4, 
        height = c(420,650),
        width = c('auto',800),
        res = c(72,80),
        add.watermark = WATERMARK
    )
    
    ##-----------------------------------------------------------
    ## Correlation network
    ##-----------------------------------------------------------
    
    wgcna_corGraph.RENDER %<a-% shiny::reactive({

        out <- wgcna.compute()

        k = "ME1"
        k <- input$selected_module
        shiny::req(k)
        message("[wgcna_corGraph.RENDER] k = ",k)
        genes <- out$me.genes[[k]]
        
        dim(out$datExpr)
        xx <- cbind( out$net$MEs[,k,drop=FALSE], out$datExpr[,genes])
        rho1 <- cor(xx, out$net$MEs[,k] )[,1]
        ntop <- min(nrow(xx)-1,20)
        topgg <- names(sort(rho1,decreasing=TRUE))

        if(0) {
            gs0 <- out$gse[out$gse$module==k,]
            gs0 <- head(gs0[order(gs0$p.value),],10)
            gs.genes <- unique(unlist(strsplit(gs0$genes,split="\\|")))
            gs.genes
            topgg <- intersect(topgg, gs.genes) ## only GSET genes???
        }
        topgg <- head(topgg,ntop)
        topgg

        message("[wgcna_corGraph.RENDER] ncol.xx = ",ncol(xx))
        message("[wgcna_corGraph.RENDER] nrow.xx = ",nrow(xx))
        message("[wgcna_corGraph.RENDER] len.topgg = ",length(topgg))
        
        ##rho <- cor(xx[,topgg])
        ##rho <- try( qgraph::cor_auto(xx[,topgg], forcePD=TRUE) )
        ##rho <- qgraph::cor_auto(xx, forcePD=TRUE)
        rho <- Matrix::nearPD(cor(xx[,topgg]))$mat
        me.color <- out$me.colors[k]
        color1 <- me.color
        color1 <- c("white", me.color)[1 + 1*(colnames(rho)==k)]
        size1  <- c(7,10)[1 + 1*(colnames(rho)==k)]

        
        message("[wgcna_corGraph.RENDER] nrow(rho) = ",nrow(rho))
        

        qgraph::qgraph(rho, graph="glasso", layout="spring", sampleSize=nrow(xx),
               labels = rownames(rho), color=color1,
               tuning = 0,  ## gamma for EBIClasso. 0.5=default, 0=BIC
               vsize = size1, cut=0, maximum=.45,
               border.width=1.5)


        
    })

    wgcna_corGraph.opts <- shiny::tagList(
        ##selectInput(ns('wgcna_corGraph_selmodule'),'module:', choices=NULL)
    )

    wgcna_corGraph_info <- "<b>Correlation network.</b> Partial correlation graph centered on module eigen-gene with top most correlated features. Green edges correspond to positive (partial) correlation, red edges to negative (partial) correlation. Width of the edges is proportional to the correlation strength of the gene pair. The regularized partial correlation matrix is computed using the 'graphical lasso' (Glasso) with BIC model selection."
    
    shiny::callModule(
        plotModule,
        id = "wgcna_corGraph", label = "b",
        func = wgcna_corGraph.RENDER,
        func2 = wgcna_corGraph.RENDER,        
        info.text = wgcna_corGraph_info,
        options = wgcna_corGraph.opts,
        title = "Correlation network",
        ## caption = topEnriched_caption
        caption2 = wgcna_corGraph_info,
        ##pdf.width = 14, pdf.height = 4, 
        height = c(420,650),
        width = c('auto',1000),
        res = c(72,80),
        add.watermark = WATERMARK
    )

    ##-----------------------------------------------------------
    ## TOM UMAP/t-SNE
    ##-----------------------------------------------------------
    
    wgcna_umap.RENDER %<a-% shiny::reactive({

        out <- wgcna.compute()

        method="tsne2d"
        method <- input$clust_method
        
        par(mfrow=c(1,1), mar=c(2,3,1,1))
        me1 <- paste0("ME", out$net$colors)
        pos <- out$clust[[method]]        
        pgx.scatterPlotXY.BASE(pos, var=me1, col=out$me.colors)
        ##WGCNA::TOMplot(dissTOM, net$dendrograms[[1]], net$colors)

    })

    wgcna_umap.opts <- shiny::tagList(
        shiny::selectInput(ns('clust_method'),'method:', choices=c("tsne2d","umap2d","pca2d"))
    )

    wgcna_umap_info <- "<b>TOM umap.</b> UMAP visualization of TOM correlation of genes."
    
    shiny::callModule(
        plotModule,
        id = "wgcna_umap", label = "d",
        func = wgcna_umap.RENDER,
        func2 = wgcna_umap.RENDER,        
        info.text = wgcna_umap_info,
        options = wgcna_umap.opts,
        title = "Gene clustering",
        ## caption = topEnriched_caption
        caption2 = wgcna_umap_info,
        ##pdf.width = 14, pdf.height = 4, 
        pdf.width=5, pdf.height=10, 
        height = c(420,650),
        width = c('auto',650),
        res = c(72,80),
        add.watermark = WATERMARK
    )
    
    ##----------------------------------------
    ##--------- enrichment plot --------------
    ##----------------------------------------
    
    wgcna_eigenClustering.RENDER %<a-% shiny::reactive({
        
        out <- wgcna.compute()
        
        ## MET = orderMEs(cbind(MEs, weight))
        MET = out$net$MEs
        if(NCOL(MET)<=2) MET <- cbind(MET,MET)  ## error if ncol(MET)<=2 !!!!
        
        ## Plot the relationships among the eigengenes and the trait
        ##sizeGrWindow(5,7.5);
        ## par(cex=0.9)
        WGCNA::plotEigengeneNetworks(
                   MET, "", marDendro = c(0,4,1,2),
                   marHeatmap = c(3,4,1,2), cex.lab = 0.8,
                   xLabelsAngle = 90)

    })

    wgcna_eigenClustering_info = "eigenClustering."
    wgcna_eigenClustering_caption = "<b>eigenClustering</b> <b>(a)</b> ..."    
    wgcna_eigenClustering.opts <- shiny::tagList(
        ##selectInput(ns('wgcna_enrichPlot_selmodule'),'module:', choices=NULL)
    )
    
    shiny::callModule(
        plotModule,
        id = "wgcna_eigenClustering", label = "a",
        func = wgcna_eigenClustering.RENDER,
        func2 = wgcna_eigenClustering.RENDER,        
        info.text = wgcna_eigenClustering_info,
        options = wgcna_eigenClustering.opts,
        title = "Eigengene clustering",
        ## caption = topEnriched_caption
        caption2 = wgcna_eigenClustering_info,
        pdf.width = 6, pdf.height = 10, 
        height = c(fullH,700),
        width = c('auto',450),
        res = c(80,90),
        add.watermark = WATERMARK
    )
    
    ##----------------------------------------
    ##--------- eigen barplots ---------------
    ##----------------------------------------
    
    wgcna_eigenCorrelation.RENDER %<a-% shiny::reactive({
    ##wgcna_eigenCorrelation.RENDER <- shiny::reactive({

        message("[wgcna_eigenCorrelation.RENDER] reacted")
        
        out <- wgcna.compute()
        
        MEs <- out$net$MEs
        rho <- cor(MEs, out$datExpr)
        rho[is.na(rho) | is.infinite(rho)] <- 0

        ylab0 = "ME correlation"
        if(input$eigen_cov) {
            sdx <- apply( out$datExpr,2,sd,na.rm=TRUE)
            rho <- t(t(rho) * sdx**2)
            ylab0 = "ME covariance"
        }
        
        n=6
        n  <- nrow(rho)
        nr <- ceiling(sqrt(n))
        nc <- ceiling(n / nr)
        nr
        nc

        message("[wgcna_eigenCorrelation.RENDER] n  = ", n)
        message("[wgcna_eigenCorrelation.RENDER] nr = ", nr)
        message("[wgcna_eigenCorrelation.RENDER] nc = ", nc)
        
        ntop = 15        
        par(mfrow=c(nr,nc), mar=c(6,3.1,2.3,1), oma=c(1,1,1,1)*0, mgp=c(2.1,0.8,0))
        ##par(mfrow=c(nr,nc), mar=c(8,4,3,1)*0, oma=c(1,1,1,1)*0)
        k=1
        me <- names(out$me.colors)  ## sorted
        for(m in me) {
            message("[wgcna_eigenCorrelation.RENDER] plotting ", m," ...")
            i1 <- head(order(rho[m,]),ntop)
            i2 <- tail(order(rho[m,]),ntop)
            barplot( sort(rho[k,c(i1,i2)]),
                    ylab = ylab0, las=3, cex.names=0.90, main=NULL)
            title(m, line=0.3)
        }
        
    })

    wgcna_eigenCorrelation_opts = shiny::tagList(
        shiny::checkboxInput(ns("eigen_cov"),"covariance", FALSE)
    )

    wgcna_eigenCorrelation_info =
        "<b>WGCNA Module membership (eigengene correlation).</b> For each module, we also define a quantitative measure of module membership (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module."

    shiny::callModule(
        plotModule, 
        id = "wgcna_eigenCorrelation", ##ns=ns,
        title="Module membership (eigengene correlation)", label="b",
        func  = wgcna_eigenCorrelation.RENDER,
        func2 = wgcna_eigenCorrelation.RENDER, 
        download.fmt = c("png","pdf"),
        options = wgcna_eigenCorrelation_opts,
        info.text = wgcna_eigenCorrelation_info,        
        height = c(fullH,720), width = c('auto',1050),
        pdf.width=10, pdf.height=6, res=c(90,105),
        add.watermark = WATERMARK
    )

    ##----------------------------------------
    ##------ intramodular analysis -----------
    ##----------------------------------------
    
    wgcna_intraHeatmap.RENDER %<a-% shiny::reactive({
    ##wgcna_intraHeatmap.RENDER <- shiny::reactive({

        message("[wgcna_intraHeatmap.RENDER] reacted")
        
        out <- wgcna.compute()
        
        MEs <- out$net$MEs
        rho1 <- cor(MEs, out$datExpr, use="pairwise")
        rho1[is.na(rho1) | is.infinite(rho1)] <- 0
        
        rho2 <- cor(out$datTraits, out$datExpr, use="pairwise")
        rho2[is.na(rho2) | is.infinite(rho2)] <- 0

        rho3 <- cor( t(rho2), t(rho1), use="pairwise")
        rho3[is.na(rho3) | is.infinite(rho3)] <- 0
        
        gx.heatmap(rho3, nmax=50, mar=c(5,10),
                   keysize=0.5, scale="none", key=FALSE)
        
    })

    wgcna_intraHeatmap_opts = shiny::tagList(
        shiny::checkboxInput(ns("eigen_cov"),"covariance", FALSE)
    )

    wgcna_intraHeatmap_info =
        "<b>WGCNA Module membership (eigengene correlation).</b> For each module, we also define a quantitative measure of module membership (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module."

    shiny::callModule(
        plotModule, 
        id = "wgcna_intraHeatmap", ##ns=ns,
        title="Membership-trait heatmap", label="a",
        func  = wgcna_intraHeatmap.RENDER,
        func2 = wgcna_intraHeatmap.RENDER, 
        download.fmt = c("png","pdf"),
        options = wgcna_intraHeatmap_opts,
        info.text = wgcna_intraHeatmap_info,        
        height = c(fullH,720), width = c('auto',1050),
        pdf.width=6, pdf.height=9, res=c(85,100),
        add.watermark = WATERMARK
    )


    ##----------------------------------------
    ##------ intramodular scatter ------------
    ##----------------------------------------
    
    wgcna_intraScatter.RENDER %<a-% shiny::reactive({
    ##wgcna_intraScatter.RENDER <- shiny::reactive({

        message("[wgcna_intraScatter.RENDER] reacted")
        
        out <- wgcna.compute()
        
        MEs <- out$net$MEs
        rho1 <- cor(MEs, out$datExpr, use="pairwise")
        rho1[is.na(rho1) | is.infinite(rho1)] <- 0
        
        rho2 <- cor(out$datTraits, out$datExpr, use="pairwise")
        rho2[is.na(rho2) | is.infinite(rho2)] <- 0
        
        rho3 <- cor( t(rho2), t(rho1), use="pairwise")
        rho3[is.na(rho3) | is.infinite(rho3)] <- 0

        k="ME1"
        k = input$selected_module
        in.mod <- colnames(rho1) %in% out$me.genes[[k]]
        table(in.mod)
        col1 <- c("grey60",out$me.colors[k])[1 + 1*in.mod]

        ntop <- ifelse(nrow(rho3)>=20, 20, 12)
        top.px <- head(order(-abs(rho3[,k])), ntop)
        if(ntop==20) mfrow0 <- c(4,5)
        if(ntop==12) mfrow0 <- c(3,4)
        
        par(mfrow=mfrow0, mar=c(4,4,2,1), mgp=c(2.0,0.8,0))
        i=top.px[1]
        for(i in top.px) {
            KEGGgraph::plot( rho1[k,], rho2[i,], pch=20, cex=0.7, col=col1,
                 xlab = "Module membership (eigengene cor)",
                 ylab = "Gene significance (trait cor)")
            title(paste(k,"vs.",paste(rownames(rho2)[i])), cex=1)                  
        }
        
    })

    wgcna_intraScatter_opts = shiny::tagList(
        ## shiny::checkboxInput(ns("eigen_cov"),"covariance", FALSE)
    )

    wgcna_intraScatter_info =
        "<b>WGCNA Module membership (eigengene correlation).</b> For each module, we also define a quantitative measure of module membership (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module."

    shiny::callModule(
        plotModule, 
        id = "wgcna_intraScatter", ##ns=ns,
        title="Membership vs. trait correlation", label="b",
        func  = wgcna_intraScatter.RENDER,
        func2 = wgcna_intraScatter.RENDER, 
        download.fmt = c("png","pdf"),
        options = wgcna_intraScatter_opts,
        info.text = wgcna_intraScatter_info,        
        height = c(fullH,720), width = c('auto',1150),
        pdf.width=12, pdf.height=9, res=c(85,90),
        add.watermark = WATERMARK
    )
    
    
    
} ## end of Board
