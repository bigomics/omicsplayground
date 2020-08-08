ClusteringInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
    )
}

ClusteringUI <- function(id) {
    ns <- NS(id)  ## namespace
    fillRow(
        flex = c(1.4,0.10,1),
        height = 780,
        tabsetPanel(
            id = ns("tabs1"),
            tabPanel("Heatmap",uiOutput(ns("hm_heatmap_UI"))),
            tabPanel("PCA/tSNE",uiOutput(ns("hm_pcaUI"))),
            tabPanel("Parallel",uiOutput(ns("hm_parcoordUI")))
        ),
        br(),
        tabsetPanel(
            id = ns("tabs2"),
            tabPanel("Annotate clusters",uiOutput(ns("hm_annotateUI"))),
            tabPanel("Phenotypes",uiOutput(ns("hm_phenoplotUI"))),
            tabPanel("Feature ranking",uiOutput(ns("hm_featurerankUI")))      
        )
    )
}


ClusteringModule <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE

    inputData <- env[["load"]][["inputData"]]

    fullH = 750  ## full height of page
    
    description = "<b>Cluster Analysis.</b> Discover clusters of similar genes or samples using unsupervised machine learning."
    output$description <- renderUI(HTML(description))
    
    clust_infotext = paste('
The <strong>Cluster Analysis</strong> module performs unsupervised clustering analysis of the data. After having done the QC, it is probably the first way to explore your data. The main purpose is to discover patterns and subgroups in the data, show correlation with known phenotypes, detect outliers, or investigate batch effects.

<br><br>In the <strong>Heatmap</strong> panel hierarchical clustering can be performed on gene level or gene set level (selected under the {Level} dropdown list). During the heatmap generation, the platform provides a functional annotation for each feature cluster in <strong>Annotate cluster</strong> panel. Users can select from a variety of annotation databases from the literature, such as ',a_MSigDB,', ',a_KEGG,' and ',a_GO,'. The <strong>PCA/tSNE</strong> panel shows unsupervised clustering of the samples in 2D/3D as obtained by ',a_PCA,' or ',a_tSNE,' algorithms. The <strong>Phenotypes</strong> panel on the right, shows the phenotype distribution as colors on the t-SNE plot. 

<br><br>EXPERT MODE ONLY: The <strong>Feature ranking</strong> panel computes a discriminant score for gene (or geneset) families. This allows to investigate what family of genes (or gene sets) can best discriminate the groups. 

<br><br><br><br>
<center><iframe width="500" height="333" src="https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=2" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>

')

    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================

    output$inputsUI <- renderUI({
        ui <- tagList(
            ##tipify( actionLink(ns("clust_info"), "Info", icon = icon("info-circle")),
            ##       "Show more information about this module."),
            tipify( actionLink(ns("clust_info"), "Tutorial", icon = icon("youtube")),
                   "Show more information and video tutorial about this module."),
            hr(), br(),             
            tipify( selectInput(ns("hm_level"),"Level:", choices=c("gene","geneset")),
                   "Specify the level analysis: gene or geneset level.", placement="top"),
            tipify( selectInput(ns("hm_features"),"Features:", choices=NULL, multiple=FALSE),
                   "Select a family of features.", placement="top"),
            conditionalPanel(
                "input.hm_features == '<custom>'", ns=ns,
                tipify( textAreaInput(ns("hm_customfeatures"), NULL, value = NULL,
                                      rows=5, placeholder="Paste your custom gene list"),
                       "Paste a custom list of genes to be used as features.", placement="bottom")
            ),
            tipify( selectInput(ns("hm_samplefilter"),"Filter samples:", choices=NULL, multiple=TRUE),
                   "Filter the relevant samples for the analysis.", placement="top"),
            tipify( actionLink(ns("hm_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Toggle advanced options.", placement="top"),
            br(),
            conditionalPanel(
                "input.hm_options % 2 == 1", ns=ns,
                tagList(
                    tipify( checkboxInput(ns('hm_group'),'group by condition',FALSE),
                           "Group the samples by condition.", placement="bottom")
                ),
                tagList(
                    tipify( checkboxInput(ns('hm_excludedown'),'exclude down genesets',FALSE),
                           "Exclude down genesets.", placement="bottom")
                )
            )
        )
        if(USER_MODE=="DEV") {
            ui1 <- tagList(
                br(),hr(),
                h5("Developer options:"),
                radioButtons(ns('hm_gsetmatrix'),'Gene set matrix (dev):',
                             choices=c("meta"), inline=TRUE)
            )
            ui <- c(ui, ui1)
        }
        ui
    })
    outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================

    observeEvent( input$clust_info, {
        showModal(modalDialog(
            title = HTML("<strong>Clustering Module</strong>"),
            HTML(clust_infotext),
            easyClose = TRUE, size="l" ))
    })
       
    ## update filter choices upon change of data set 
    observe({
        ngs <- inputData()
        req(ngs)
        
        levels = getLevels(ngs$Y)

        dbg("[ClusteringModule::observe] updating hm_samplefilter")
        dbg("[ClusteringModule::observe] head(levels)=",head(levels))

        updateSelectInput(session, "hm_samplefilter", choices=levels)
        
        if(DEV.VERSION && !is.null(ngs$gset.meta$matrices) ) {
            jj = which(!sapply(ngs$gset.meta$matrices,is.null))
            mat.names = names(ngs$gset.meta$matrices)[jj]
            updateRadioButtons(session, "hm_gsetmatrix", choices=mat.names,
                               selected="meta", inline=TRUE)    
            ##updateRadioButtons(session, "hm_gsetmatrix", choices=c("meta"), inline=TRUE)    
        }
        
        updateRadioButtons(session, "hm_splitby", selected='none')

        ## update defaults??
        lab <- c("row","column")
        if(ncol(ngs$X) > 80) lab <- c("row")
        ##if(nrow(ngs$X)<100) lab <- c(lab,"row")
        updateCheckboxGroupInput(session,"hm_showlabel", selected=lab)        
    })

    observeEvent( input$hm_splitby, {
        ngs <- inputData()
        req(ngs)
        if(input$hm_splitby=='none') return()
        if(input$hm_splitby=='gene') {
            xgenes <- sort(rownames(ngs$X))
            updateSelectizeInput(session, "hm_splitvar", choices=xgenes, server=TRUE)
        }
        if(input$hm_splitby=='phenotype') {
            cvar <- sort(pgx.getCategoricalPhenotypes(ngs$samples, min.ncat=2, max.ncat=20))
            sel <- cvar[1]
            cvar0 <- grep("^[.]",cvar,value=TRUE,invert=TRUE) ## no estimated vars
            sel <- head(c(grep("type|family|class|stat",cvar0,ignore.case=TRUE,value=TRUE),
                          cvar0,cvar),1)
            updateSelectInput(session, "hm_splitvar", choices=cvar, selected=sel)
        }
    })
    
    input_hm_samplefilter <- reactive({
        input$hm_samplefilter
    }) %>% debounce(2000)
    
    ## update choices upon change of level
    ##observeEvent(input$hm_level, {
    observe({
        ngs = inputData()
        req(ngs)
        req(input$hm_level)
        ###if(is.null(input$hm_level)) return(NULL)
        choices = names(ngs$families)

        dbg("[ClusteringModule::observe] input$hm_level=",input$hm_level)
        dbg("[ClusteringModule::observe] dim(ngs$counts)=",dim(ngs$counts))
        dbg("[ClusteringModule::observe] dim(ngs$gsetX)=",dim(ngs$gsetX))
        dbg("[ClusteringModule::observe] length(choices)=",length(choices))
        dbg("[ClusteringModule::observe] length(COLLECTIONS)=",length(COLLECTIONS))
        
        if(input$hm_level=="geneset") {
            dbg("[ClusteringModule::observe] change to COLLECTIONS")            
            nk <- sapply(COLLECTIONS, function(k) sum(k %in% rownames(ngs$gsetX)))
            choices = names(COLLECTIONS)[nk>=5]
        }
        choices <- c("<custom>",choices)
        choices <- sort(unique(choices))

        dbg("[ClusteringModule::observe] updating hm_features")
        dbg("[ClusteringModule::observe] length(choices)=",length(choices))
        dbg("[ClusteringModule::observe] head(choices)=",head(choices))
        updateSelectInput(session, "hm_features", choices=choices)
    })
    

    ##================================================================================
    ##========================= REACTIVE FUNCTIONS ===================================
    ##================================================================================
    
    getFilteredMatrix <- reactive({
        ## Returns filtered matrix ready for clustering. Filtering based
        ## on user selected geneset/features or custom list of genes.
        ##
        ##
        ##
        ngs <- inputData()
        req(ngs)
        genes = as.character(ngs$genes$gene_name)
        genesets = rownames(ngs$gsetX)

        if(input$hm_level=="geneset") {
            ##-----------------------------------
            ## Gene set level features
            ##-----------------------------------
            gsets = rownames(ngs$gsetX)
            gsets = unique(unlist(COLLECTIONS[input$hm_features]))
            zx = ngs$gsetX
            if(DEV.VERSION && !is.null(ngs$gset.meta$matrices) ) {
                k <- input$hm_gsetmatrix
                zx <- as.matrix(ngs$gset.meta$matrices[[k]])
            }            
            if(input$hm_customfeatures!="") {
                gsets1 = genesets[grep(input$hm_customfeatures, genesets,ignore.case=TRUE)]
                if(length(gsets1)>2) gsets = gsets1
            }        
            zx = zx[intersect(gsets,rownames(zx)),]
        } else {
            ##-----------------------------------
            ## Gene level features
            ##-----------------------------------
            pp = head(rownames(ngs$X),400)
            gg = ngs$families[[1]]
            if(input$hm_features =="<all>") {
                gg = rownames(ngs$X)
            } else if(input$hm_features %in% names(ngs$families)) {
                gg = ngs$families[[input$hm_features]]
            } else {
                gg = rownames(ngs$X)
            }
            if(input$hm_customfeatures!="") {
                gg1 = strsplit(input$hm_customfeatures,split="[, ;\n\t]")[[1]]
                if(length(gg1)==1) gg1 <- paste0(gg1,"*")
                gg1 = gsub("[ \n\t]","",gg1)
                starred = grep("[*]",gg1)
                if(length(starred)>0) {
                    gg2 = lapply(gg1[starred], function(a)
                        genes[grep(paste0("^",sub("[*]","",a)),genes,ignore.case=TRUE)])
                    gg1 = unique(c(gg1,unlist(gg2)))
                }
                gg1 = intersect(gg1,ngs$genes$gene_name)            
                if(length(gg1)>1) gg = gg1
            }
            jj <- match(toupper(gg), toupper(ngs$genes$gene_name))
            pp = rownames(ngs$genes)[setdiff(jj,NA)]
            zx = ngs$X[intersect(pp,rownames(ngs$X)),]
        }
        if(nrow(zx)==0) return(NULL)

        dim(zx)
        kk <- selectSamplesFromSelectedLevels(ngs$Y, input_hm_samplefilter() )
        zx <- zx[,kk,drop=FALSE]    
        
        if( input$hm_level=="gene" && "chr" %in% names(ngs$genes)
           && input$hm_filterXY ) {
            ## Filter out X/Y chromosomes before clustering
            chr.col <- grep("^chr$|^chrom$",colnames(ngs$genes))
            chr <- ngs$genes[rownames(zx),chr.col]
            not.xy <- !(chr %in% c("X","Y",23,24)) & !grepl("^X|^Y|chrX|chrY",chr)
            table(not.xy)
            zx <- zx[which(not.xy), ]
            dim(zx)
        }

        if(input$hm_excludedown && input$hm_level=="geneset") {
            sel <- grep("dn|down",rownames(zx),ignore.case=TRUE,invert=TRUE)
            zx <- zx[sel,,drop=FALSE]
        }
        
        return(zx)
    })

    
    getTopMatrix <- reactive({

        ngs <- inputData()
        req(ngs)
        zx <- getFilteredMatrix()

        dbg("[ClusteringModule:getTopMatrix] dim(ngs$counts) = ",dim(ngs$counts))
        dbg("[ClusteringModule:getTopMatrix] dim(zx) = ",dim(zx))
        
        nmax = 4000
        nmax = as.integer(input$hm_ntop)
        idx <- NULL
        splitvar ="none"
        splitvar <- input$hm_splitvar
        splitby  <- input$hm_splitby
        do.split <- splitby!='none'
        
        if(splitby=="gene" && !splitvar %in% rownames(ngs$X)) return(NULL)
        if(splitby=="phenotype" && !splitvar %in% colnames(ngs$samples)) return(NULL)

        grp <- NULL
        if(do.split && splitvar %in% colnames(ngs$samples)) {
            dbg("[ClusteringModule:getTopMatrix] splitting by phenotype: ",splitvar)
            grp <- ngs$samples[colnames(zx),splitvar]
        }
        table(grp)

        if(do.split && splitvar %in% rownames(ngs$X)) {
            dbg("[ClusteringModule:getTopMatrix] splitting by gene: ",splitvar)
            ##xgene <- rownames(ngs$X)[1]
            gx <- ngs$X[splitvar,]
            if(ncol(zx)>50) {
                qq <- quantile(gx, probs=c(0.33,0.66), na.rm=TRUE)
                grp <- cut(gx, breaks=c(-999,qq,999), label=c("low","medium","high"))
            } else {
                qq <- quantile(gx, probs=c(0.50), na.rm=TRUE)
                grp <- cut(gx, breaks=c(-999,qq,999), label=c("low","high"))
            }
            grp <- paste0(splitvar,":",grp)
        }
        
        dbg("[ClusteringModule:getTopMatrix] len.grp=",length(grp))
        dbg("[ClusteringModule:getTopMatrix] splitby=",splitby)
        
        ##if(length(grp)==0) splitby <- 'none'
        if(do.split && length(grp)==0) return(NULL)        
        
        topmode="specific"
        topmode <- input$hm_topmode
        if(topmode == "specific" && length(table(grp))<=1) {
            topmode <- "sd"
        }
        if(!do.split && topmode=="specific") topmode <- "sd"
        
        if(topmode=="pca") {
            dbg("[ClusteringModule:getTopMatrix] splitting by PCA")
            require(irlba)
            NPCA=5
            svdres <- irlba(zx - rowMeans(zx), nv=NPCA)
            ntop = 12
            ntop <- as.integer(input$hm_ntop) / NPCA
            gg <- rownames(zx)
            sv.top <- lapply(1:NPCA,function(i) gg[head(order(-abs(svdres$u[,i])),ntop)] )
            gg.top <- unlist(sv.top)
            for(i in 1:length(sv.top)) {
                sv.top[[i]] <- paste0("PC",i,":",sv.top[[i]])
            }
            sv.top1 <- unlist(sv.top)
            zx <- zx[gg.top,,drop=FALSE]
            ##rownames(zx) <- sv.top1
            dim(zx)
            ##idx <- paste0("PC",sub(":.*","",sv.top1))
            idx <- sub(":.*","",sv.top1)
            table(idx)
        } else if(topmode=="specific" && splitby!="none") {
            ##
            ## sample cluster specifice gene prioritazion
            ##
            ##grp <- ngs$samples[colnames(zx),"cluster"]
            table(grp)

            dbg("[ClusteringModule:getTopMatrix] head.grp=",head(grp))
            dbg("[ClusteringModule:getTopMatrix] dim.zx=",dim(zx))

            grp.zx <- t(apply(zx, 1, function(x) tapply(x, grp, mean)))
            if(length(table(grp))==1) {
                grp.zx <- t(grp.zx)
                colnames(grp.zx) <- paste0(splitvar,":",grp[1])                
            }
            grp.dx <- grp.zx*0
            nc <- ncol(grp.dx)
            for(i in 1:nc) {
                grp.dx[,i] <- grp.zx[,i] - rowMeans(grp.zx[,-i,drop=FALSE])
            }
            gg <- rownames(zx)
            ntop = 12
            ntop <- ceiling( as.integer(input$hm_ntop) / ncol(grp.dx) )
            grp.top <- lapply(1:nc,function(i) gg[head(order(-grp.dx[,i]),ntop)] )
            ##idx <- unlist(sapply(1:nc,function(i) rep(i,length(grp.top[[i]]))))
            idx <- unlist(mapply(rep,1:nc,sapply(grp.top,length)))
            idx <- paste0("M",idx)
            table(idx)        
            gg.top <- unlist(grp.top)
            zx <- zx[gg.top,,drop=FALSE]
            dim(zx)
        } else {
            ## Order by SD
            ## 
            dbg("[ClusteringModule:getTopMatrix] order by SD")
            ii <- order(-apply(zx,1,sd,na.rm=TRUE))
            zx = zx[ii,,drop=FALSE] ## order
            zx = head(zx,nmax)
        }
        ##zx = zx / apply(zx,1,sd,na.rm=TRUE)  ## scale??
        
        ## ------------- cluster the genes???
        if(is.null(idx)) {
            D <- as.dist(1 - cor(t(zx),use="pairwise"))
            system.time( hc <- fastcluster::hclust(D, method="ward.D2" ) )
            ## system.time( hc <- flashClust::hclust(D, method="ward" ) )
            ## system.time( hc <- nclust(D, link="ward") )
            ngrp = 4  ## how many default groups???
            idx = paste0("S",cutree(hc, ngrp))
        }

        ## ------------- matched annotation
        annot = ngs$Y[colnames(zx),,drop=FALSE]
        kk = grep("sample|patient",colnames(annot),invert=TRUE)
        annot = annot[,kk,drop=FALSE]  ## no group??    
        samples = colnames(zx) ## original sample list
        
        ## ------------- calculate group summarized (for what???)
        if(1) { 
            grpvar = "group"
            gg.grp = ngs$samples[colnames(zx),grpvar]
            if(0) {
                grp.zx = t(apply( zx, 1, function(x) tapply(x, gg.grp, mean)))
                grp.annot = annot[match(colnames(zx),gg.grp),,drop=FALSE] ## NEED RETHINK!!!
                rownames(grp.annot) = colnames(zx)
            } else {
                grp.zx = tapply( colnames(zx), gg.grp, function(k) rowMeans(zx[,k,drop=FALSE],na.rm=TRUE))
                grp.zx = do.call( cbind, grp.zx)
                ## take most frequent term as group annotation value
                most.freq <- function(x) names(sort(-table(x)))[1]
                grp.annot = tapply( rownames(annot), gg.grp, function(k) {
                    f <- apply(annot[k,,drop=FALSE],2,function(x) most.freq(x))
                    w.null <- sapply(f,is.null)
                    if(any(w.null))  f[which(w.null)] <- NA
                    unlist(f)
                })
                grp.annot = data.frame(do.call( rbind, grp.annot))
                grp.annot = grp.annot[colnames(grp.zx),]
            }
        }

        dbg("[ClusteringModule:getTopMatrix] done!")
        
        ##input$top_terms
        filt <- list(mat=zx, annot=annot, grp=grp,
                     grp.mat = grp.zx, grp.annot = grp.annot,
                     clust=idx, samples=samples)
        return(filt)
    })


    ##================================================================================
    ##========================= FUNCTIONS ============================================
    ##================================================================================


    ## Heatmap

    hm_splitmap_text = tagsub("Under the <strong>Heatmap</strong> panel, hierarchical clustering can be performed on gene level or gene set level expression in which users have to specify it under the {Level} dropdown list. <p>Under the plot configuration {{Settings}}, users can split the samples by a phenotype class (e.g., tissue, cell type, or gender) using the {split by} setting. In addition, users can specify the top N = (50, 150, 500) features to be used in the heatmap. The ordering of top features is selected under {top mode}. The criteria to select the top features are: <ol><li>SD - features with the highest standard deviation across all the samples, </li><li>specific - features that are overexpressed in each phenotype class compared to the rest, or by </li><li>PCA - by principal components.<br></ol> <br><p>Users can also choose between 'relative' or 'absolute' expression scale. Under the {CexCol} and {CexRow} settings, it is also possible to adjust the cex for the column and row labels.")


    require(plotly)
    require(scatterD3)
    ##require(tippy)
    
    observe({
        ngs <- inputData()
        req(ngs)

        cvar <- sort(pgx.getCategoricalPhenotypes(ngs$samples, min.ncat=2, max.ncat=20))
        updateSelectInput(session, "hm_pcvar", choices=cvar, selected="group")
        updateSelectInput(session, "hm2_pcvar", choices=cvar, selected="group")        
    })
    
    ##hm1_splitmap.RENDER %<a-% reactive({
    hm1_splitmap.RENDER <- reactive({    
        
        ngs <- inputData()
        alertDataLoaded(session,ngs)
        req(ngs)
        
        filt <- getTopMatrix()        
        req(filt)
        ##if(is.null(filt)) return(NULL)
        
        if(input$hm_group) {
            zx <- filt$grp.mat    
            annot = filt$grp.annot
        } else {
            zx <- filt$mat    
            annot = filt$annot        
        }
        zx.clust <- filt$clust    

        dbg("hm1_splitmap.RENDER: 1: dim(zx)=", dim(zx))
        dbg("hm1_splitmap.RENDER: 1: sum.dup.rownames.zx=", sum(duplicated(rownames(zx))))
        
        show_rownames = TRUE
        ##if(nrow(zx) > 100) rownames(zx) = rep("",nrow(zx))
        if(nrow(zx) > 100) show_rownames = FALSE   

        main = "  "
        main = input$hm_features
        
        cex1 = ifelse(ncol(zx)>50,0.75,1)
        cex1 = ifelse(ncol(zx)>100,0.5,cex1)
        cex1 = ifelse(ncol(zx)>200,0,cex1)
        
        scale.mode = "row.center"
        scale.mode = ifelse(input$hm_scale=="relative","row.center","none")
        scale.mode
        
        ## split genes dimension in 5 groups
        splity = 5
        splity = 6
        if(!is.null(zx.clust)) splity = zx.clust
        
        ## split samples
        splitx = NULL    
        splitx = filt$grp

        show_legend=show_colnames=TRUE
        if(input$hm_level=="geneset" || !is.null(splitx)) show_legend = FALSE
        annot$group = NULL  ## no group in annotation
        show_colnames <- ("column" %in% input$hm_showlabel)
        ##if(ncol(zx) > 200) show_colnames <- FALSE ## never...    
        
        if(input$hm_level=="gene") {
            rownames(zx) = sub(".*:","",rownames(zx))
        }
        rownames(zx) <- sub("HALLMARK:HALLMARK_","HALLMARK:",rownames(zx))
        rownames(zx) = gsub(GSET.PREFIX.REGEX,"",rownames(zx))
        rownames(zx) = substring(rownames(zx),1,50)  ## cut long names...
        if(input$hm_level=="geneset")  rownames(zx) <- tolower(rownames(zx))
        
        cex2 <- ifelse( nrow(zx) > 60, 0.8, 0.9)    
        cex1 <- as.numeric(input$hm_cexCol)*0.85
        cex2 <- as.numeric(input$hm_cexRow)*0.75
        
        crot <- 0
        totnchar <- nchar(paste0(unique(splitx),collapse=""))
        totnchar
        nx <- length(unique(splitx))
        if(!is.null(splitx) & (totnchar > 44 || nx>=6) ) crot=90
        if(main=="<all>") main <- "  "

        nrownames = 60
        nrownames = ifelse("row" %in% input$hm_showlabel,60,0)
        
        if(0) {
            split=splity;splitx=splitx;mar=c(5,25); scale=scale.mode; show_legend=show_legend;
            show_colnames = show_colnames; column_title_rot=crot;
            show_rownames = nrownames; softmax=0;
            ## side.height.fraction=0.03+0.055*NCOL(annot); 
            labRow=NULL; cexCol=cex1; cexRow=cex2; 
            col.annot=annot; row.annot=NULL; annot.ht=2.2;
            main=main; nmax=-1
        }

        dbg('rendering heatmap: dim(zx)=',dim(zx))
        dbg('rendering heatmap: table(splitx)=',table(splitx))
        dbg('rendering heatmap: table(splity)=',table(splity))
        dbg('rendering heatmap: calling gx.splitmap() ... ')
        
        showNotification('rendering heatmap...')

        gx.splitmap( zx, 
                    split=splity, splitx=splitx,
                    mar=c(15,25),
                    scale=scale.mode, show_legend=show_legend,
                    show_colnames = show_colnames, column_title_rot=crot,
                    show_rownames = nrownames, softmax=0,
                    ## side.height.fraction=0.03+0.055*NCOL(annot), 
                    labRow=NULL, cexCol=cex1, cexRow=cex2, 
                    col.annot=annot, row.annot=NULL, annot.ht=2.2,
                    main=main, nmax=-1)
    })

    hm2_splitmap.RENDER <- reactive({

        ngs <- inputData()
        req(ngs)
        
        cat("<module-clustering:hm2_splitmap.RENDER> reacted\n")

        ## -------------- variable to split samples        
        scale = ifelse(input$hm_scale=="relative","row.center","none")    
        plt <- NULL

        filt <- getTopMatrix()
        ##if(is.null(filt)) return(NULL)
        req(filt)
        
        if(input$hm_group) {
            X <- filt$grp.mat    
            annot = filt$grp.annot
        } else {
            X <- filt$mat    
            annot = filt$annot        
        }
        idx <- filt$clust    
        
        ## sample clustering index
        splitx <- NULL
        splitx <- filt$grp
        
        ## iheatmapr needs factors for sharing between groups
        annotF <- data.frame(as.list(annot),stringsAsFactors=TRUE)
        rownames(annotF) = rownames(annot)
        
        colcex <- as.numeric(input$hm_cexCol)
        rowcex = as.numeric(input$hm_cexRow)
        rowcex = ifelse("row" %in% input$hm_showlabel, rowcex, 0)
        colcex = ifelse("column" %in% input$hm_showlabel, colcex, 0)
        
        tooltips = NULL
        if(input$hm_level=="gene") {
            getInfo <- function(g) {
                aa = paste0("<b>",ngs$genes[g,"gene_name"],"</b>. ",
                            ## ngs$genes[g,"map"],". ",
                            ngs$genes[g,"gene_title"],".")
                breakstring2(aa, 50, brk="<br>")
            }
            tooltips = sapply(rownames(X), getInfo)
        } else {
            aa = gsub("_"," ",rownames(X)) ## just geneset names
            tooltips = breakstring2(aa, 50, brk="<br>")
        }
        ##genetips = rownames(X)
        
        cat("<module-clustering:hm2_splitmap.RENDER> plotting\n")            

        showNotification('rendering iHeatmap...')

        plt <- pgx.splitHeatmapFromMatrix(
            X=X, annot=annotF, ytips=tooltips,
            idx=idx, splitx=splitx, scale=scale,
            row_annot_width=0.03, rowcex=rowcex,
            colcex=colcex )
       
        cat("<module-clustering:hm2_splitmap.RENDER> done\n")
        
        ## DOES NOT WORK...
        ##plt <- plt %>%
        ## config(toImageButtonOptions = list(format='svg', height=800, width=800))
        
        return(plt)
    })

    topmodes <- c("specific","sd","pca")
    ##if(DEV.VERSION) topmodes <- c("sd","specific","pca")
    
    hm_splitmap_opts = tagList(
        tipify( radioButtons(ns("hm_plottype"), "Plot type:",
                             choices=c("ComplexHeatmap","iHeatmap"),
                             selected="ComplexHeatmap", inline=TRUE, width='100%'),
               "Choose plot type: ComplexHeatmap (static) or iHeatmap (interactive)",
               placement="right",options = list(container = "body")),
        tipify( radioButtons(
            ns("hm_splitby"), "split by:", inline=TRUE,
            ## selected="phenotype",
            choices=c("none","phenotype","gene")),
               "Split the samples by phenotype or expression level of a gene.",
               placement="right",options = list(container = "body")),
        conditionalPanel(
            "input.hm_splitby != 'none'", ns=ns,
            tipify( selectInput(ns("hm_splitvar"), NULL, choices=""),
                   "Specify phenotype or gene for splitting the columns of the heatmap.",
                   placement="right",options = list(container = "body")),
        ),
        fillRow(
            height = 50,
            tipify( selectInput(ns('hm_topmode'),'Top mode:',topmodes, width='100%'),
                   "Specify the criteria for selecting top features to be shown in the heatmap.",
                   placement="right", options = list(container = "body")),
            tipify( selectInput(ns('hm_ntop'),'Top N:',c(50,150,500),selected=50),
                   "Select the number of top features in the heatmap.", placement="right")
        ),
        tipify( checkboxInput(ns('hm_filterXY'),'Exclude X/Y genes',FALSE),
               "Exclude genes on X/Y chromosomes.", placement="right"),
        tipify( radioButtons(
            ns('hm_scale'), 'Scale:', choices=c('relative','absolute'), inline=TRUE),
            "Show relative (i.e. mean-centered) or absolute expression values.",
            placement="bottom"),
        tipify( checkboxGroupInput(ns('hm_showlabel'), 'Show labels:', choices=c('row','column'),
                                   selected=c('row','column'), inline=TRUE),
               "Show/hide row or column names.", placement="bottom"),
        fillRow(
            height = 50,
            tipify( selectInput(ns("hm_cexCol"), "CexCol:", choices=seq(0.1,1.2,0.1), selected=1, width='100%'),
                   "Specify the column label cex.", placement="right",options = list(container = "body")),
            tipify( selectInput(ns("hm_cexRow"), "CexRow:", choices=seq(0.1,1.2,0.1), selected=1, width='100%'),
                   "Specify the row label cex.", placement="right",options = list(container = "body"))
        ),
        br()
    )
    
    hm_splitmap_caption = "<b>Clustered heatmap.</b> Heatmap showing gene expression sorted by 2-way hierarchical clustering. Red corresponds to overexpression, blue to underexpression of the gene. At the same time, gene clusters are functionally annotated in the 'Annotate clusters' panel on the right."

    hm_splitmap_captionDN <- reactive({
        text1 = "<b>Clustered heatmap.</b> Heatmap showing gene expression sorted by 2-way hierarchical clustering. Red corresponds to overexpression, blue to underexpression of the gene. At the same time, gene clusters are functionally annotated in the 'Annotate clusters' panel on the right."
        return(text1)
    })
    
    output$hm1_splitmap <- renderPlot({
        hm1_splitmap.RENDER()
        ## plot(sin)
    }, res=90)
    
    output$hm2_splitmap <- renderIheatmap({
        hm2_splitmap.RENDER()
    })

    hm_splitmap.switchRENDER <- reactive({    
        cat("<module_intersect::hm_splitmap.switchRENDER\n")
        ##req(input$hm_plottype)
        p = NULL
        if(input$hm_plottype %in% c("ComplexHeatmap","static") ) {
            p = plotOutput(ns("hm1_splitmap"), height=fullH-100)  ## height defined here!!
        } else {
            p = iheatmaprOutput(ns("hm2_splitmap"), height=fullH-100) ## height defined here!!
        }
        return(p)
    })
    
    ##output$hm_splitmap_pdf <- downloadHandler(
    hm_splitmap_downloadPDF <- downloadHandler(
        filename = "plot.pdf",
        content = function(file) {
            ##PDFFILE = hm_splitmap_module$.tmpfile["pdf"]  ## from above!
            PDFFILE = paste0(gsub("file","plot",tempfile()),".pdf")            
            dbg("hm_splitmap_pdf:: exporting SWITCH to PDF...")
            ##showNotification("exporting to PDF")            
            if(input$hm_plottype %in% c("ComplexHeatmap","static")) {
                pdf(PDFFILE, width=10, height=8)
                print(hm1_splitmap.RENDER())
                dev.off()
            } else {
                save_iheatmap(hm2_splitmap.RENDER(), filename=PDFFILE,
                              vwidth=1000, vheight=800)
            }
            if(WATERMARK) {
                cat("adding watermark to PDF...\n")
                addWatermark.PDF(PDFFILE) ## from pgx-modules.R
            }            
            dbg("hm_splitmap_pdf:: exporting done...")
            file.copy(PDFFILE,file)        
        }
    )

    hm_splitmap_downloadPNG <- downloadHandler(
        filename = "plot.png",
        content = function(file) {
            PNGFILE = paste0(gsub("file","plot",tempfile()),".png")            
            dbg("hm_splitmap_pdf:: exporting SWITCH to PNG...")
            ##showNotification("exporting to PNG")
            if(input$hm_plottype %in% c("ComplexHeatmap","static")) {
                png(PNGFILE, width=600, height=600, pointsize=14)
                print(hm1_splitmap.RENDER())  ## should be done inside render for base plot...
                dev.off()
            } else {
                save_iheatmap(hm2_splitmap.RENDER(), filename=PNGFILE,  vwidth=1000, vheight=800)
            }
            dbg("hm_splitmap_png:: exporting done...")
            file.copy(PNGFILE,file)        
        }
    )

    hm_splitmap_downloadHTML <- downloadHandler(
        filename = "plot.html",
        content = function(file) {
            ##HTMLFILE = hm_splitmap_module$.tmpfile["html"]  ## from above!
            HTMLFILE = paste0(gsub("file","plot",tempfile()),".html")            
            dbg("renderIheatmap:: exporting SWITCH to HTML...")
            withProgress({
                ##write("<body>HTML export error</body>", file=HTMLFILE)    
                p <- hm2_splitmap.RENDER()
                incProgress(0.5)
                save_iheatmap(p, filename=HTMLFILE)
            }, message="exporting to HTML", value=0 )
            dbg("renderIheatmap:: ... exporting done")
            file.copy(HTMLFILE,file)        
        }
    )
    
    ##    hm_splitmap_module <- plotModule(
    ##        "hm_splitmap",
    hm_splitmap_module <- callModule(
        plotModule,
        id = "hm_splitmap",
        func  = hm_splitmap.switchRENDER, ## ns=ns,
        ## func2 = hm_splitmap.switchRENDER, ## ns=ns,
        show.maximize = FALSE,
        plotlib = "generic",
        renderFunc="renderUI",
        outputFunc="uiOutput",
        download.fmt = c("pdf","png"),
        options = hm_splitmap_opts,
        height = fullH-100, width='100%',
        pdf.width=10, pdf.height=8, 
        title="Clustered Heatmap",
        info.text = hm_splitmap_text,
        info.width="350px",
        download.pdf = hm_splitmap_downloadPDF,
        download.png = hm_splitmap_downloadPNG,
        download.html = hm_splitmap_downloadHTML
        ##caption = hm_splitmap_caption
    )

    output$hm_heatmap_UI <- renderUI({
        fillCol(
            flex = c(NA,0.025,1),
            height = fullH,
            div(HTML(hm_splitmap_caption), class="caption"),
            br(),
            plotWidget(ns("hm_splitmap"))
        )
    })
    outputOptions(output, "hm_heatmap_UI", suspendWhenHidden=FALSE) ## important!!!
    
    ##================================================================================
    ##================================ PCA/tSNE ======================================
    ##================================================================================
    
    hm_PCAplot_text = tagsub(paste0(' The <b>PCA/tSNE</b> panel visualizes unsupervised clustering obtained by the principal components analysis (',a_PCA,') or t-distributed stochastic embedding (',a_tSNE,') algorithms. This plot shows the relationship (or similarity) between the samples for visual analytics, where similarity is visualized as proximity of the points. Samples that are ‘similar’ will be placed close to each other. 
<br><br>Users can customise the PCA/tSNE plot in the plot settings, including the {color} and {shape} of points using a phenotype class, choose t-SNE or PCA layout, label the points, or display 2D and 3D visualisation of the PCA/tSNE plot.'))


    require(plotly)
    require(scatterD3)
    
    observe({
        ngs <- inputData()
        req(ngs)
        ##input$menuitem  ## upon menuitem change
        var.types = colnames(ngs$Y)
        var.types = var.types[grep("sample|patient",var.types,invert=TRUE)]
        vv = c(var.types,rep("<none>",10))
        var.types0 = c("<none>","<cluster>",var.types)
        var.types0 = c("<none>",var.types)
        var.types1 = c("<none>",var.types)
        grp = vv[1]
        if("group" %in% var.types) grp = "group"
        updateSelectInput(session, "hmpca.colvar", choices=var.types0, selected=grp)
        updateSelectInput(session, "hmpca.shapevar", choices=var.types1, selected="<none>")
        ##updateSelectInput(session, "hmpca.line", choices=var.types1, selected="<none>")
        ##updateSelectInput(session, "hmpca.text", choices=var.types0, selected="group")
    })

    hm_getClusterPositions <- reactive({

        dbg("[hm_getClusterPositions] reacted")

        ngs <- inputData()
        req(ngs)

        ## take full matrix
        zx <- getFilteredMatrix()

        if(0) {
            zx <- as.matrix(ngs$X)
            kk <- selectSamplesFromSelectedLevels(ngs$Y, input_hm_samplefilter() )
            zx <- zx[,kk,drop=FALSE]
            ## --------------- apply feature filter
            gg = ngs$families[[1]]
            if(input$hm_features != "<all>") {
                gg = ngs$families[[input$hm_features]]
                jj <- match(toupper(gg), toupper(ngs$genes$gene_name))
                gg = ngs$genes$gene_name[setdiff(jj,NA)]
                pp = filterProbes(ngs$genes, gg)
                zx = zx[intersect(pp,rownames(zx)),]
            }
        }
        
        ntop = 1000
        ntop = as.integer(input$hm_ntop2)    
        zx = zx[order(-apply(zx,1,sd)),,drop=FALSE]  ## OK?
        if(nrow(zx) > ntop) {
            ##zx = head(zx,ntop)  ## OK?
            zx = zx[1:ntop,,drop=FALSE]  ## OK?
        }
        
        if("normalize" %in% input$hmpca_options) {
            zx <- scale(t(scale(t(zx))))
        }
        pdim = 2
        do3d <- ("3D" %in% input$hmpca_options)
        pdim = c(2,3)[ 1 + 1*do3d]
        
        pos = NULL
        force.compute = FALSE
        clustmethod = input$hm_clustmethod
        if(clustmethod=="fixed" && !force.compute) {
            if(pdim==2 && !is.null(ngs$tsne2d) ) {
                pos <- ngs$tsne2d[colnames(zx),]
            } else if(pdim==3 && !is.null(ngs$tsne3d) ) {
                pos <- ngs$tsne3d[colnames(zx),]
            }
        } else if("cluster" %in% names(ngs) &&
                  clustmethod %in% names(ngs$cluster$pos) )  {
            showNotification(paste("switching to ",clustmethod," layout...\n"))
            pos <- ngs$cluster$pos[[clustmethod]]
            if(pdim==2 && clustmethod=="pca") pos <- pos[,1:2]
            if(pdim==3 && clustmethod=="pca") pos <- pos[,1:3]

        } else  {
            showNotification(paste("computing ",clustmethod,"...\n"))
            
            perplexity = max(1,min((ncol(zx)-1)/3, 30))	
            perplexity

            res <- pgx.clusterSamplesFromMatrix(
                zx, perplexity=perplexity, is.logx=TRUE,
                ntop=999999, sv.rank=-1, prefix="C",         
                kclust=1, prior.count=NULL, 
                dims=pdim, find.clusters=FALSE,
                row.center=TRUE, row.scale=FALSE,
                method=clustmethod)
            if(pdim==2) pos <- res$pos2d
            if(pdim==3) pos <- res$pos3d
        }

        pos = scale(pos) ## scale
        ##colnames(pos) = paste0("dim",1:ncol(pos))
        ##rownames(pos) = colnames(zx)

        idx <- NULL
        if(FALSE) {
            ## should be relatively fast in sample space...
            cat("[hm_getClusterPositions] louvain clustering...\n")
            require(igraph)
            dist = as.dist(dist(pos))
            gr = graph_from_adjacency_matrix(
                1.0/dist, diag=FALSE, mode="undirected")
            idx <- cluster_louvain(gr)$membership
        }
        dbg("[hm_getClusterPositions] done")

        clust = list(pos=pos, clust=idx) 
        return(clust)
    })
    
    hm_PCAplot.RENDER <- reactive({

        ngs <- inputData()
        req(ngs)
        do3d = ("3D" %in% input$hmpca_options)

        dbg("[hm_PCAplot.RENDER:reactive] reacted")
        
        clust <- hm_getClusterPositions()

        dbg("[hm_PCAplot.RENDER:reactive] 1")
        
        pos <- clust$pos
        sel <- rownames(pos)
        df <- cbind(pos, ngs$Y[sel,])
        if(!is.null(clust$clust)) df[["<cluster>"]] <- clust$clust

        dbg("[hm_PCAplot.RENDER:reactive] 2")
        
        colvar = shapevar = linevar = textvar = NULL
        if(input$hmpca.colvar %in% colnames(df)) colvar <- factor(df[,input$hmpca.colvar])
        if(input$hmpca.shapevar %in% colnames(df)) shapevar <- factor(df[,input$hmpca.shapevar])
        ##if(input$hmpca.line %in% colnames(df)) linevar = factor(df[,input$hmpca.line])
        ##if(input$hmpca.text %in% colnames(df)) textvar = factor(df[,input$hmpca.text])
        mode = "markers"
        ann.text = rep(" ",nrow(df))
        if(!do3d && "label" %in% input$hmpca_options) ann.text = rownames(df)
        if(!is.null(colvar)) {
            colvar = factor(colvar)
            textvar <- factor(df[,input$hmpca.colvar])
        }
        symbols = c('circle','square','star','triangle-up','triangle-down','pentagon',
                    'bowtie','hexagon', 'asterisk','hash','cross','triangle-left',
                    'triangle-right','+',c(15:0))

        require(plotly)
        tt.info <- paste('Sample:', rownames(df),'</br>Group:', df$group)
        cex1 = c(1.0,0.8,0.6)[1 + 1*(nrow(pos)>30) + 1*(nrow(pos)>200)]

        if(do3d ) {
            ## 3D plot
            j0 = 1:nrow(df)
            j1 = NULL
            if(!is.null(linevar)) {
                linevar = factor(linevar)
                j0 = which(linevar==levels(linevar)[1])
                j1 = which(linevar!=levels(linevar)[1])
            }
            plt <- plot_ly(df, mode=mode) %>%
                add_markers(x = df[j0,1], y = df[j0,2], z = df[j0,3], type="scatter3d",
                            color = colvar[j0], ## size = sizevar, sizes=c(80,140),
                            ##marker = list(size = 5*cex1),
                            marker = list(size=5*cex1, line=list(color="grey10", width=0.1)),
                            symbol = shapevar[j0], symbols=symbols,
                            text = tt.info[j0] ) %>%
                add_annotations(x = pos[,1], y = pos[,2], z = pos[,3],
                                text = ann.text,
                                ##xref = "x", yref = "y",
                                showarrow = FALSE)
            if(!is.null(j1) & length(j1)>0) {
                plt <- plt %>%  add_markers(
                                    x = df[j1,1], y = df[j1,2], z = df[j1,3], type="scatter3d",
                                    color = colvar[j1], ## size = sizevar, sizes=c(80,140),
                                    ##marker = list(size=5*cex1, line=list(color="grey10", width=2)),
                                    symbol = shapevar[j1], symbols=symbols,
                                    text=tt.info[j1])
            }
            ## add cluster annotation labels
            if(0 && length(unique(colvar))>1) {
                ## add cluster annotation labels
                grp.pos <- apply(pos,2,function(x) tapply(x,colvar,median))
                ##grp.pos <- matrix(grp.pos, ncol=3)
                cex2 <- ifelse(length(grp.pos)>20,0.8,1)
                plt <- plt %>% add_annotations(
                                   x = grp.pos[,1], y = grp.pos[,2], z = grp.pos[,3],
                                   text = rownames(grp.pos),
                                   font=list(size=24*cex2, color='#555'),
                                   showarrow = FALSE)
            }

        } else {
            ## 2D plot
            j0 = 1:nrow(df)
            j1 = NULL
            if(!is.null(linevar)) {
                linevar = factor(linevar)
                j0 = which(linevar==levels(linevar)[1])
                j1 = which(linevar!=levels(linevar)[1])
            }
            plt <- plot_ly(df, mode=mode) %>%
                add_markers(x = df[j0,1], y = df[j0,2], type="scatter",
                            color = colvar[j0], ## size = sizevar, sizes=c(80,140),
                            marker = list(size=16*cex1, line=list(color="grey20", width=0.6)),
                            symbol = shapevar[j0], symbols=symbols,
                            text = tt.info[j0] ) %>%
                add_annotations(x = pos[,1], y = pos[,2],
                                text = ann.text,
                                ##xref = "x", yref = "y",
                                showarrow = FALSE)

            ## add node labels
            if(!is.null(j1) & length(j1)>0 ) {
                plt <- plt %>%  add_markers(
                                    x = df[j1,1], y = df[j1,2], type="scatter",
                                    color = colvar[j1], ## size = sizevar, sizes=c(80,140),
                                    marker = list(size=16*cex1, line=list(color="grey20", width=1.8)),
                                    symbol = shapevar[j1], symbols=symbols,
                                    text=tt.info[j1])
            }

            ## add group/cluster annotation labels
            if(!is.null(textvar) && length(unique(textvar))>1) {
                grp.pos <- apply(pos,2,function(x) tapply(x,as.character(textvar),median))
                cex2 <- 1
                if(length(grp.pos)>20) cex2 <- 0.8
                if(length(grp.pos)>50) cex2 <- 0.6
                plt <- plt %>% add_annotations(
                                   x = grp.pos[,1], y = grp.pos[,2],
                                   text = paste0("<b>",rownames(grp.pos),"</b>"),
                                   font = list(size=24*cex2, color='#555'),
                                   showarrow = FALSE)
            }
            plt <- plt %>%
                layout(showlegend = FALSE)
            
        }
        title = paste0("<b>PCA</b>  (",nrow(pos)," samples)")
        if(input$hm_clustmethod=="tsne") title = paste0("<b>tSNE</b>  (",nrow(pos)," samples)")
        ## plt <- plt %>% layout(title=title) %>% 
        ##     config(displayModeBar = FALSE)
        plt <- plt %>% 
            config(displayModeBar = FALSE) %>%
            ##config(modeBarButtonsToRemove = all.plotly.buttons ) %>%
            config(displaylogo = FALSE) %>% 
            config(toImageButtonOptions = list(format='svg', height=800, width=800))
        ##print(plt)
        return(plt)
    })

    hm_PCAplot_opts = tagList(
        tipify( selectInput( ns("hmpca.colvar"), "Color/label:", choices=NULL, width='100%'),
               "Set colors/labels according to a given phenotype.",
               placement="right", options = list(container = "body")),
        tipify( selectInput( ns("hmpca.shapevar"), "Shape:", choices=NULL, width='100%'),
               "Set shapes according to a given phenotype.",
               placement="right", options = list(container = "body")),
        tipify( checkboxGroupInput( ns('hmpca_options'),"Other:",
                                   choices=c('label','3D','normalize'), inline=TRUE),
               "Normalize matrix before calculating distances.",
               placement="right", options = list(container = "body")),
      
        ## NEED RETINK!!!: Would like to (dynamically) hide this in
        ## BASIC version. How to do conditional panel for BASIC/PRO in
        ## non-reactive env????
        ##
        ## conditionalPanel("input.hm_features == '<custom>'", ns=ns,
        div( class="pro-feature",
            tipify( radioButtons( ns('hm_clustmethod'),"Layout:",c("fixed","tsne","pca","umap"),inline=TRUE),
                   "Choose the layout method for clustering to visualise.",
                   placement="right", options = list(container = "body")),
            tipify( radioButtons( ns('hm_ntop2'),"Ntop:",c(100,1000,4000,9999),inline=TRUE,selected=1000),
                   "Number of top genes for dimensionality reduction.",
                   placement="right", options = list(container = "body"))
            )
    )

    hm_PCAplot_caption <- reactive({
        text1 = "The plot visualizes the similarity of samples as a scatterplot in reduced dimension (2D or 3D). Samples that are similar (in expression) are clustered near to each other, while samples with different expression are positioned farther away. Groups of samples with similar profiles will appear as <i>clusters</i> in the plot."
        if(input$hmpca.colvar!="<none>") {
            text1 <- paste(text1, "Colors correspond to the <strong>",input$hmpca.colvar,"</strong>phenotype.")
        }
        if(input$hmpca.shapevar!="<none>") {
            text1 <- paste(text1, "Shapes correspond to the <strong>",input$hmpca.shapevar,"</strong>phenotype.")
        }
        return(HTML(text1))
    })

    pca_caption_static = "<b>PCA/tSNE plot.</b> The plot visualizes the similarity in expression of samples as a scatterplot in reduced dimension (2D or 3D). Samples that are similar are clustered near to each other, while samples with different expression are positioned farther away. Groups of samples with similar profiles will appear as <i>clusters</i> in the plot."

    ##hm_PCAplot_module <- plotModule(
    callModule(
        plotModule, 
        id = "hm_PCAplot",
        func = hm_PCAplot.RENDER, ## ns=ns,
        plotlib = "plotly", 
        options = hm_PCAplot_opts,
        height = c(fullH-100,700), width=c("auto",800),
        pdf.width=8, pdf.height=8,
        title="PCA/tSNE plot",
        info.text = hm_PCAplot_text
        ##caption = pca_caption_static
    )

    output$hm_pcaUI <- renderUI({
        fillCol(
            flex = c(NA,0.02,1),
            height = fullH,
            div(HTML(pca_caption_static), class="caption"),
            br(),
            plotWidget(ns("hm_PCAplot"))
        )
    })
    outputOptions(output, "hm_pcaUI", suspendWhenHidden=FALSE) ## important!!!
    
    
    ##================================================================================
    ## Parallel coordinates
    ##================================================================================

    hm_parcoord.ranges <- reactiveValues()
    
    hm_parcoord.matrix <- reactive({
        ngs <- inputData()
        filt <- getTopMatrix()
        req(filt)
        ##zx <- filt$grp.mat[,]
        zx <- filt$mat[,]    
        grp <- input$hm_pcvar
        if(!(grp %in% colnames(ngs$Y))) grp <- "group"
        y <- ngs$Y[colnames(zx),grp]
        zx <- t(apply(zx,1,function(x) tapply(x,y,mean,na.rm=TRUE)))   
        if(input$hm_pcscale) { zx <- t(scale(t(zx))) }
        
        rr <- isolate(reactiveValuesToList(hm_parcoord.ranges))
        nrange <- length(rr)
        for(i in names(rr)) hm_parcoord.ranges[[i]] <- NULL
        zx <- round(zx, digits=3)
        list(mat=zx, clust=filt$clust)
    })
    
    hm_parcoord.RENDER <- reactive({
        
        pc <- hm_parcoord.matrix()
        req(pc)
        zx <- pc$mat
        ## build dimensions
        dimensions <- list()
        for(i in 1:ncol(zx)) {
            dimensions[[i]] <-  list(
                range = c(min(zx[,i]),max(zx[,i])),
                ## constraintrange = c(100000,150000),
                ## tickvals = c(0,0.5,1,2,3),
                ## ticktext = c('A','AB','B','Y','Z'),
                visible = TRUE,
                label = colnames(zx)[i],
                values = zx[,i]
            )
        }
        
        clust.id <- as.integer(factor(pc$clust))
        table(clust.id)
        
        df <- data.frame(clust.id=clust.id, zx)
        klrpal = rep(brewer.pal(8,"Set2"),99)
        ##klrpal = rep(c("red","blue","green","yellow","magenta","cyan","black","grey"),99)
        klrpal = klrpal[1:max(clust.id)]
        ##klrpal <- setNames(klrpal, sort(unique(clust.id)))
        klrpal2 <- lapply(1:length(klrpal),function(i) c((i-1)/(length(klrpal)-1),klrpal[i]))
        
        plt <-  plot_ly(df, source = "pcoords") %>%    
            add_trace(type = 'parcoords', 
                      line = list(color = ~clust.id,
                                  ## colorscale = list(c(0,'red'),c(0.5,'green'),c(1,'blue'))
                                  ##colorscale = 'Jet',
                                  colorscale = klrpal2,
                                  cmin = min(clust.id), cmax = max(clust.id),
                                  showscale = FALSE
                                  ##reversescale = TRUE
                                  ),
                      dimensions = dimensions)
        plt <- plt %>%
            layout(margin = list(l=60, r=60, t=0, b=30)) %>%
            ##config(displayModeBar = FALSE) %>%
            ##config(modeBarButtonsToRemove = setdiff(all.plotly.buttons,"toImage") ) %>%
            config(toImageButtonOptions = list(format='svg', width=900, height=350, scale=1.2)) %>%
            config(displaylogo = FALSE) %>%
            event_register("plotly_restyle")

        plt 

    })

    observeEvent( event_data("plotly_restyle", source = "pcoords"), {
        ## From: https://rdrr.io/cran/plotly/src/inst/examples/shiny/event_data_parcoords/app.R
        ##
        d <- event_data("plotly_restyle", source = "pcoords")
        ## what is the relevant dimension (i.e. variable)?
        dimension <- as.numeric(stringr::str_extract(names(d[[1]]), "[0-9]+"))
        ## If the restyle isn't related to a dimension, exit early.
        if (!length(dimension)) return()
        if (is.na(dimension)) return()
        
        pc <- hm_parcoord.matrix()
        req(pc)
        ## careful of the indexing in JS (0) versus R (1)!
        dimension_name <- colnames(pc$mat)[[dimension + 1]]
        ## a given dimension can have multiple selected ranges
        ## these will come in as 3D arrays, but a list of vectors 
        ## is nicer to work with
        info <- d[[1]][[1]]    
        if (length(dim(info)) == 3) {
            hm_parcoord.ranges[[dimension_name]] <- lapply(seq_len(dim(info)[2]), function(i) info[,i,])
        } else {
            hm_parcoord.ranges[[dimension_name]] <- list(as.numeric(info))
        }
    })

    hm_parcoord.selected <- reactive({

        mat <- hm_parcoord.matrix()$mat
        clust <- hm_parcoord.matrix()$clust
        req(mat)
        keep <- TRUE
        for (i in names(hm_parcoord.ranges)) {
            range_ <- hm_parcoord.ranges[[i]]
            range_ <- range_[sapply(range_,length)>0]
            if(length(range_)>0) {
                keep_var <- FALSE
                for (j in seq_along(range_)) {
                    rng <- range_[[j]]
                    keep_var <- keep_var | dplyr::between(mat[,i], min(rng), max(rng))
                }
                keep <- keep & keep_var
            }
        }
        list(mat=mat[keep,,drop=FALSE], clust=clust[keep])
    })


    hm_parcoord_opts = tagList(
        tipify( selectInput(ns("hm_pcvar"), "Variable:", choices="group", width='100%'),
               "Split the parallel axis by this phenotype.",
               placement="right",options = list(container = "body")),
        tipify( checkboxInput(ns('hm_pcscale'),'Scale values',TRUE),
               "Scale expression values to mean=0 and SD=1.",
               placement="right",options = list(container = "body"))
    )


    hm_parcoord_text = tagsub("The <strong>Parallel Coordinates</strong> panel
displays the expression levels of selected genes across all conditions in the analysis. On the x-axis the experimental conditions are plotted. The y-axis shows the expression level of the genes grouped by condition. The colors correspond to the gene groups as defined by the hierarchical clustered heatmap.")

    callModule(
        plotModule,     
        ## hm_parcoord_module <- plotModule(
        "hm_parcoord",
        func = hm_parcoord.RENDER, ## ns = ns,
        plotlib="plotly", ## renderFunc="renderPlotly",
        ## download.fmt = c("png","pdf","html"),  ## PNG & PDF do not work!!! 
        ## download.fmt = c("html"),
        options = hm_parcoord_opts,
        height = c(0.45*fullH,600), width = c("100%",1000),
        pdf.width=10, pdf.height=6, info.width="350px",
        title="Parallel coordinates", label="a",
        info.text = hm_parcoord_text
        ## caption = hm_parcoord_text,
    )

    hm_parcoord_table.RENDER <- reactive({

        mat = hm_parcoord.selected()$mat
        clust = hm_parcoord.selected()$clust
        df <- data.frame(cluster=clust, mat, check.names=FALSE)
        numeric.cols <- 2:ncol(df)
        DT::datatable(
                df, rownames=TRUE, ## escape = c(-1,-2),
                extensions = c('Buttons','Scroller'),
                selection=list(mode='single', target='row', selected=NULL),
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
            formatSignif(numeric.cols,3) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') 
    })

    hm_parcoord_table_info = "In this table, users can check mean expression values of features across the conditions for the selected genes."

    parcoord_caption = "<b>Parallel Coordinates plot.</b> <b>(a)</b>The Parallel Coordinates plot displays the expression levels of selected genes across all conditions. On the x-axis the experimental conditions are plotted. The y-axis shows the expression level of the genes grouped by condition. The colors correspond to the gene groups as defined by the hierarchical clustered heatmap. <b>(b)</b> Average expression of selected genes across conditions."
    
    ##hm_parcoord_table_module <- tableModule(
    hm_parcoord_table_module <- callModule(
        tableModule, id = "hm_parcoord_table",
        func = hm_parcoord_table.RENDER, ## ns=ns,
        info.text = hm_parcoord_table_info,
        title = "Selected genes", label="b",
        height = c(270,700)
        ## caption = parcoord_caption
    )

    output$hm_parcoordUI <- renderUI({
        fillCol(
            flex = c(NA,0.02,1.2,0.05,1),
            height = fullH,
            div(HTML(parcoord_caption), class="caption"),
            br(),
            plotWidget(ns("hm_parcoord")),
            br(),
            tableWidget(ns("hm_parcoord_table"))
        )
    })
    outputOptions(output, "hm_parcoordUI", suspendWhenHidden=FALSE) ## important!!!

    ##================================================================================
    ## Annotate clusters
    ##================================================================================

    clustannot_plots_text = paste0('The top features of the heatmap in the <code>Heatmap</code> panel are divided into gene (or gene set) clusters based on their expression profile patterns. For each cluster, the platform provides a functional annotation in the <code>Annotate cluster</code> panel by correlating annotation features from more than 42 published reference databases, including well-known databases such as ',a_MSigDB,', ',a_KEGG,' and ',a_GO,'. In the plot settings, users can specify the level and reference set to be used under the <code>Reference level</code> and <code>Reference set</code> settings, respectively.')
    
    observe({
        ngs <- inputData()    
        if(is.null(input$xann_level)) return(NULL)
        ann.types=sel=NULL
        if(input$xann_level!="phenotype") {
            if(input$xann_level=="geneset") {
                ann.types <- names(COLLECTIONS)
                cc = sapply(COLLECTIONS,function(s) length(intersect(s,rownames(ngs$gsetX))))
                ann.types <- ann.types[cc>=3]
            }
            if(input$xann_level=="gene") {
                ann.types <- names(ngs$families)
                cc = sapply(ngs$families,function(g) length(intersect(g,rownames(ngs$X))))
                ann.types <- ann.types[cc>=3]
            }
            ann.types <- setdiff(ann.types,"<all>")  ## avoid slow...
            ann.types <- grep("^<",ann.types,invert=TRUE,value=TRUE)  ## remove special groups
            sel = ann.types[1]
            if("H" %in% ann.types) sel = "H"
            j <- grep("^transcription",ann.types,ignore.case=TRUE)
            if(input$xann_level=="geneset") j <- grep("hallmark",ann.types,ignore.case=TRUE)
            if(length(j)>0) sel = ann.types[j[1]]
            ann.types <- sort(ann.types)
        } else {
            ann.types = sel = "<all>"
        }
        updateSelectInput(session, "xann_refset", choices=ann.types, selected=sel)    
    })


    getClustannotCorrelation <- reactive({
        
        ngs <- inputData()
        req(ngs)
        filt <- getTopMatrix()
        req(filt)
        
        zx  <- filt$mat
        idx <- filt$clust
        samples <- filt$samples

        ann.level="geneset"
        ann.refset="Hallmark collection"
        ann.level = input$xann_level
        ##if(is.null(ann.level)) return(NULL)
        ann.refset = input$xann_refset
        ##if(is.null(ann.refset)) return(NULL)
        req(input$xann_level, input$xann_refset)
        
        ref = NULL
        ref = ngs$gsetX[,]    
        ref = ngs$X[,]    
        if(ann.level=="gene" && ann.refset %in% names(ngs$families) ) {
            gg = ngs$families[[ann.refset]]
            jj = match(toupper(gg), toupper(ngs$genes$gene_name))
            jj <- setdiff(jj,NA)
            pp = rownames(ngs$genes)[jj]
            ref = ngs$X[intersect(pp,rownames(ngs$X)),]    
        }
        if(ann.level=="geneset" && ann.refset %in% names(COLLECTIONS)) {
            ss = COLLECTIONS[[ann.refset]]
            ss = intersect(ss, rownames(ngs$gsetX))
            length(ss)
            ref = ngs$gsetX[ss,]    
        }
        if(ann.level=="phenotype") {
            ref = t(expandAnnotationMatrix(ngs$Y))
        }
        if(is.null(ref)) {
            cat("<clustering:getClustannotCorrelation> WARNING:: ref error\n")
            return(NULL)
        }
        
        ##-----------  restrict to top??
        dim(ref)
        if(nrow(ref)>1000) {
            ref = head(ref[order(-apply(ref,1,sd)),],1000)
        }

        ##-----------  get original data level
        X = ngs$X
        if(input$hm_level=="geneset") X <- ngs$gsetX
        
        ##----------- for each gene cluster compute average correlation
        hm_topmode = "sd"
        hm_topmode <- input$hm_topmode
        idxx = setdiff(idx, c(NA," ","   "))    
        rho <- matrix(NA, nrow(ref), length(idxx))
        colnames(rho) <- idxx
        rownames(rho) <- rownames(ref)    
        i=1
        if(nrow(ref)>0) {
            for(i in 1:length(idxx)) {
                gg = rownames(zx)[which(idx==idxx[i])]        
                aa <- t(X[gg,samples,drop=FALSE])
                bb <- t(ref[,samples,drop=FALSE])
                ##rr = cor(aa , bb, use="pairwise", method="spearman")
                rr = cor(apply(aa,2,rank), apply(bb,2,rank), use="pairwise")
                if(hm_topmode=="pca") rr <- abs(rr)
                rho[,i] <- colMeans(rr,na.rm=TRUE)
            }
        }

        if(input$hm_level=="gene" && ann.level=="geneset" && input$xann_odds_weighting ) {
            table(idx)
            grp <- tapply( toupper(rownames(zx)), idx, list)  ## toupper for mouse!!
            gmt <- GSETS[rownames(rho)]
            bg.genes <- toupper(rownames(X))
            dbg("[getClustannotCorrelation] bg.genes=",head(bg.genes))
            dbg("[getClustannotCorrelation] head(gmt[[1]])=",head(gmt[[1]]))

            P <- c()
            for(i in 1:ncol(rho)) {
                k <- colnames(rho)[i]
                res <- gset.fisher(
                    grp[[k]], gmt, fdr=1, min.genes=0, max.genes=Inf,
                    background = bg.genes )
                res <- res[rownames(rho),]
                r <- res[,"odd.ratio"]
                odd.prob <- r / (1+r)
                ##odd.1mpv <- 1 - res[,"p.value"]
                ##P <- cbind(P,odd.1mpv)
                P <- cbind(P,odd.prob)
            }
            colnames(P) <- colnames(rho)
            rownames(P) <- rownames(rho)
            rho <- rho * (P/max(P))
        }

        ##rho = round(rho, digits=3)
        dim(rho)
        return(rho)
    })

    clustannot_plots.PLOTLY <- reactive({
        require(RColorBrewer)
        rho = getClustannotCorrelation()
        ##if(is.null(rho)) return(NULL)
        req(rho)
        
        ##par(mfrow=c(2,3), mar=c(3.5,2,2,1), mgp=c(2,0.8,0))
        NTERMS = 6
        NTERMS = 12
        slen=40
        if(ncol(rho)>=5)  {
            slen=20
        }
        if(ncol(rho)>6)  {
            NTERMS=6
        }
        if(ncol(rho)<=2) {
            NTERMS=22
        }

        klrpal = rep(brewer.pal(8,"Set2"),2)
        ##klrpal = paste0(klrpal,"88")
        col.addalpha <- function(clr,a=100)
            paste0("rgba(",paste(col2rgb(clr)[,1],collapse=","),",",a,")")
        ##klrpal = as.character(sapply(klrpal, col.addalpha, a=50))
        klrpal <- paste0(klrpal,"55")
        
        plot_list <- list()
        i=1    
        for(i in 1:min(9,ncol(rho))) {
            
            x = rev(head(sort(rho[,i],decreasing=TRUE),NTERMS))
            names(x) = sub(".*:","",names(x))
            names(x) = gsub(GSET.PREFIX.REGEX,"",names(x))

            y = names(x)
            y = factor(y, levels=y)
            anntitle <- function(tt) {
                list(text=tt, font = list(size=13),
                     xref="paper", yref="paper",
                     yanchor = "bottom", xanchor = "center",
                     align = "center", x=0.5, y=1.02 , showarrow = FALSE )
            }
            
            plot_list[[i]] <- plot_ly(
                x=x, y=y, type='bar',  orientation='h',
                text=y, hoverinfo = 'text',
                hovertemplate = paste0("%{y}<extra>",colnames(rho)[i],"</extra>"),
                ##hovertemplate = "%{y}",
                marker = list(color=klrpal[i])) %>%
                layout(
                    showlegend = FALSE,
                    annotations = anntitle(colnames(rho)[i]),
                    ## annotations = list(text="TITLE"),
                    ## margin = c(0, 0.0, 0.05, 0.05),
                    margin = list(l=5, r=0, t=25, b=15),
                    xaxis = list(range = c(0,0.9),
                                 titlefont = list(size=11),
                                 tickfont = list(size=10),
                                 showgrid=FALSE,
                                 title = "\ncorrelation (R)" ),
                    yaxis = list(title = "",
                                 showgrid = FALSE,
                                 showline = FALSE,
                                 showticklabels = FALSE,
                                 showgrid=FALSE,
                                 zeroline = FALSE)
                ) %>%
                ## labeling the y-axis inside bars
                add_annotations(xref = 'paper', yref = 'y',
                                x = 0.01, y = y, xanchor='left',
                                text = shortstring(y,slen), 
                                font = list(size = 10),
                                showarrow = FALSE, align = 'right')
            ##layout(margin = c(0,0,0,0))
        }

        if(length(plot_list) <= 4) {
            nrows = ceiling(length(plot_list)/2 )
        } else {
            nrows = ceiling(length(plot_list)/3 )
        }
        
        subplot( plot_list, nrows=nrows, shareX=TRUE,
                ## template = "plotly_dark",
                margin = c(0, 0.0, 0.05, 0.05) ) %>%
            config(displayModeBar = FALSE) 
    })

    clustannot_plots_opts = tagList(
        tipify( selectInput(ns("xann_level"), "Reference level:",
                            choices=c("gene","geneset","phenotype"),
                            selected="geneset", width='80%'),
               "Select the level of an anotation analysis.",
               placement="left", options = list(container = "body")),
        conditionalPanel(
            "input.xann_level == 'geneset'", ns=ns,
            tipify( checkboxInput(ns("xann_odds_weighting"), "Fisher test weighting"),
                   "Enable weighting with Fisher test probability for gene sets. This will effectively penalize small clusters and increase robustness.",
                   placement="right", options = list(container = "body"))
        ),
        tipify( selectInput( ns("xann_refset"), "Reference set:", choices="", width='80%'),
               "Specify a reference set to be used in the annotation.",
               placement="left",options = list(container = "body"))
    )
    
    ##clustannot_plots_module <- plotModule(
    callModule(
        plotModule, 
        id="clustannot_plots", ##ns=ns,
        ##func=clustannot_plots.RENDER, plotlib = "base",
        func = clustannot_plots.PLOTLY, plotlib="plotly",
        download.fmt = c("png","pdf"),
        options = clustannot_plots_opts,
        height = c(360,600), width = c(500,1000),
        pdf.width=8, pdf.height=5, res=80,
        title="Functional annotation of clusters", label="a",
        info.text = clustannot_plots_text        
    )
    
    clustannot_table.RENDER <- reactive({
        
        rho = getClustannotCorrelation()
        if(is.null(rho)) return(NULL)
        
        ##rownames(rho) = shortstring(rownames(rho),50)
        rho.name = shortstring(sub(".*:","",rownames(rho)),60)
        ##rho = data.frame(cbind( name=rho.name, rho))
        df = data.frame( feature=rho.name, round(as.matrix(rho),digits=3))
        rownames(df) = rownames(rho)
        if(input$xann_level=="geneset") {
            df$feature <- wrapHyperLink(df$feature, rownames(df))
        }

        DT::datatable(
                df, rownames=FALSE, escape = c(-1,-2),
                extensions = c('Buttons','Scroller'),
                selection=list(mode='single', target='row', selected=c(1)),
                class = 'compact cell-border stripe hover',
                fillContainer = TRUE,
                options=list(
                    dom = 'lfrtip', buttons = c('copy','csv','pdf'),
                    ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                    scrollX = TRUE, ##scrollY = TRUE,
                    ##scrollY = 170,
                    scrollY = '70vh',
                    scroller=TRUE,
                    deferRender=TRUE
                )  ## end of options.list 
            ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') 
    })

    clustannot_table_info_text = "In this table, users can check mean correlation values of features in the clusters with respect to the annotation references database selected in the settings."
    
    ##clustannot_table_module <- tableModule(
    clustannot_table_module <- callModule(
        tableModule, 
        id = "clustannot_table", 
        func = clustannot_table.RENDER,
        ##options = clustannot_table_opts,
        info.text = clustannot_table_info_text,
        title="Annotation scores", label="b",
        height = c(240,700), width=c('auto',1000),
        ##caption = clustannot_caption
    )

    clustannot_caption = "<b>Cluster annotation.</b> <b>(a)</b> Top ranked annotation features (by correlation) for each gene cluster as defined  in the heatmap. <b>(b)</b> Table of average correlation values of annotation features, for each gene cluster."
    
    output$hm_annotateUI <- renderUI({
        fillCol(
            flex = c(NA,0.03,1.4,0.03,1),
            height = fullH,
            div(HTML(clustannot_caption), class="caption"),
            br(),
            plotWidget(ns("clustannot_plots")),
            br(),
            plotWidget(ns("clustannot_table"))
        )
    })
    outputOptions(output, "hm_annotateUI", suspendWhenHidden=FALSE) ## important!!!
    
    ##================================================================================
    ## Phenotypes {data-height=800}
    ##================================================================================

    clust_phenoplot.RENDER %<a-% reactive({
        ##if(!input$tsne.all) return(NULL)
        require(RColorBrewer)
        
        ngs <- inputData()
        req(ngs)
        
        ## get t-SNE positions
        clust <- hm_getClusterPositions()
        ##pos = ngs$tsne2d
        pos = clust$pos
        Y <- ngs$Y[rownames(pos),,drop=FALSE]
        pheno = colnames(Y)    

        ## don't show these...
        pheno <- grep("batch|sample|donor|repl|surv",pheno,
                      invert=TRUE, ignore.case=TRUE,value=TRUE)
        
        ## layout
        par(mfrow = c(3,2), mar=c(0.3,0.7,2.8,0.7))
        if(length(pheno)>=6) par(mfrow = c(4,3), mar=c(0.3,0.4,2.8,0.4)*0.8)
        if(length(pheno)>=12) par(mfrow = c(5,4), mar=c(0.2,0.2,2.5,0.2)*0.8)
        i=1    

        cex1 <- 1.3*c(1.8,1.3,0.8,0.5)[cut(nrow(pos),breaks=c(-1,40,200,1000,1e10))]    
        cex1 = cex1 * ifelse(length(pheno)>6, 0.8, 1)
        cex1 = cex1 * ifelse(length(pheno)>12, 0.8, 1)
        
        require(RColorBrewer)
        for(i in 1:min(20,length(pheno))) {

            ## ------- set colors
            colvar = factor(Y[,1])
            colvar = factor(Y[,pheno[i]])
            colvar[which(colvar %in% c(NA,""," ","NA","na"))] <- NA
            colvar = factor(as.character(colvar))
            klrpal = COLORS  
            klr1 = klrpal[colvar]
            klr1 = paste0(col2hex(klr1),"99")
            jj = which(is.na(klr1))
            if(length(jj)) klr1[jj] <- "#AAAAAA22"
            tt = tolower(pheno[i])

            ## ------- start plot
            plot( pos[,], pch=19, cex=cex1, col=klr1,
                 fg = gray(0.5), bty = "o", xaxt='n', yaxt='n',
                 xlab="tSNE1", ylab="tSNE2")
            title( tt, cex.main=1.3, line=0.5, col="grey40")
            if(input$clust_phenoplot_labelmode=="legend") {
                legend("bottomright", legend=levels(colvar), fill=klrpal,
                       cex=0.95, y.intersp=0.85, bg="white")
            } else {
                grp.pos <- apply(pos,2,function(x) tapply(x,colvar,mean,na.rm=TRUE))
                grp.pos <- apply(pos,2,function(x) tapply(x,colvar,median,na.rm=TRUE))
                nvar <- length(setdiff(colvar,NA))
                if(nvar==1) {
                    grp.pos <- matrix(grp.pos,nrow=1)
                    rownames(grp.pos) <- setdiff(colvar,NA)[1]
                }
                labels = rownames(grp.pos)
                boxes = sapply(nchar(labels),function(n) paste(rep("\u2588",n),collapse=""))
                cex2 = 0.99*cex1**0.33
                text( grp.pos, labels=boxes, cex=cex2*0.95, col="#CCCCCC99")
                text( grp.pos, labels=labels, font=2, cex=cex2)
            }
        }    
    })

    clust_phenoplot.opts = tagList(
        radioButtons(ns('clust_phenoplot_labelmode'),"Label",c("groups","legend"),inline=TRUE)
    )

    clust_phenoplot_info = tagsub("<strong>Phenotype distribution.</strong> This figure visualizes the distribution of the available phenotype data. You can choose to put the group labels in the figure or as separate legend in the {Label} setting, in the plot {{settings}}")
    
    clust_phenoplot_caption = "<b>Phenotype distribution.</b> The plots show the distribution of the phenotypes superposed on the t-SNE clustering. Often, we can expect the t-SNE distribution to be driven by the particular phenotype that is controlled by the experimental condition or unwanted batch effects."

    ## clust_phenoplot.module <- plotModule(
    callModule(
        plotModule,         
        "clust_phenoplot", ## ns=ns,
        func  = clust_phenoplot.RENDER, ## plotlib="base",
        func2 = clust_phenoplot.RENDER, ## plotlib="base",
        options = clust_phenoplot.opts,
        height = c(fullH-100,700), res=85,
        pdf.width=6, pdf.height=9, 
        info.text = clust_phenoplot_info
        ## caption = clust_phenoplot_caption
    )

    output$hm_phenoplotUI <- renderUI({
        fillCol(
            flex = c(NA,0.025,1),
            height = fullH,
            div(HTML(clust_phenoplot_caption),class="caption"),
            br(),
            plotWidget(ns("clust_phenoplot"))
        )
    })
    ## outputOptions(output, "hm_phenoplotUI", suspendWhenHidden=FALSE) ## important!!!
    
    ##=============================================================================
    ## Feature ranking
    ##=============================================================================
    
    require(plotly)
    
    calcFeatureRanking <- reactive({
        ngs <- inputData()
        req(ngs)

        features=X=NULL
        if(input$hm_level=="geneset") {
            features = COLLECTIONS
            X = ngs$gsetX
        } else {
            features = ngs$families
            X = ngs$X
        }

        ## ------------ intersect features, set minimum set size
        rownames(X) <- toupper(rownames(X))
        genes <- toupper(rownames(X))
        features <- lapply(features, toupper)
        features <- lapply(features, function(f) intersect(toupper(f), genes))
        features <- features[sapply(features,length) >=10 ]

        cat("[calcFeatureRanking] length(features)=",length(features),"\n")
        
        ## ------------ Just to get current samples
        ##samples = getTopMatrix()$samples    
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input_hm_samplefilter() )
        X = X[,samples]
        cvar <- pgx.getCategoricalPhenotypes(ngs$Y)
        cvar = grep("group|sample|patient|years|days|months|gender",
                    cvar,invert=TRUE,value=TRUE) ## no sample IDs
        cvar
        Y = ngs$Y[colnames(X),cvar,drop=FALSE]
        kk = which(apply(Y,2,function(y) length(unique(y))>1))
        Y = Y[,kk,drop=FALSE]
        dim(Y)

        cat("[calcFeatureRanking] dim(X)=",dim(X),"\n")
        cat("[calcFeatureRanking] dim(Y)=",dim(Y),"\n")
        
        ## ------------ Note: this takes a while. Maybe better precompute off-line...
        sdx = apply(X,1,sd)
        names(sdx) = rownames(X)
        S = matrix(NA, nrow=length(features), ncol=ncol(Y))
        rownames(S) = names(features)
        colnames(S) = colnames(Y)

        ## ------------ Create a Progress object
        progress <- shiny::Progress$new()
        on.exit(progress$close())    
        progress$set(message = "Calculating feature-set scores", value = 0)
        
        gene.level = TRUE
        gene.level = (input$hm_level=="gene")
        i=1
        for(i in 1:ncol(Y)) {

            progress$inc(1/ncol(Y))

            grp = Y[,i]
            grp = as.character(grp)
            
            cat("[calcFeatureRanking] head(grp)=",head(grp),"\n")
            
            score = rep(NA, length(features))
            names(score) = names(features)
            j=1
            for(j in 1:length(features)) {

                pp = features[[j]]
                
                if(gene.level) {
                    pp = filterProbes(ngs$genes, features[[j]])
                }
                pp = head(pp[order(-sdx[pp])],1000)  ## how many top SD??
                pp = intersect(pp, rownames(X))
                X1 = X[pp,,drop=FALSE]
                dim(X1)
                ##cat("<clust_featureRank> dim(X1)=",dim(X1),"\n")
                ##if( nrow(X1) 
                
                s1 = s2 = 1
                method = input$clust_featureRank_method
                if(method %in% c("correlation","meta")) {
                    mx = t(apply(X1, 1, function(x) tapply(x,grp,mean)))
                    if(nrow(mx)==0 || ncol(mx)==0) next
                    D = 1 - cor(mx, use="pairwise")
                    diag(D) = NA
                    s1 = mean(D,na.rm=TRUE)
                }
                
                if(method %in% c("p-value","meta")) {
                    jj  <- which(!is.na(grp))
                    design = model.matrix( ~ grp[jj])                    
                    suppressWarnings( fit <- eBayes(lmFit( X1[,jj], design)) )
                    suppressWarnings( suppressMessages( top <- topTable(fit) ))
                    ##s2 = mean(-log10(top$P.Value))  ## as score
                    s2 = mean(-log10(top$adj.P.Val),na.rm=TRUE)  ## as score
                }
                
                f = 1
                f <- (1 - exp(-(length(pp)/20)**2)) ## penalize smaller sets
                score[j] = f * (s1 * s2) ** ifelse(method=="meta",0.5,1)
                
            }
            S[,i] = score
        }
        return(S)
    })
    
    clust_featureRank.RENDER %<a-% reactive({

        S <- calcFeatureRanking()
        
        if(is.null(S) || nrow(S)==0 || ncol(S)==0 ) return(NULL)

        ## top scoring
        S = tail( S[order(rowSums(S)),,drop=FALSE], 35)  

        par(mfrow=c(2,1), mar=c(1,5,3,3) )
        par(mfrow=c(1,2), mar=c(5,5,3,2), oma=c(6,0,3,0)); frame()
        ## par(mfrow=c(1,1), mar=c(10,5,3,3) )
        rownames(S) = substring(rownames(S),1,80)
        bpos = barplot( t(S), beside=FALSE, las=1,
                       cex.names=0.9, horiz=TRUE,
                       xlab="discriminant score" )
        ##title("feature-set score", cex=1.3)
        cc1 = grey.colors(ncol(S))
        legend("bottomright",legend=colnames(S), fill=cc1,
               cex=0.8, y.intersp=0.8, inset=c(0,0.035), bg="white")
        
    })


    clust_featureRank_info = "Ranked discriminant score for top feature sets. The plot ranks the discriminitive power of the feature set (genes) as a cumulative discriminant score for all phenotype variables. In this way, we can find which feature set (or gene family/set) can explain the variance in the data the best. <p>Correlation-based discriminative power is calculated as the average '(1-cor)' between the groups. Thus, a feature set is highly discriminative if the between-group correlation is low. P-value based scoring is computed as the average negative log p-value from the ANOVA. The 'meta' method combines the score of the former methods in a multiplicative manner."

    clust_featureRank_caption = "<b>Feature-set ranking.</b> Ranked discriminant score for top feature sets. The plot ranks the discriminative power of feature sets (or gene sets) as the cumulative discriminant score for all phenotype variables."

    clust_featureRank.opts =  tagList(
        tipify( radioButtons( ns('clust_featureRank_method'),'Method:',
                             choices=c("p-value","correlation","meta"),
                             inline=TRUE),
               "Choose ranking method: p-value based or correlation-based.",
               placement="right", options = list(container = "body") )
    )

##    clust_featureRank_module <- plotModule(
##        title="Feature-set ranking", ns=ns,
    callModule(
        plotModule, 
        id="clust_featureRank",
        title="Feature-set ranking", 
        func = clust_featureRank.RENDER,
        func2 = clust_featureRank.RENDER,
        options = clust_featureRank.opts,
        pdf.width=8, pdf.height=8,
        height = c(fullH-80,700), width=c("auto",800),
        res = c(72,90),
        info.text = clust_featureRank_info
        ## caption = clust_featureRank_caption
    )
    
    output$hm_featurerankUI <- renderUI({
        fillCol(
            flex = c(NA,0.035,1),
            height = fullH,
            div(HTML(clust_featureRank_caption),class="caption"),
            br(),
            plotWidget(ns("clust_featureRank"))
        )
    })


} ## end of Module
