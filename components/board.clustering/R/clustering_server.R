##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

ClusteringBoard <- function(id, pgx)
{
  moduleServer(id, function(input, output, session) {

    ns <- session$ns ## NAMESPACE
    fullH = 850  ## full height of page
    clust_infotext = paste('
The <strong>Clustering Analysis</strong> module performs unsupervised clustering analysis of the data. After having done the QC, it is probably the first way to explore your data. The main purpose is to discover patterns and subgroups in the data, show correlation with known phenotypes, detect outliers, or investigate batch effects.

<br><br>In the <strong>Heatmap</strong> panel hierarchical clustering can be performed on gene level or gene set level (selected under the {Level} dropdown list). During the heatmap generation, the platform provides functional annotation for each feature cluster in <strong>Annotate cluster</strong> panel. Users can select from a variety of annotation databases from the literature, such as ',a_MSigDB,', ',a_KEGG,' and ',a_GO,'. The <strong>PCA/tSNE</strong> panel shows unsupervised clustering of the samples in 2D/3D as obtained by ',a_PCA,' or ',a_tSNE,' algorithms. The <strong>Phenotypes</strong> panel on the right, shows the phenotype distribution as colors on the t-SNE plot. 

<br><br>EXPERT MODE ONLY: The <strong>Feature ranking</strong> panel computes a discriminant score for gene (or geneset) families. This allows to investigate what family of genes (or gene sets) can best discriminate the groups. 

<br><br><br>
<center><iframe width="560" height="315" src="https://www.youtube.com/embed/hyDEk_MCaTk" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>

') 
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    shiny::observeEvent( input$clust_info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Clustering Board</strong>"),
            shiny::HTML(clust_infotext),
            easyClose = TRUE, size="l" ))
    })
       
    ## update filter choices upon change of data set 
    shiny::observe({

        shiny::req(pgx$X)
        
        levels = getLevels(pgx$Y)
        shiny::updateSelectInput(session, "hm_samplefilter", choices=levels)
        
        if(DEV && !is.null(pgx$gset.meta$matrices) ) {
            jj = which(!sapply(pgx$gset.meta$matrices,is.null))
            mat.names = names(pgx$gset.meta$matrices)[jj]
            shiny::updateRadioButtons(session, "hm_gsetmatrix", choices=mat.names,
                               selected="meta", inline=TRUE)
        }
        
        shiny::updateRadioButtons(session, "hm_splitby", selected='none')

        ## update defaults??
        ##if(ncol(pgx$X) > 80) shiny::updateNumericInput(session,"hm_cexCol", value=0)

        ## update defaults??
        n1 <- nrow(pgx$samples)-1
        groupings <- colnames(pgx$samples)
        ## groupings <- pgx.getCategoricalPhenotypes(pgx$samples, min.ncat=2, max.ncat=n1)
        groupings <- c("<ungrouped>",sort(groupings))
        shiny::updateSelectInput(session,"hm_group", choices=groupings)
        contrasts <- pgx.getContrasts(pgx)
        shiny::updateSelectInput(session,"hm_contrast", choices=contrasts)                

    })

    shiny::observeEvent( input$hm_splitby, {

        shiny::req(pgx$X, pgx$samples)
        
        if(input$hm_splitby=='none') return()
        if(input$hm_splitby=='gene') {
            xgenes <- sort(rownames(pgx$X))
            shiny::updateSelectizeInput(session, "hm_splitvar", choices=xgenes, server=TRUE)
        }
        if(input$hm_splitby=='phenotype') {
            cvar <- sort(pgx.getCategoricalPhenotypes(pgx$samples, min.ncat=2, max.ncat=999))
            sel <- cvar[1]
            cvar0 <- grep("^[.]",cvar,value=TRUE,invert=TRUE) ## no estimated vars
            sel <- head(c(grep("type|family|class|stat",cvar0,ignore.case=TRUE,value=TRUE),
                          cvar0,cvar),1)
            shiny::updateSelectInput(session, "hm_splitvar", choices=cvar, selected=sel)
        }
    })
    
    input_hm_samplefilter <- shiny::reactive({
        input$hm_samplefilter
    }) %>% shiny::debounce(3000)
    
    ## update choices upon change of level
    shiny::observe({

        shiny::req(pgx$families, pgx$gsetX)
        shiny::req(input$hm_level)
        ###if(is.null(input$hm_level)) return(NULL)
        choices = names(pgx$families)        
        if(input$hm_level=="geneset") {
            nk <- sapply(COLLECTIONS, function(k) sum(k %in% rownames(pgx$gsetX)))
            choices = names(COLLECTIONS)[nk>=5]
        }
        choices <- c("<custom>","<contrast>",choices)
        choices <- sort(unique(choices))
        shiny::updateSelectInput(session, "hm_features", choices=choices)
    })
    

    ##================================================================================
    ##========================= REACTIVE FUNCTIONS ===================================
    ##================================================================================
    
    getFilteredMatrix <- shiny::reactive({
        ## Returns filtered matrix ready for clustering. Filtering based
        ## on user selected geneset/features or custom list of genes.
        ##
        ##
        ##
        shiny::req(pgx$X, pgx$Y, pgx$gsetX, pgx$families, pgx$genes)

        genes = as.character(pgx$genes[rownames(pgx$X),"gene_name"])
        genesets = rownames(pgx$gsetX)
        ft <- input$hm_features
        shiny::req(ft)
        
        if(input$hm_level=="geneset") {
            ##-----------------------------------
            ## Gene set level features
            ##-----------------------------------
            gsets = rownames(pgx$gsetX)
            ##gsets = unique(unlist(COLLECTIONS[ft]))
            gsets = unique(COLLECTIONS[[ft]])
            zx = pgx$gsetX
            if(input$hm_customfeatures!="") {
                gsets1 = genesets[grep(input$hm_customfeatures, genesets,ignore.case=TRUE)]
                if(length(gsets1)>2) gsets = gsets1
            }        
            zx = zx[intersect(gsets,rownames(zx)),]
        }

        idx <- NULL
        if(input$hm_level=="gene") {
            
            ##-----------------------------------
            ## Gene level features
            ##-----------------------------------
            gg = pgx$families[[1]]
            if(ft =="<all>") {
                gg = rownames(pgx$X)
            } else if(ft =="<contrast>") {
                ct <- input$hm_contrast
                shiny::req(ct)
                shiny::req(input$hm_ntop)
                fc <- names(sort(pgx.getMetaMatrix(pgx)$fc[,ct]))
                n1 <- floor(as.integer(input$hm_ntop)/2)
                gg <- unique(c(head(fc,n1),tail(fc,n1)))
            } else if(ft %in% names(pgx$families)) {
                gg = pgx$families[[ft]]
            } else if(ft =="<custom>" && ft!="") {
                message("[getFilteredMatrix] selecting for <custom> features")
                customfeatures = "ADORA2A ARHGEF5 BTLA CD160 CD244 CD27 CD274 CD276 CD47 CD80 CEACAM1 CTLA4 GEM HAVCR2 ICOS IDO1 LAG3"
                customfeatures = "CTLA4 GEM HAVCR2 ICOS IDO1 LAG3 PDCD1 TNFSF4 VISTA VTCN1 TIGIT PVR --- CD28 CD40 CD40LG ICOSLG TNFRSF9 TNFSF9 CD70 TNFRSF4 TNFRSF18 --- TNFSF18 SIRPA LGALS9 ARG1 CD86 IDO2 PDCD1LG2 KIR2DL3"
                customfeatures <- input$hm_customfeatures
                gg1 = strsplit(customfeatures,split="[, ;\n\t]")[[1]]

                is.regx <- grepl("[*.?\\[]",gg1[1])
                if(length(gg1)==1 && is.regx) {
                    gg1 <- grep(gg1,genes,ignore.case=TRUE,value=TRUE)
                }
                if(length(gg1)==1 && !is.regx) {
                    gg1 <- c(gg1,gg1) ## heatmap does not like single gene
                }

                gg1 = gg1[toupper(gg1) %in% toupper(genes) | grepl("---",gg1)]
                idx <- NULL
                if(any(grepl("^---",gg1))) {
                    message("[getFilteredMatrix] <custom> groups detected")
                    idx <- rep("F1",length(gg1))
                    names(idx) <- gg1
                    kk <- c(1,grep("^---",gg1),length(gg1)+1)
                    for(i in 1:(length(kk)-1)) {
                        ii <- kk[i]:(kk[i+1]-1)
                        idx[ii] <- paste0("F",i)
                    }
                    gg1 <- gg1[grep("---",gg1,invert=TRUE)]
                    idx <- idx[gg1]
                }
                gg <- gg1
                ##if(length(gg1)>1) gg = gg1
            } else {
                message("[getFilteredMatrix] ERROR!!:: switch error : ft= ",ft)
                gg <- NULL
                return(NULL)
            }
            
            gg <- gg[which(toupper(gg) %in% toupper(genes))]
            jj <- match(toupper(gg), toupper(genes))
            pp  <- rownames(pgx$X)[jj]
            zx = pgx$X[pp,,drop=FALSE]
            if(!is.null(idx)) {
                idx <- idx[gg]
                names(idx) <- rownames(zx)
            }
        }
        if(nrow(zx)==0) return(NULL)

        dim(zx)
        kk <- selectSamplesFromSelectedLevels(pgx$Y, input_hm_samplefilter() )
        zx <- zx[,kk,drop=FALSE]    
        
        if( input$hm_level=="gene" &&
            "chr" %in% names(pgx$genes) &&
            input$hm_filterXY )
        {
            ## Filter out X/Y chromosomes before clustering
            chr.col <- grep("^chr$|^chrom$",colnames(pgx$genes))
            chr <- pgx$genes[rownames(zx),chr.col]
            not.xy <- !(chr %in% c("X","Y",23,24)) & !grepl("^X|^Y|chrX|chrY",chr)
            table(not.xy)
            zx <- zx[which(not.xy), ]
            if(!is.null(idx)) idx <- idx[rownames(zx)]
        }

        if( input$hm_level=="gene" && input$hm_filterMitoRibo )
        {
            ## Filter out X/Y chromosomes before clustering
            is.ribomito <- grepl("^RP[LS]|^MT-",rownames(zx),ignore.case=TRUE)
            table(is.ribomito)
            zx <- zx[which(!is.ribomito),,drop=FALSE]
            if(!is.null(idx)) idx <- idx[rownames(zx)]
        }

        flt <- list(zx=zx, idx=idx)

        return(flt)
    })
    

    getTopMatrix <- shiny::reactive({

        ##pgx <- inputData()
        shiny::req(pgx$X, pgx$samples)

        flt <- getFilteredMatrix()
        zx <- flt$zx
        if(is.null(flt)) return(NULL)
        if(is.null(zx) || nrow(zx)==0) return(NULL)
        
        nmax = 4000
        nmax = as.integer(input$hm_ntop)
        idx <- NULL
        splitvar ="none"
        splitvar <- input$hm_splitvar
        splitby  <- input$hm_splitby
        do.split <- splitby!='none'
        
        if(splitby=="gene" && !splitvar %in% rownames(pgx$X)) return(NULL)
        if(splitby=="phenotype" && !splitvar %in% colnames(pgx$samples)) return(NULL)

        grp <- NULL
        ## split on a phenotype variable
        if(do.split && splitvar %in% colnames(pgx$samples)) {
            dbg("[ClusteringBoard:getTopMatrix] splitting by phenotype: ",splitvar)
            grp <- pgx$samples[colnames(zx),splitvar]
            table(grp)
        }

        ## split on gene expression value: hi vs. low
        if(do.split && splitvar %in% rownames(pgx$X)) {

            ##xgene <- rownames(pgx$X)[1]
            gx <- pgx$X[1,]
            gx <- pgx$X[splitvar,colnames(zx)]
            
            ## estimate best K
            within.ssratio <- sapply(1:4, function(k) {km=kmeans(gx,k);km$tot.withinss/km$totss})
            within.ssratio
            diff(within.ssratio)
            k.est <- min(which(within.ssratio < 0.10))
            k.est <- min(which(abs(diff(within.ssratio)) < 0.10))
            k.est
            k.est <- pmax(pmin(k.est,3),2)
            k.est = 2 ## for now...
            
            if (k.est==2) {
                km <- kmeans(gx, centers=2)
                km.rnk <- rank(km$centers,ties.method="random")
                grp.labels <- c("low","high")[]
                grp <- grp.labels[km$cluster]
            } else if(k.est==3) {
                km  <- kmeans(gx, centers=3)
                km.rnk <- rank(km$centers,ties.method="random")                
                grp.labels <- c("low","mid","high")[km.rnk]
                grp <- grp.labels[km$cluster]                
            } else if(k.est==4) {
                km  <- kmeans(gx, centers=4)
                km.rnk <- rank(km$centers,ties.method="random")
                grp.labels <- c("low","mid-low","mid-high","high")[km.rnk]
                grp <- grp.labels[km$cluster]                
            }
            grp <- paste0(splitvar,":",grp)
            names(grp) <- colnames(zx)
        }        
        ##if(length(grp)==0) splitby <- 'none'
        if(do.split && length(grp)==0) return(NULL)        

        ##------------------------------------------------------------
        ## Any BMC scaling??
        ##------------------------------------------------------------
        if(do.split && input$hm_scale=="BMC") {
            dbg("[ClusteringBoard:getTopMatrix] batch-mean centering...")            
            for(g in unique(grp)) {
                jj <- which(grp == g)
                zx1 <- zx[,jj,drop=FALSE]
                zx[,jj] <- zx1 - rowMeans(zx1,na.rm=TRUE)
            }
        }        
        
        ##------------------------------------------------------------
        ## Create reduced matrix according to topmode
        ##------------------------------------------------------------
        
        topmode="specific"
        topmode="sd"
        topmode <- input$hm_topmode
        if(topmode == "specific" && length(table(grp))<=1) {
            topmode <- "sd"
        }
        
        addsplitgene <- function(gg) {
            if(do.split && splitvar %in% rownames(pgx$X)) {
                gg <- unique(c(splitvar,gg))
            }
            gg
        }
        
        if(!do.split && topmode=="specific") topmode <- "sd"
        grp.zx <- NULL        
        if(topmode=="pca") {
            dbg("[ClusteringBoard:getTopMatrix] splitting by PCA")

            NPCA=5
            svdres <- irlba::irlba(zx - rowMeans(zx), nv=NPCA)
            ntop = 12
            ntop <- as.integer(input$hm_ntop) / NPCA
            gg <- rownames(zx)
            sv.top <- lapply(1:NPCA,function(i) gg[head(order(-abs(svdres$u[,i])),ntop)] )
            gg.top <- unlist(sv.top)
            ##gg.top <- addsplitgene(gg.top) 
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
            ##grp <- pgx$samples[colnames(zx),"cluster"]
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
            dbg("[ClusteringBoard:getTopMatrix] order by SD")
            ii <- order(-apply(zx,1,sd,na.rm=TRUE))
            zx = zx[ii,,drop=FALSE] ## order
            zx = head(zx,nmax)
            ##gg <- addsplitgene(rownames(zx))
            ##zx = zx[gg,,drop=FALSE]
            
        }
        ##zx = zx / apply(zx,1,sd,na.rm=TRUE)  ## scale??
        
        ## ------------- cluster the genes???
        if(!is.null(flt$idx))  {
            idx <- flt$idx[rownames(zx)]  ## override
        }

        CLUSTK = 4  ## number of gene groups (NEED RETHINK)
        CLUSTK <- as.integer(input$hm_clustk)
        if(is.null(idx)) {
            D <- as.dist(1 - cor(t(zx),use="pairwise"))
            system.time( hc <- fastcluster::hclust(D, method="ward.D2" ) )
            ## system.time( hc <- flashClust::hclust(D, method="ward" ) )
            ## system.time( hc <- nclust(D, link="ward") )
            ngrp = min(CLUSTK, nrow(zx))  ## how many default groups???
            idx = paste0("S",cutree(hc, ngrp))                            
        }
        
        ## ------------- matched annotation
        ##annot = pgx$Y[colnames(zx),,drop=FALSE]  ## Y or full matrix??
        annot = pgx$samples[colnames(zx),,drop=FALSE]  ## Y or full matrix??        
        kk = grep("sample|patient",colnames(annot),invert=TRUE)
        annot = annot[,kk,drop=FALSE]  ## no group??    
        samples = colnames(zx) ## original sample list

        ## ----------------------------------------------------
        ## ------------ calculate group summarized ------------
        ## ----------------------------------------------------
        grp.zx <- NULL
        grp.var = "group"
        grp.var <- input$hm_group
        
        if(grp.var %in% colnames(pgx$samples)) { 
            gg.grp <- pgx$samples[colnames(zx),grp.var]
            ## take most frequent term as group annotation value
            grp.zx = tapply( colnames(zx), gg.grp, function(k)
                rowMeans(zx[,k,drop=FALSE],na.rm=TRUE))
            grp.zx = do.call( cbind, grp.zx)
            most.freq <- function(x) names(sort(-table(x)))[1]
            grp.annot = tapply( rownames(annot), gg.grp, function(k) {
                f <- apply(annot[k,,drop=FALSE],2,function(x) most.freq(x))
                w.null <- sapply(f,is.null)
                if(any(w.null))  f[which(w.null)] <- NA
                unlist(f)
            })
            grp.annot = data.frame(do.call( rbind, grp.annot))
            grp.annot = grp.annot[colnames(grp.zx),,drop=FALSE]
        } else {

            grp.zx <- zx
            grp.annot <- annot

        }

        ##input$top_terms
        filt <- list(
            ## mat=zx, annot=annot, 
            mat = grp.zx,
            annot = grp.annot,
            grp = grp,
            idx = idx,
            samples = samples)
        return(filt)
    })


    ##================================================================================
    ##========================= FUNCTIONS ============================================
    ##================================================================================

    hm_splitmap_text = tagsub("Under the <strong>Heatmap</strong> panel, hierarchical clustering can be performed on gene level or gene set level expression in which users have to specify it under the {Level} dropdown list. <p>Under the plot configuration {{Settings}}, users can split the samples by a phenotype class (e.g., tissue, cell type, or gender) using the {split by} setting. In addition, users can specify the top N = (50, 150, 500) features to be used in the heatmap. The ordering of top features is selected under {top mode}. The criteria to select the top features are: <ol><li>SD - features with the highest standard deviation across all the samples, </li><li>specific - features that are overexpressed in each phenotype class compared to the rest, or by </li><li>PCA - by principal components.<br></ol> <br><p>Users can also choose between 'relative' or 'absolute' expression scale. Under the {cexCol} and {cexRow} settings, it is also possible to adjust the cex for the column and row labels.")

    hm1_splitmap.RENDER <- shiny::reactive({    
        
        ##------------------------------------------------------------
        ## ComplexHeatmap based splitted heatmap
        ##------------------------------------------------------------
                
        filt <- getTopMatrix()        
        shiny::req(filt)
        ##if(is.null(filt)) return(NULL)
        
        ##if(input$hm_group) {
        zx <- filt$mat    
        annot = filt$annot                
        zx.idx <- filt$idx

        if(nrow(zx) <= 1) return(NULL)
        
        show_rownames = TRUE
        if(nrow(zx) > 100) show_rownames = FALSE   
        
        cex1 = ifelse(ncol(zx)>50,0.75,1)
        cex1 = ifelse(ncol(zx)>100,0.5,cex1)
        cex1 = ifelse(ncol(zx)>200,0,cex1)
        
        scale.mode = "none"
        if(input$hm_scale=="relative") scale.mode <- "row.center"
        if(input$hm_scale=="BMC") scale.mode <- "row.bmc"        
        scale.mode
        
        ## split genes dimension in 5 groups
        splity = 5
        splity = 6
        if(!is.null(zx.idx)) splity = zx.idx
        
        ## split samples
        splitx = NULL    
        splitx = filt$grp

        show_legend=show_colnames=TRUE
        show_legend <- input$hm_legend
        if(input$hm_level=="geneset" || !is.null(splitx)) show_legend = FALSE
        
        annot$group = NULL  ## no group in annotation??
        show_colnames <- (input$hm_cexCol != 0)
        ##if(ncol(zx) > 200) show_colnames <- FALSE ## never...    
        
        if(input$hm_level=="gene") {
            ## strip any prefix
            rownames(zx) = sub(".*:","",rownames(zx)) 
        }
        rownames(zx) <- sub("HALLMARK:HALLMARK_","HALLMARK:",rownames(zx))
        rownames(zx) = gsub(GSET.PREFIX.REGEX,"",rownames(zx))
        rownames(zx) = substring(rownames(zx),1,50)  ## cut long names...
        if(input$hm_level=="geneset")  rownames(zx) <- tolower(rownames(zx))
        
        cex2 <- ifelse( nrow(zx) > 60, 0.8, 0.9)    
        cex1 <- as.numeric(input$hm_cexCol)*0.85
        cex2 <- as.numeric(input$hm_cexRow)*0.75
        cex0 <- ifelse(!is.null(splitx) && length(splitx)<=10, 1.05, 0.85)  ## title
                       
        crot <- 0
        totnchar <- nchar(paste0(unique(splitx),collapse=""))
        totnchar
        nx <- length(unique(splitx))
        if(!is.null(splitx) & (totnchar > 44 || nx>=6) ) crot=90
        
        nrownames = 60
        nrownames = 9999
        if(input$hm_cexRow==0) nrownames <- 0
        
        if(0) {
            split=splity;splitx=splitx;mar=c(5,25); scale=scale.mode; show_legend=show_legend;
            show_colnames = show_colnames; column_title_rot=crot;
            show_rownames = nrownames; softmax=0;
            ## side.height.fraction=0.03+0.055*NCOL(annot); 
            cexCol=cex1; cexRow=cex2;title_cex=1.0 
            col.annot=annot; row.annot=NULL; annot.ht=2.2;
            nmax=-1
        }

        if(0) {
            dbg("[hm1_splitmap.RENDER] rendering heatmap...")
            dbg("[hm1_splitmap.RENDER] dim(annot) = ", paste(dim(annot),collapse="x"))
            dbg("[hm1_splitmap.RENDER] rownames(annot) = ", rownames(annot))
            dbg("[hm1_splitmap.RENDER] colnames(annot) = ", colnames(annot))
            dbg("[hm1_splitmap.RENDER] splitx = ", paste(splitx,collapse=" "))
            dbg("[hm1_splitmap.RENDER] splity = ", paste(splity,collapse=" "))
        }
        shiny::showNotification('rendering heatmap...')
        plt <- grid::grid.grabExpr(
                         gx.splitmap(
                             zx, 
                             split = splity, splitx = splitx,
                             scale = scale.mode, show_legend = show_legend,
                             show_colnames = show_colnames, column_title_rot = crot,
                             column_names_rot = 45,                             
                             show_rownames = nrownames, rownames_width = 40,
                             softmax = 0,
                             ## side.height.fraction=0.03+0.055*NCOL(annot), 
                             title_cex = cex0, cexCol = cex1, cexRow = cex2, 
                             col.annot = annot, row.annot = NULL, annot.ht = 2.3,
                             key.offset = c(0.89,1.01),
                             main=" ", nmax = -1, mar = c(8,16)
                         )
                     )
        plt
    })

    hm2_splitmap.RENDER <- shiny::reactive({

        ##------------------------------------------------------------
        ## iHeatmap based splitted heatmap
        ##------------------------------------------------------------
        shiny::req(pgx$genes)
        
        ## -------------- variable to split samples        
        ##scale = ifelse(input$hm_scale=="relative","row.center","none")    
        scale = "none"
        if(input$hm_scale=="relative") scale <- "row.center"
        if(input$hm_scale=="BMC") scale <- "row.bmc"        
        scale

        plt <- NULL

        filt <- getTopMatrix()
        ##if(is.null(filt)) return(NULL)
        shiny::req(filt)
        
        ##if(input$hm_group) {
        X <- filt$mat    
        annot = filt$annot        
        idx <- filt$idx
        
        ## sample clustering index
        splitx <- NULL
        splitx <- filt$grp
        
        ## iheatmapr needs factors for sharing between groups
        annotF <- data.frame(as.list(annot),stringsAsFactors=TRUE)
        rownames(annotF) = rownames(annot)
        
        colcex <- as.numeric(input$hm_cexCol)
        rowcex = as.numeric(input$hm_cexRow)
        
        tooltips = NULL
        if(input$hm_level=="gene") {
            getInfo <- function(g) {
                aa = paste0("<b>",pgx$genes[g,"gene_name"],"</b>. ",
                            ## pgx$genes[g,"map"],". ",
                            pgx$genes[g,"gene_title"],".")
                breakstring2(aa, 50, brk="<br>")
            }
            tooltips = sapply(rownames(X), getInfo)
        } else {
            aa = gsub("_"," ",rownames(X)) ## just geneset names
            tooltips = breakstring2(aa, 50, brk="<br>")
        }
        ##genetips = rownames(X)
                
        shiny::showNotification('rendering iHeatmap...')

        plt <- pgx.splitHeatmapFromMatrix(
            X=X, annot=annotF, ytips=tooltips,
            idx=idx, splitx=splitx, scale=scale,
            row_annot_width=0.03, rowcex=rowcex,
            colcex=colcex )
       
        ## DOES NOT WORK...
        ##plt <- plt %>%
        ## plotly::config(toImageButtonOptions = list(format='svg', height=800, width=800))
        
        return(plt)
    })

    topmodes <- c("sd","pca","specific")
    ##if(DEV) topmodes <- c("sd","specific","pca")
    
    hm_splitmap_opts = shiny::tagList(
        withTooltip( shiny::radioButtons(ns("hm_plottype"), "Plot type:",
                             choices=c("ComplexHeatmap","iHeatmap"),
                             selected="ComplexHeatmap", inline=TRUE, width='100%'),
               "Choose plot type: ComplexHeatmap (static) or iHeatmap (interactive)",
               placement="right",options = list(container = "body")),
        withTooltip( shiny::radioButtons(
            ns("hm_splitby"), "Split samples by:", inline=TRUE,
            ## selected="phenotype",
            choices=c("none","phenotype","gene")),
               "Split the samples by phenotype or expression level of a gene.",
               placement="right",options = list(container = "body")),
        shiny::conditionalPanel(
            "input.hm_splitby != 'none'", ns=ns,
            withTooltip( shiny::selectInput(ns("hm_splitvar"), NULL, choices=""),
                   "Specify phenotype or gene for splitting the columns of the heatmap.",
                   placement="right",options = list(container = "body")),
        ),
        shiny::fillRow(
            height = 50,
            withTooltip( shiny::selectInput(ns('hm_topmode'),'Top mode:',topmodes, width='100%'),
                   "Specify the criteria for selecting top features to be shown in the heatmap.",
                   placement = "right", options = list(container = "body")),
            withTooltip( shiny::selectInput(ns('hm_ntop'),'Top N:',c(50,150,500),selected=50),
                   "Select the number of top features in the heatmap.",
                   placement="right", options = list(container = "body")),
            withTooltip( shiny::selectInput(ns('hm_clustk'),'K:',1:6,selected=4),
                   "Select the number of gene clusters.",
                   placement="right", options = list(container = "body"))
        ),
        ##br(),
        withTooltip( shiny::radioButtons(
            ns('hm_scale'), 'Scale:', choices=c('relative','absolute','BMC'), inline=TRUE),
            ## ns('hm_scale'), 'Scale:', choices=c('relative','absolute'), inline=TRUE),
            "Show relative (i.e. mean-centered), absolute expression values or batch-mean-centered.",
            placement="right", options = list(container = "body")),
        withTooltip( shiny::checkboxInput(
            ns('hm_legend'), 'show legend', value=TRUE), "Show or hide the legend.",
            placement="right", options = list(container = "body")),
        shiny::fillRow(
            height = 50,
            ## shiny::checkboxInput(ns("hm_labRow"),NULL),
            withTooltip( shiny::numericInput(ns("hm_cexRow"), "cexRow:", 1, 0, 1.4, 0.1, width='100%'),
                   "Specify the row label size. Set to 0 to suppress row labels.",
                   placement="right",options = list(container = "body")),
            withTooltip( shiny::numericInput(ns("hm_cexCol"), "cexCol:", 1, 0, 1.4, 0.1, width='100%'),
                   "Specify the column label size. Set to 0 to suppress column labels.",
                   placement="right", options = list(container = "body"))            
        ),
        shiny::br()
    )
    
    hm_splitmap_caption = "<b>Clustered heatmap.</b> Heatmap showing gene expression sorted by 2-way hierarchical clustering. Red corresponds to overexpression, blue to underexpression of the gene. At the same time, gene clusters are functionally annotated in the 'Annotate clusters' panel on the right."
    
    output$hm1_splitmap <- shiny::renderPlot({
        plt <- hm1_splitmap.RENDER()
        grid::grid.draw(plt, recording=FALSE)
    }, res=90)
    
    output$hm2_splitmap <- renderIheatmap({
        hm2_splitmap.RENDER()
    })

    hm_splitmap.switchRENDER <- shiny::reactive({    
        ##req(input$hm_plottype)
        p = NULL
        if(input$hm_plottype %in% c("ComplexHeatmap","static") ) {
            p = shiny::plotOutput(ns("hm1_splitmap"), height=fullH-80)  ## height defined here!!
        } else {
            p = iheatmaprOutput(ns("hm2_splitmap"), height=fullH-80) ## height defined here!!
        }
        return(p)
    })
    
    ##output$hm_splitmap_pdf <- shiny::downloadHandler(
    hm_splitmap_downloadPDF <- shiny::downloadHandler(
        filename = "plot.pdf",
        content = function(file) {
            ##PDFFILE = hm_splitmap_module$.tmpfile["pdf"]  ## from above!
            PDFFILE = paste0(gsub("file","plot",tempfile()),".pdf")            
            dbg("[ClusteringBoard] hm_splitmap_downloadPDF: exporting SWITCH to PDF...")
            ##showNotification("exporting to PDF")            
            ##wd <- input$hm_pdfwidth
            ##ht <- input$hm_pdfheight
            ##wd <- input$pdf_width
            ##ht <- input$pdf_height
            wd <- input[["hm_splitmap-pdf_width"]]  ## ugly!!
            ht <- input[["hm_splitmap-pdf_height"]] ## ugly!!
            
            if(1 && input$hm_plottype %in% c("ComplexHeatmap","static")) {
                pdf(PDFFILE, width=wd, height=ht)
                grid::grid.draw(hm1_splitmap.RENDER())
                ##print(hm1_splitmap.RENDER())
                ##hm1_splitmap.RENDER()
                dev.off()
            } else {
                save_iheatmap(hm2_splitmap.RENDER(), filename=PDFFILE,
                              vwidth=wd*100, vheight=ht*100)
            }
            if(WATERMARK) {
                dbg("[ClusteringBoard] adding watermark to PDF...\n")
                addWatermark.PDF(PDFFILE) ## from pgx-modules.R
            }            
            dbg("[ClusteringBoard] hm_splitmap_downloadPDF: exporting done...")
            file.copy(PDFFILE,file)        
        }
    )
    
    hm_splitmap_downloadPNG <- shiny::downloadHandler(
        filename = "plot.png",
        content = function(file) {
            PNGFILE = paste0(gsub("file","plot",tempfile()),".png")            
            dbg("[ClusteringBoard] hm_splitmap_downloadPDF:: exporting SWITCH to PNG...")
            ##showNotification("exporting to PNG")
            wd <- 100*as.integer(input[["hm_splitmap-pdf_width"]])
            ht <- 100*as.integer(input[["hm_splitmap-pdf_height"]])
            if(1 && input$hm_plottype %in% c("ComplexHeatmap","static")) {
                png(PNGFILE, width=wd, height=ht, pointsize=24)
                grid::grid.draw(hm1_splitmap.RENDER())
                ##print(hm1_splitmap.RENDER())  ## should be done inside render for base plot...
                ##hm1_splitmap.RENDER()  ## should be done inside render for base plot...
                ##plot(sin)
                dev.off()
            } else {
                save_iheatmap(hm2_splitmap.RENDER(), filename=PNGFILE,  vwidth=wd, vheight=ht)
            }
            dbg("[ClusteringBoard] hm_splitmap_downloadPNG: exporting done...")            
            file.copy(PNGFILE,file)        
        }
    )

    hm_splitmap_downloadHTML <- shiny::downloadHandler(
        filename = "plot.html",
        content = function(file) {
            ##HTMLFILE = hm_splitmap_module$.tmpfile["html"]  ## from above!
            HTMLFILE = paste0(gsub("file","plot",tempfile()),".html")            
            dbg("renderIheatmap:: exporting SWITCH to HTML...")
            shiny::withProgress({
                ##write("<body>HTML export error</body>", file=HTMLFILE)    
                p <- hm2_splitmap.RENDER()
                shiny::incProgress(0.5)
                save_iheatmap(p, filename=HTMLFILE)
            }, message="exporting to HTML", value=0 )
            dbg("renderIheatmap:: ... exporting done")
            file.copy(HTMLFILE,file)        
        }
    )
    
    
    ## call plotModule
    hm_splitmap_module <- shiny::callModule(
        plotModule,
        id = "hm_splitmap",
        func = hm_splitmap.switchRENDER, ## ns=ns,
        ## func2 = hm_splitmap.switchRENDER, ## ns=ns,
        show.maximize = FALSE,
        plotlib = "generic",
        renderFunc = "renderUI",
        outputFunc = "uiOutput",
        download.fmt = c("pdf","png"),
        options = hm_splitmap_opts,
        height = 2*fullH-80, ##???
        width = '100%',
        pdf.width = 10, pdf.height = 8, 
        title ="Clustered Heatmap",
        info.text = hm_splitmap_text,
        info.width = "350px",
        download.pdf = hm_splitmap_downloadPDF,
        download.png = hm_splitmap_downloadPNG,
        download.html = hm_splitmap_downloadHTML,
        add.watermark = WATERMARK
    )

    ##================================================================================
    ##================================ PCA/tSNE ======================================
    ##================================================================================
    
    hm_PCAplot_text = tagsub(paste0(' The <b>PCA/tSNE</b> panel visualizes unsupervised clustering obtained by the principal components analysis (',a_PCA,') or t-distributed stochastic embedding (',a_tSNE,') algorithms. This plot shows the relationship (or similarity) between the samples for visual analytics, where similarity is visualized as proximity of the points. Samples that are ‘similar’ will be placed close to each other. 
<br><br>Users can customise the PCA/tSNE plot in the plot settings, including the {color} and {shape} of points using a phenotype class, choose t-SNE or PCA layout, label the points, or display 2D and 3D visualisation of the PCA/tSNE plot.'))
    
    shiny::observe({
        shiny::req(pgx$Y)
        ##input$menuitem  ## upon menuitem change
        var.types = colnames(pgx$Y)
        var.types = var.types[grep("sample|patient",var.types,invert=TRUE)]
        vv = c(var.types,rep("<none>",10))
        var.types0 = c("<none>","<cluster>",var.types)
        var.types0 = c("<none>",var.types)
        var.types1 = c("<none>",var.types)
        grp = vv[1]
        if("group" %in% var.types) grp = "group"
        shiny::updateSelectInput(session, "hmpca.colvar", choices=var.types0, selected=grp)
        shiny::updateSelectInput(session, "hmpca.shapevar", choices=var.types1, selected="<none>")
        ##updateSelectInput(session, "hmpca.line", choices=var.types1, selected="<none>")
        ##updateSelectInput(session, "hmpca.text", choices=var.types0, selected="group")
    })

    hm_getClusterPositions <- shiny::reactive({

        dbg("[hm_getClusterPositions] reacted")

        ##pgx <- inputData()
        shiny::req(pgx$tsne2d,pgx$tsne3d,pgx$cluster)

        ## take full matrix
        flt <- getFilteredMatrix()
        zx <- flt$zx
        
        clustmethod="tsne";pdim=2
        do3d <- ("3D" %in% input$hmpca_options)
        pdim = c(2,3)[ 1 + 1*do3d]
        
        pos = NULL
        force.compute = FALSE
        clustmethod = input$hm_clustmethod
        clustmethod0 <- paste0(clustmethod,pdim,"d")
        
        if(clustmethod=="default" && !force.compute) {
            if(pdim==2 && !is.null(pgx$tsne2d) ) {
                pos <- pgx$tsne2d[colnames(zx),]
            } else if(pdim==3 && !is.null(pgx$tsne3d) ) {
                pos <- pgx$tsne3d[colnames(zx),]
            }
        } else if( clustmethod0 %in% names(pgx$cluster$pos))  {
            shiny::showNotification(paste("switching to ",clustmethod0," layout...\n"))
            pos <- pgx$cluster$pos[[clustmethod0]]
            if(pdim==2) pos <- pos[colnames(zx),1:2]
            if(pdim==3) pos <- pos[colnames(zx),1:3]
        } else  {
            ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ## This should not be necessary anymore as we prefer to
            ## precompute all clusterings.
            shiny::showNotification(paste("computing ",clustmethod,"...\n"))            

            ntop = 1000
            ## ntop = as.integer(input$hm_ntop2)    
            zx = zx[order(-apply(zx,1,sd)),,drop=FALSE]  ## OK?
            if(nrow(zx) > ntop) {
                ##zx = head(zx,ntop)  ## OK?
                zx = zx[1:ntop,,drop=FALSE]  ## OK?
            }            
            if("normalize" %in% input$hmpca_options) {
                zx <- scale(t(scale(t(zx))))
            }
            perplexity = max(1,min((ncol(zx)-1)/3, 30))	
            perplexity
            res <- pgx.clusterMatrix(
                zx, dims = pdim, perplexity = perplexity, 
                ntop = 999999, prefix = "C",         
                find.clusters = FALSE, kclust = 1,                 
                row.center = TRUE, row.scale = FALSE,
                method = clustmethod)
            if(pdim==2) pos <- res$pos2d
            if(pdim==3) pos <- res$pos3d
        }

        pos <- pos[colnames(zx),]
        pos = scale(pos) ## scale
        ##colnames(pos) = paste0("dim",1:ncol(pos))
        ##rownames(pos) = colnames(zx)
        
        idx <- NULL
        dbg("[hm_getClusterPositions] done")

        clust = list(pos=pos, clust=idx) 
        return(clust)
    })
    
    hm_PCAplot.RENDER <- shiny::reactive({

        ##pgx <- inputData()
        shiny::req(pgx$Y)

        do3d = ("3D" %in% input$hmpca_options)        
        clust <- hm_getClusterPositions()
        pos <- clust$pos
        sel <- rownames(pos)
        df <- cbind(pos, pgx$Y[sel,])
        if(!is.null(clust$clust)) df[["<cluster>"]] <- clust$clust
        
        colvar = shapevar = linevar = textvar = NULL
        if(input$hmpca.colvar %in% colnames(df)) colvar <- factor(df[,input$hmpca.colvar])
        if(input$hmpca.shapevar %in% colnames(df)) shapevar <- factor(df[,input$hmpca.shapevar])
        ##if(input$hmpca.line %in% colnames(df)) linevar = factor(df[,input$hmpca.line])
        ##if(input$hmpca.text %in% colnames(df)) textvar = factor(df[,input$hmpca.text])
        mode = "markers"
        ann.text = rep(" ",nrow(df))
        if(!do3d && "sample label" %in% input$hmpca_options) ann.text = rownames(df)
        if(!is.null(colvar)) {
            colvar = factor(colvar)
            textvar <- factor(df[,input$hmpca.colvar])
        }
        symbols = c('circle','square','star','triangle-up','triangle-down','pentagon',
                    'bowtie','hexagon', 'asterisk','hash','cross','triangle-left',
                    'triangle-right','+',c(15:0))


        Y <- cbind("sample"=rownames(pos), pgx$Y[sel,])
        ##tt.info <- paste('Sample:', rownames(df),'</br>Group:', df$group)
        tt.info <- apply(Y, 1, function(y) paste0(colnames(Y),": ",y,"</br>",collapse=""))
        tt.info <- as.character(tt.info)
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
            plt <- plotly::plot_ly(df, mode=mode) %>%
                plotly::add_markers(x = df[j0,1], y = df[j0,2], z = df[j0,3], type="scatter3d",
                            color = colvar[j0], ## size = sizevar, sizes=c(80,140),
                            ##marker = list(size = 5*cex1),
                            marker = list(size=5*cex1, line=list(color="grey10", width=0.1)),
                            symbol = shapevar[j0], symbols=symbols,
                            text = tt.info[j0] ) %>%
                plotly::add_annotations(x = pos[,1], y = pos[,2], z = pos[,3],
                                text = ann.text,
                                ##xref = "x", yref = "y",
                                showarrow = FALSE)
            if(!is.null(j1) & length(j1)>0) {
                plt <- plt %>%  plotly::add_markers(
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
                plt <- plt %>% plotly::add_annotations(
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
            plt <- plotly::plot_ly(df, mode=mode) %>%
                plotly::add_markers(x = df[j0,1], y = df[j0,2], type="scatter",
                            color = colvar[j0], ## size = sizevar, sizes=c(80,140),
                            marker = list(size=16*cex1, line=list(color="grey20", width=0.6)),
                            symbol = shapevar[j0], symbols=symbols,
                            text = tt.info[j0] ) %>%
                plotly::add_annotations(x = pos[,1], y = pos[,2],
                                text = ann.text,
                                ##xref = "x", yref = "y",
                                showarrow = FALSE)

            ## add node labels
            if(!is.null(j1) & length(j1)>0 ) {
                plt <- plt %>%  plotly::add_markers(
                                    x = df[j1,1], y = df[j1,2], type="scatter",
                                    color = colvar[j1], ## size = sizevar, sizes=c(80,140),
                                    marker = list(size=16*cex1, line=list(color="grey20", width=1.8)),
                                    symbol = shapevar[j1], symbols=symbols,
                                    text=tt.info[j1])
            }

            ## add group/cluster annotation labels
            req(input$hmpca_legend)
            if(input$hmpca_legend == 'inside') {
                plt <- plt %>%
                    plotly::layout(legend = list(x=0.05, y=0.95))
            } else if(input$hmpca_legend == 'bottom') {
                plt <- plt %>%
                    plotly::layout(legend = list(orientation='h'))
            } else {
                if(!is.null(textvar) && length(unique(textvar))>1) {
                    grp.pos <- apply(pos,2,function(x) tapply(x,as.character(textvar),median))
                    cex2 <- 1
                    if(length(grp.pos)>20) cex2 <- 0.8
                    if(length(grp.pos)>50) cex2 <- 0.6
                    plt <- plt %>% plotly::add_annotations(
                                               x = grp.pos[,1], y = grp.pos[,2],
                                               text = paste0("<b>",rownames(grp.pos),"</b>"),
                                               font = list(size=24*cex2, color='#555'),
                                               showarrow = FALSE)
                }
                plt <- plt %>%
                    plotly::layout(showlegend = FALSE)
            }


        }
        title = paste0("<b>PCA</b>  (",nrow(pos)," samples)")
        if(input$hm_clustmethod=="tsne") title = paste0("<b>tSNE</b>  (",nrow(pos)," samples)")
        ## plt <- plt %>% plotly::layout(title=title) %>% 
        ##     plotly::config(displayModeBar = FALSE)
        plt <- plt %>% 
            ##config(displayModeBar = FALSE) %>%
            plotly::config(displayModeBar = TRUE) %>%
            ##config(modeBarButtonsToRemove = all.plotly.buttons ) %>%
            plotly::config(displaylogo = FALSE) %>% 
            plotly::config(toImageButtonOptions = list(format='svg', height=800, width=800))
        ##print(plt)
        return(plt)
    })

    hm_PCAplot_opts = shiny::tagList(
        tipifyR( shiny::selectInput( ns("hmpca.colvar"), "Color/label:", choices=NULL, width='100%'),
               "Set colors/labels according to a given phenotype."),
        tipifyR( shiny::selectInput( ns("hmpca.shapevar"), "Shape:", choices=NULL, width='100%'),
               "Set shapes according to a given phenotype."),
        tipifyR( shiny::radioButtons(
                            ns('hmpca_legend'), label =  "Legend:",
                            choices = c('group label','bottom'), inline=TRUE),
                "Normalize matrix before calculating distances."),
        tipifyR( shiny::checkboxGroupInput( ns('hmpca_options'),"Other:",
                                   choices=c('sample label','3D','normalize'), inline=TRUE),
                "Normalize matrix before calculating distances."),
        tipifyR( shiny::radioButtons( ns('hm_clustmethod'),"Layout:",
                              c("default","tsne","pca","umap"),inline=TRUE),
                "Choose the layout method for clustering to visualise.")
    )
    
    hm_PCAplot_caption <- shiny::reactive({
        text1 = "The plot visualizes the similarity of samples as a scatterplot in reduced dimension (2D or 3D). Samples that are similar (in expression) are clustered near to each other, while samples with different expression are positioned farther away. Groups of samples with similar profiles will appear as <i>clusters</i> in the plot."
        if(input$hmpca.colvar!="<none>") {
            text1 <- paste(text1, "Colors correspond to the <strong>",input$hmpca.colvar,"</strong>phenotype.")
        }
        if(input$hmpca.shapevar!="<none>") {
            text1 <- paste(text1, "Shapes correspond to the <strong>",input$hmpca.shapevar,"</strong>phenotype.")
        }
        return(shiny::HTML(text1))
    })

    pca_caption_static = "<b>PCA/tSNE plot.</b> The plot visualizes the similarity in expression of samples as a scatterplot in reduced dimension (2D or 3D). Samples that are similar are clustered near to each other, while samples with different expression are positioned farther away. Groups of samples with similar profiles will appear as <i>clusters</i> in the plot."

    shiny::callModule(
        plotModule, 
        id = "hm_PCAplot",
        func = hm_PCAplot.RENDER, ## ns=ns,
        plotlib = "plotly", 
        options = hm_PCAplot_opts,
        height = c(fullH-80,700), width=c("auto",800),
        pdf.width=8, pdf.height=8,
        title="PCA/tSNE plot",
        info.text = hm_PCAplot_text,
        add.watermark = WATERMARK        
    )
    
    ##================================================================================
    ## Parallel coordinates
    ##================================================================================

    hm_parcoord.ranges <- shiny::reactiveValues()
    
    hm_parcoord.matrix <- shiny::reactive({

        filt <- getTopMatrix()
        shiny::req(filt)
        zx <- filt$mat[,]
        if(input$hm_pcscale) {
            zx <- t(scale(t(zx)))
        }
        rr <- shiny::isolate(shiny::reactiveValuesToList(hm_parcoord.ranges))
        nrange <- length(rr)
        for(i in names(rr)) hm_parcoord.ranges[[i]] <- NULL
        zx <- round(zx, digits=3)
        list(mat=zx, clust=filt$idx)
    })
    
    hm_parcoord.RENDER <- shiny::reactive({
        
        pc <- hm_parcoord.matrix()
        shiny::req(pc)
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
        klrpal = rep(RColorBrewer::brewer.pal(8,"Set2"),99)
        ##klrpal = rep(c("red","blue","green","yellow","magenta","cyan","black","grey"),99)
        klrpal = klrpal[1:max(clust.id)]
        ##klrpal <- setNames(klrpal, sort(unique(clust.id)))
        klrpal2 <- lapply(1:length(klrpal),function(i) c((i-1)/(length(klrpal)-1),klrpal[i]))
        
        plt <-  plotly::plot_ly(df, source = "pcoords") %>%    
            plotly::add_trace(type = 'parcoords', 
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
            plotly::layout(margin = list(l=60, r=60, t=0, b=30)) %>%
            ##config(displayModeBar = FALSE) %>%
            ##config(modeBarButtonsToRemove = setdiff(all.plotly.buttons,"toImage") ) %>%
            plotly::config(toImageButtonOptions = list(format='svg', width=900, height=350, scale=1.2)) %>%
            plotly::config(displaylogo = FALSE) %>%
            plotly::event_register("plotly_restyle")

        plt 

    })

    shiny::observeEvent( plotly::event_data("plotly_restyle", source = "pcoords"), {
        ## From: https://rdrr.io/cran/plotly/src/inst/examples/shiny/event_data_parcoords/app.R
        ##
        d <- plotly::event_data("plotly_restyle", source = "pcoords")
        ## what is the relevant dimension (i.e. variable)?
        dimension <- as.numeric(stringr::str_extract(names(d[[1]]), "[0-9]+"))
        ## If the restyle isn't related to a dimension, exit early.
        if (!length(dimension)) return()
        if (is.na(dimension)) return()
        
        pc <- hm_parcoord.matrix()
        shiny::req(pc)
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

    hm_parcoord.selected <- shiny::reactive({

        mat <- hm_parcoord.matrix()$mat
        clust <- hm_parcoord.matrix()$clust
        shiny::req(mat)
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


    hm_parcoord_opts = shiny::tagList(
        withTooltip( shiny::checkboxInput(ns('hm_pcscale'),'Scale values',TRUE),
               "Scale expression values to mean=0 and SD=1.",
               placement="right",options = list(container = "body"))
    )


    hm_parcoord_text = tagsub("The <strong>Parallel Coordinates</strong> panel
displays the expression levels of selected genes across all conditions in the analysis. On the x-axis the experimental conditions are plotted. The y-axis shows the expression level of the genes grouped by condition. The colors correspond to the gene groups as defined by the hierarchical clustered heatmap.")

    shiny::callModule(
        plotModule,     
        ## hm_parcoord_module <- plotModule(
        "hm_parcoord",
        func = hm_parcoord.RENDER, ## ns = ns,
        plotlib = "plotly", ## renderFunc="renderPlotly",
        ## download.fmt = c("png","pdf","html"),  ## PNG & PDF do not work!!! 
        ## download.fmt = c("html"),
        options = hm_parcoord_opts,
        height = c(0.45*fullH,600), width = c("100%",1000),
        pdf.width=10, pdf.height=6, info.width="350px",
        title = "Parallel coordinates", label = "a",
        info.text = hm_parcoord_text,
        add.watermark = WATERMARK        
        ## caption = hm_parcoord_text,
    )

    hm_parcoord_table.RENDER <- shiny::reactive({

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
            DT::formatSignif(numeric.cols,3) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') 
    })

    hm_parcoord_table_info = "In this table, users can check mean expression values of features across the conditions for the selected genes."


    hm_parcoord_table_module <- shiny::callModule(
        tableModule, id = "hm_parcoord_table",
        func = hm_parcoord_table.RENDER, ## ns=ns,
        info.text = hm_parcoord_table_info,
        title = "Selected genes", label="b",
        height = c(270,700)
    )

    ##================================================================================
    ## Annotate clusters
    ##================================================================================

    clustannot_plots_text = paste0('The top features of the heatmap in the <code>Heatmap</code> panel are divided into gene (or gene set) clusters based on their expression profile patterns. For each cluster, the platform provides a functional annotation in the <code>Annotate cluster</code> panel by correlating annotation features from more than 42 published reference databases, including well-known databases such as ',a_MSigDB,', ',a_KEGG,' and ',a_GO,'. In the plot settings, users can specify the level and reference set to be used under the <code>Reference level</code> and <code>Reference set</code> settings, respectively.')
    
    shiny::observe({

        ##pgx <- inputData()    
        shiny::req(pgx$X,pgx$gsetX,pgx$families)

        if(is.null(input$xann_level)) return(NULL)
        ann.types=sel=NULL
        if(input$xann_level!="phenotype") {
            if(input$xann_level=="geneset") {
                ann.types <- names(COLLECTIONS)
                cc = sapply(COLLECTIONS,function(s) length(intersect(s,rownames(pgx$gsetX))))
                ann.types <- ann.types[cc>=3]
            }
            if(input$xann_level=="gene") {
                ann.types <- names(pgx$families)
                cc = sapply(pgx$families,function(g) length(intersect(g,rownames(pgx$X))))
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
        shiny::updateSelectInput(session, "xann_refset", choices=ann.types, selected=sel)    
    })


    getClustAnnotCorrelation <- shiny::reactive({
        
        ##pgx <- inputData()
        shiny::req(pgx$X,pgx$Y,pgx$gsetX,pgx$families)

        filt <- getTopMatrix()
        shiny::req(filt)
        
        zx  <- filt$mat
        idx <- filt$idx
        samples <- filt$samples

        if(nrow(zx) <= 1) return(NULL)
        
        ann.level="geneset"
        ann.refset="Hallmark collection"
        ann.level = input$xann_level
        ##if(is.null(ann.level)) return(NULL)
        ann.refset = input$xann_refset
        ##if(is.null(ann.refset)) return(NULL)
        shiny::req(input$xann_level, input$xann_refset)
        
        ref = NULL
        ref = pgx$gsetX[,,drop=FALSE]    
        ref = pgx$X[,,drop=FALSE]    
        if(ann.level=="gene" && ann.refset %in% names(pgx$families) ) {
            gg = pgx$families[[ann.refset]]
            jj = match(toupper(gg), toupper(pgx$genes$gene_name))
            jj <- setdiff(jj,NA)
            pp = rownames(pgx$genes)[jj]
            ref = pgx$X[intersect(pp,rownames(pgx$X)),,drop=FALSE]    
        }
        if(ann.level=="geneset" && ann.refset %in% names(COLLECTIONS)) {
            ss = COLLECTIONS[[ann.refset]]
            ss = intersect(ss, rownames(pgx$gsetX))
            length(ss)
            ref = pgx$gsetX[ss,]    
        }
        if(ann.level=="phenotype") {
            ref = t(expandAnnotationMatrix(pgx$Y))
        }
        if(is.null(ref)) {
            cat("<clustering:getClustAnnotCorrelation> WARNING:: ref error\n")
            return(NULL)
        }
        
        ##-----------  restrict to top??
        dim(ref)
        if(nrow(ref)>1000) {
            ref = head(ref[order(-apply(ref,1,sd)),],1000)
        }

        ##-----------  get original data level
        X = pgx$X
        if(input$hm_level=="geneset") X <- pgx$gsetX
        
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
            ##gmt <- GSETS[rownames(rho)]
            gmt <- getGSETS(rownames(rho))
            bg.genes <- toupper(rownames(X))
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

    clustannot_plots.PLOTLY <- shiny::reactive({

        rho = getClustAnnotCorrelation()
        ##if(is.null(rho)) return(NULL)
        shiny::req(rho)
        
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

        klrpal = rep(RColorBrewer::brewer.pal(8,"Set2"),2)
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
            
            plot_list[[i]] <- plotly::plot_ly(
                x=x, y=y, type='bar',  orientation='h',
                ## text=y,
                hoverinfo = 'text',
                hovertemplate = paste0("%{y}<extra>",colnames(rho)[i],"</extra>"),
                ##hovertemplate = "%{y}",
                marker = list(color=klrpal[i])) %>%
                plotly::layout(
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
                plotly::add_annotations(xref = 'paper', yref = 'y',
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
        
        plotly::subplot( plot_list, nrows=nrows, shareX=TRUE,
                ## template = "plotly_dark",
                margin = c(0, 0.0, 0.05, 0.05) ) %>%
            plotly::config(displayModeBar = FALSE) 
    })

    clustannot_plots_opts = shiny::tagList(
        withTooltip( shiny::selectInput(ns("xann_level"), "Reference level:",
                            choices=c("gene","geneset","phenotype"),
                            selected="geneset", width='80%'),
               "Select the level of an anotation analysis.",
               placement="left", options = list(container = "body")),
        shiny::conditionalPanel(
            "input.xann_level == 'geneset'", ns=ns,
            withTooltip( shiny::checkboxInput(ns("xann_odds_weighting"), "Fisher test weighting"),
                   "Enable weighting with Fisher test probability for gene sets. This will effectively penalize small clusters and increase robustness.",
                   placement="left", options = list(container = "body"))
        ),
        withTooltip( shiny::selectInput( ns("xann_refset"), "Reference set:", choices="", width='80%'),
               "Specify a reference set to be used in the annotation.",
               placement="left",options = list(container = "body"))
    )
    
    ##clustannot_plots_module <- plotModule(
    shiny::callModule(
        plotModule, 
        id="clustannot_plots", ##ns=ns,
        ##func=clustannot_plots.RENDER, plotlib = "base",
        func = clustannot_plots.PLOTLY, plotlib="plotly",
        download.fmt = c("png","pdf"),
        options = clustannot_plots_opts,
        height = c(360,600), width = c(500,1000),
        pdf.width=8, pdf.height=5, res=80,
        title="Functional annotation of clusters", label="a",
        info.text = clustannot_plots_text,
        add.watermark = WATERMARK                
    )
    
    clustannot_table.RENDER <- shiny::reactive({
        
        rho = getClustAnnotCorrelation()
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
    clustannot_table_module <- shiny::callModule(
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
    
    output$hm_annotateUI <- shiny::renderUI({
        shiny::fillCol(
            flex = c(1.4,1,NA),
            height = fullH,
            plotWidget(ns("clustannot_plots")),
            plotWidget(ns("clustannot_table")),
            shiny::div(shiny::HTML(clustannot_caption), class="caption"),
        )
    })
    
    ##================================================================================
    ## Phenotypes {data-height=800}
    ##================================================================================

    clust_phenoplot.RENDER <- shiny::reactive({
        
        ##pgx <- inputData()
        shiny::req(pgx$Y)
        
        ## get t-SNE positions
        clust <- hm_getClusterPositions()
        ##pos = pgx$tsne2d
        pos = clust$pos

        Y <- pgx$Y[rownames(pos),,drop=FALSE]
        pheno = colnames(Y)    

        ## don't show these...
        pheno <- grep("batch|sample|donor|repl|surv",pheno,
                      invert=TRUE, ignore.case=TRUE,value=TRUE)
        
        ## layout
        par(mfrow = c(3,2), mar=c(0.3,0.7,2.8,0.7))
        if(length(pheno)>=6) par(mfrow = c(4,3), mar=c(0.3,0.4,2.8,0.4)*0.8)
        if(length(pheno)>=12) par(mfrow = c(5,4), mar=c(0.2,0.2,2.5,0.2)*0.8)
        i=1    

        cex1 <- 1.1*c(1.8,1.3,0.8,0.5)[cut(nrow(pos),breaks=c(-1,40,200,1000,1e10))]    
        cex1 = cex1 * ifelse(length(pheno)>6, 0.8, 1)
        cex1 = cex1 * ifelse(length(pheno)>12, 0.8, 1)

        for(i in 1:min(20,length(pheno))) {

            ## ------- set colors
            colvar = factor(Y[,1])
            colvar = factor(Y[,pheno[i]])
            colvar[which(colvar %in% c(NA,""," ","NA","na"))] <- NA
            colvar = factor(as.character(colvar))
            klrpal = COLORS  
            klr1 = klrpal[colvar]
            klr1 = paste0(gplots::col2hex(klr1),"99")
            jj = which(is.na(klr1))
            if(length(jj)) klr1[jj] <- "#AAAAAA22"
            tt = tolower(pheno[i])

            ## ------- start plot
            base::plot( pos[,], pch=19, cex=cex1, col=klr1,
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
                cex2 = 0.9*cex1**0.33
                text( grp.pos, labels=boxes, cex=cex2*0.95, col="#CCCCCC99")
                text( grp.pos, labels=labels, font=2, cex=cex2)
            }
        }    
    })

    clust_phenoplot.opts = shiny::tagList(
        shiny::radioButtons(ns('clust_phenoplot_labelmode'),"Label",c("groups","legend"),inline=TRUE)
    )

    clust_phenoplot_info = tagsub("<strong>Phenotype distribution.</strong> This figure visualizes the distribution of the available phenotype data. You can choose to put the group labels in the figure or as separate legend in the {Label} setting, in the plot {{settings}}")
    
    ## clust_phenoplot.module <- plotModule(
    shiny::callModule(
        plotModule,         
        "clust_phenoplot", ## ns=ns,
        func  = clust_phenoplot.RENDER, ## plotlib="base",
        func2 = clust_phenoplot.RENDER, ## plotlib="base",
        options = clust_phenoplot.opts,
        height = c(fullH-80,700), res = 85,
        pdf.width = 6, pdf.height = 9, 
        info.text = clust_phenoplot_info,
        add.watermark = WATERMARK
    )

    
    ##=============================================================================
    ## Feature ranking
    ##=============================================================================
    
    calcFeatureRanking <- shiny::reactive({

        shiny::req(pgx$X, pgx$Y, pgx$gsetX, pgx$genes)

        features=X=NULL
        if(input$hm_level=="geneset") {
            features = COLLECTIONS
            X = pgx$gsetX
        } else {
            features = pgx$families
            X = pgx$X
        }

        ## ------------ intersect features, set minimum set size
        rownames(X) <- toupper(rownames(X))
        genes <- toupper(rownames(X))
        features <- lapply(features, toupper)
        features <- lapply(features, function(f) intersect(toupper(f), genes))
        features <- features[sapply(features,length) >=10 ]

        dbg("[calcFeatureRanking] length(features)=",length(features))
        
        ## ------------ Just to get current samples
        ##samples = colnames(X)
        samples <- selectSamplesFromSelectedLevels(pgx$Y, input_hm_samplefilter() )
        X = X[,samples]
        cvar <- pgx.getCategoricalPhenotypes(pgx$Y, max.ncat=999)
        cvar <- grep("sample|patient|years|days|months|gender",
                     cvar,invert=TRUE,value=TRUE) ## no sample IDs
        cvar
        Y = pgx$Y[colnames(X),cvar,drop=FALSE]
        kk = which(apply(Y,2,function(y) length(unique(y))>1))
        Y = Y[,kk,drop=FALSE]
        dim(Y)

        dbg("[calcFeatureRanking] dim(X)=",dim(X))
        dbg("[calcFeatureRanking] dim(Y)=",dim(Y))
        
        ## ------------ Note: this takes a while. Maybe better precompute off-line...
        sdx = apply(X,1,sd)
        names(sdx) = rownames(X)
        S = matrix(NA, nrow=length(features), ncol=ncol(Y))
        rownames(S) = names(features)
        colnames(S) = colnames(Y)

        ## ------------ Create a Progress object
        if(!interactive())  {
            progress <- shiny::Progress$new()
            on.exit(progress$close())    
            progress$set(message = "Calculating feature-set scores", value = 0)
        }
        
        gene.level = TRUE
        gene.level = (input$hm_level=="gene")
        i=1
        for(i in 1:ncol(Y)) {

            if(!interactive()) progress$inc(1/ncol(Y))

            grp = Y[,i]
            grp = as.character(grp)
            
            cat("[calcFeatureRanking] head(grp)=",head(grp),"\n")
            
            score = rep(NA, length(features))
            names(score) = names(features)
            j=1
            for(j in 1:length(features)) {

                pp = features[[j]]                
                if(gene.level) {
                    pp = filterProbes(pgx$genes, features[[j]])
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
                    suppressWarnings( fit <- limma::eBayes( limma::lmFit( X1[,jj], design)) )
                    suppressWarnings( suppressMessages( top <- limma::topTable(fit) ))
                    ##s2 = mean(-log10(top$P.Value))  ## as score
                    s2 = mean(-log10(1e-99 + top$adj.P.Val),na.rm=TRUE)  ## as score
                }
                
                f = 1
                f <- (1 - exp(-(length(pp)/20)**2)) ## penalize smaller sets
                score[j] = f * (s1 * s2) ** ifelse(method=="meta",0.5,1)
                
            }
            S[,i] = score
        }
        S[is.na(S)] <- 0 ## missing values
        return(S)
    })
    
    clust_featureRank.RENDER <- shiny::reactive({

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


    clust_featureRank.opts =  shiny::tagList(
        withTooltip( shiny::radioButtons( ns('clust_featureRank_method'),'Method:',
                             choices=c("p-value","correlation","meta"),
                             inline=TRUE),
               "Choose ranking method: p-value based or correlation-based.",
               placement="right", options = list(container = "body") )
    )

    shiny::callModule(
        plotModule, 
        id="clust_featureRank",
        title="Feature-set ranking", 
        func  = clust_featureRank.RENDER,
        func2 = clust_featureRank.RENDER,
        options = clust_featureRank.opts,
        pdf.width=8, pdf.height=10,
        height = c(fullH-80,700),
        width=c("auto",800),
        res = c(72,90),
        info.text = clust_featureRank_info,
        add.watermark = WATERMARK                
    )

  })  ## end of moduleServer
} ## end of Board
