##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##----------------------------------------------------------------------
## mixOmics related functions
##----------------------------------------------------------------------


mixHivePlot <- function(res, ngs, ct, showloops=FALSE, numlab=6, cex=1)
{

    cat("<mixHivePlot> called\n")
    if(is.null(showloops)) showloops <- FALSE
    
    ##install.packages("HiveR")
    
    gr <- res$graph
    ##edges <- data.frame(igraph::get.edgelist(gr), rho=abs(igraph::E(gr)$weight))
    ##colnames(edges)[1:2] <- c("from","to")
    ##edges$rho <- abs(edges$rho)

    ##-------------------------------------------------------------
    ## Prepare the HivePlot data
    ##-------------------------------------------------------------
    df <- data.frame(res$edges)
    Matrix::head(df)
    ##cat("<mixHivePlot> 1: Matrix::head(df)=",head(df),"\n")
    ##sum(df[,1]==df[,2])    
    hpd <- edge2HPD(df, axis.cols=rep("grey",3))        

    hpd$edges$looping <- res$edges$looping   
    hpd$edges$weight  <- res$edges$importance
    loop.nodes <- unique(
        c(as.character(res$edges$from)[res$edges$looping],
          as.character(res$edges$to)[res$edges$looping]))
    hpd$nodes$looping <- hpd$nodes$lab %in% loop.nodes
    
    hpd <- mineHPD(hpd, option = "rad <- tot.edge.count")
    hpd$nodes$degree <- hpd$nodes$radius
    ##hpd <- mineHPD(hpd, option = "axis <- source.man.sink")
    ##hpd <- mineHPD(hpd, option = "remove zero edge")    
    ##axis <- sub(":.*","",igraph::V(res$graph)$name)        
    hpd$edges$from <- hpd$nodes$lab[hpd$edges$id1]
    hpd$edges$to   <- hpd$nodes$lab[hpd$edges$id2]
    
    if(is.null(ct) || is.null(ngs)) {
        fx <- tapply( igraph::V(gr)$importance, sub(".*:","",igraph::V(gr)$name), max)        
    } else if(ct %in% names(ngs$gx.meta$meta)) {
        ## use fold change as radial layout
        fc <- ngs$gx.meta$meta[[ct]]$meta.fx
        names(fc) <- rownames(ngs$gx.meta$meta[[ct]])        
        gs <- ngs$gset.meta$meta[[ct]]$meta.fx
        names(gs) <- rownames(ngs$gset.meta$meta[[ct]])        
        fc <- fc / max(abs(fc),na.rm=TRUE)
        gs <- gs / max(abs(gs),na.rm=TRUE)
        fx <- c(fc, gs)
    } else if(ct %in% colnames(ngs$samples)) {
        group <- ngs$samples[,ct]
        ## use groupwise-SD as radial layout
        ##matx <- do.call(cbind, tapply(1:ncol(X), group,
        ##           function(i) rowMeans(X[,i,drop=FALSE])))
        ## use F-statistics as radial layout       
        design <- model.matrix(~ 0+group)
        colnames(design) <- sub("group","",colnames(design))
        fit <- limma::eBayes(limma::lmFit(ngs$X,design))
        stat <- limma::topTable(fit, number=Inf)
        fx <- stat$F
        names(fx) <- rownames(ngs$X)
    } else {
        stop("FATAL:: mixHivePlot: unknown contrast/conditions=",ct,"\n")
        ##fx <- tapply( igraph::V(gr)$importance, sub(".*:","",igraph::V(gr)$name), max)        
    }
    Matrix::head(fx)
    fx <- fx / max(abs(fx),na.rm=TRUE)
    g <- sub("[1-9]:","",hpd$nodes$lab)
    hpd$nodes$radius <- rank(fx[g],na.last="keep") 
    hpd$nodes$radius <- 100 * hpd$nodes$radius / max(hpd$nodes$radius,na.rm=TRUE)
    
    maxgrp <- unlist(lapply(res$W,function(w) max.col(w)))
    ##maxgrp <- colnames(res$W[[1]])[maxgrp]
    names(maxgrp) <- as.vector(sapply(res$W,rownames))
    ##maxgrp <- maxgrp[igraph::V(gr)$name]
    
    ## use importance as node size
    importance <- igraph::V(gr)$importance
    names(importance) <- igraph::V(gr)$name
    ##hpd$nodes$size  <- hpd$nodes$degree
    hpd$nodes$size  <- abs(importance[hpd$nodes$lab])
    hpd$nodes$size  <- 1.6*(hpd$nodes$size/max(hpd$nodes$size))**0.5 
    hpd$nodes$axis  <- as.integer(sub(":.*","",hpd$nodes$lab))
    hpd$nodes$color <- c("red3","blue2")[1 + 1*(fx[g]>0)]
    ## hpd$nodes$color <- RColorBrewer::brewer.pal(8,"Set2")[maxgrp[hpd$nodes$lab]]
    
    wt1 <- hpd$edges$weight ## edge.importance
    wt1 <- rank(abs(wt1),na.last="keep") * sign(wt1)
    hpd$edges$weight <- 3 * abs(wt1/max(abs(wt1)))**2
    hpd$edges$color  <- psych::alpha("grey70",0.3)

    ##hpd$edges$color  <- c("grey80","grey40")[1 + 1*(hpd$edges$rho>0)]
    jj = which(hpd$edges$looping)
    if(showloops && length(jj)) {
        ##jj = which(paste(hpd$edges$id1,hpd$edges$id2,sep=">") %in% top.edges)
        hpd$edges[jj,]
        hpd$edges$color  <- psych::alpha("grey70",0.2)
        hpd$edges$color[jj] <- psych::alpha("red3",0.3)
    }
    
    axis.names <- names(res$X)
    makeAcronym <- function(x) {
        x <- gsub("[)(]","",x)
        sapply(strsplit(x,split="[_ -]"), function(s) {
            if(length(s)==1) return(substring(s,1,2))
            toupper(paste(substring(s,1,1),collapse=""))
        })
    }
    axis.names <- sapply(names(res$X), makeAcronym)  ## see pgx-functions
    axis.names

    ##-------------------------------------------------------------
    ## Finally do the plotting
    ##-------------------------------------------------------------
    hpd$nodes$size = cex * hpd$nodes$size
    hpd$edges$weight = cex * hpd$edges$weight

    
    
    mr = max(hpd$nodes$radius)    
    plotHive(hpd, ch=5, bkgnd="white",
             axLabs = axis.names,
             axLab.pos = c(1,1.2,1.2)*0.15*mr, ## rot = c(0,30,-30),
             ## arrow = c("radius units", 0, 20, 60, 25, 40),
             axLab.gpar = grid::gpar(col = "black", fontsize = 18*cex,
                               lwd = 4, fontface="bold")
             )

    ##title="HivePlot"
    tt = paste( "edge.width = edge.importance",
               "node.size = variable importance",
               "axis = fold-change or F-stat", sep="\n")
    grid::grid.text(tt,  x = 0.2*mr, y = -1.0*mr, default.units = "native",
              just="left", gp = grid::gpar(fontsize = 9, col="black"))   
    
    ## axis 1
    rot.xy <- function(x,y,deg) {
        a = deg*pi/180
        rx = cos(a)*x - sin(a)*y
        ry = sin(a)*x + cos(a)*y
        cbind(rx,ry)
    }        
    rot <- c(0, 120, 240)
    mr = max(hpd$nodes$radius)
    yoff <- c(0, -0, +0)

    k=1
    for(k in 1:3) {

        kk <- which(hpd$nodes$axis==k)
        if(showloops) {
            kk <- which(hpd$nodes$axis==k & hpd$nodes$looping)
        }
        jj <- Matrix::head(kk[order(-hpd$nodes$size[kk])],numlab) ## number of labels
        rr <- hpd$nodes$radius
        rx <- rot.xy(0, rr[jj] + 5, rot[k])

        
        lab = sub(".*:","",hpd$nodes$lab[jj])
        pt <- maptools::pointLabel(rx[,1],rx[,2],labels=lab,cex=cex*2,doPlot=FALSE)
        px <- cbind(pt$x, pt$y)
        px[,1] <- px[,1] + 4
        grid::grid.text( lab,
                  ##x = 10 + rx[,1], y = rx[,2],
                  x = px[,1], y = px[,2],
                  ##rot=rot[k] + ifelse(k==1,0,180),
                  default.units = "native", just="left",
                  ##check.overlap=TRUE,
                  gp = grid::gpar(fontsize = 12*cex, col="black"))   
                      
        grid::grid.segments(rx[,1], rx[,2], px[,1] - 1, px[,2],
                      default.units = "native")
        
    }

}

mixPlotLoadings <- function(res, showloops=FALSE, cex=1)
{

    cat("<mixPlotLoadings> called\n")
    levels <- levels(res$Y)
    levels
    ny <- length(levels)
    ny

    
    klrpal <- c("blue2","orange2")
    klrpal <- rep(RColorBrewer::brewer.pal(n=8, "Set2"),10)[1:ny]
    ## if(ny==2) klrpal <- RColorBrewer::brewer.pal(n=3, name = "RdBu")[c(3,1)]
    names(klrpal) <- levels
    klrpal
    
    
    plotly::layout(matrix(1:6,1,6), widths=c(1,0.5,1,0.5,1,0.5))
    k=1
    for(k in 1:3) {

        W <- res$W[[k]]        
        par(mar=c(5,8*cex,4,0), mgp=c(2.2,0.8,0) )
        barplot( t(W), horiz=TRUE, las=1,
                border=NA, col = klrpal,
                names.arg = sub(".*:","",rownames(W)),
                xlim = c(0,1.1)*max(rowSums(W,na.rm=TRUE)),
                cex.axis=1*cex, cex.names=1.1*cex, 
                cex.lab=1*cex, xlab="importance")       
        title( names(res$loadings)[k], cex.main=1.3*cex,
              adj=0.33, xpd=NA)
        legend("topright", legend=names(klrpal),
               cex=1.1*cex, pch=15, col=klrpal, ## fill=klrpal, 
               y.intersp=0.85, inset=c(0.15,0.03))
        
        if(k<99) {
            ## add correlation lines
            par(mar=c(5,0,4,0))
            ##frame()
            plot(0,type="n",xlim=c(0,1),ylim=c(0,nrow(W)),
                 xaxt="n", yaxt="n", bty="n", xlab="")                    
            
            g1 <- rownames(res$W[[k]])
            g2 <- rownames(res$W[[ifelse(k<3,k+1,1)]])
            length(g1)
            length(g2)

            sel <- which(res$edges[,"from"] %in% c(g1,g2) &
                         res$edges[,"to"] %in% c(g1,g2) )
            sel
            ee <- res$edges[sel,]
            ii <- apply(ee[,1:2],1,function(e) which(e%in%g1))
            jj <- apply(ee[,1:2],1,function(e) which(e%in%g2))
            ee$from <- res$edges[sel,][cbind(1:nrow(ee),ii)]
            ee$to   <- res$edges[sel,][cbind(1:nrow(ee),jj)]
            
            ##lwd = exp(abs(ee$rho)/mean(abs(ee$rho)))
            lwd <- ee$importance
            ##lwd = exp(abs(ewt)/mean(abs(ewt)))
            lwd <- rank(abs(lwd),na.last="keep")**1.5
            lwd <- 3.0 * cex * (lwd/max(lwd))
            lty <- 1 + 1*(sign(ee$rho)<0)
            xy <- cbind(match(ee$from,g1), match(ee$to,g2))            
            xy[,2] <- (xy[,2]-0.5) / length(g2) * length(g1)
            klr <- rep(psych::alpha("grey70",0.3), nrow(ee))
            if(showloops) {
                klr <- rep(psych::alpha("grey70",0.2), nrow(ee))
                klr[which(ee$looping)] <- psych::alpha("red3",0.3)
            }
            table(klr)            
            segments(0, xy[,1]-0.5, 1, xy[,2], lwd=lwd, col=klr, lty=lty)
            rr <- paste(round(range(abs(ee$rho)),2),collapse=",")
            ##title(sub=paste0("|rho| = [",rr,"]"),line=-1)
            title(sub=paste0("[",rr,"]"),line=-1.2,cex.sub=cex)
        }
    }
    
}

pgx.makeTriSystemGraph <- function(data, Y, nfeat=25, numedge=100, posonly=FALSE)
{
    if(is.null(names(data))) stop("X.list must be named")   
    if(!all(sapply(data,ncol)==length(Y))) {
        stop("data columns must match Y")
    }
    if(is.null(posonly)) posonly <- FALSE
    
    .getLoadingImportance <- function(k,res.pls,nfeat) {
        U <- res.pls$loadings[[k]]
        V <- res.pls$variates[[k]]
        dim(U)
        U <- U[which(rowSums(abs(U))>0),]
        U <- U[order(-rowSums(U*U)),]
        y <- res.pls$Y        
        vm <- apply(V, 2, function(x) tapply(x,y,mean))
        R0  <- U %*% t(V)
        R <- U %*% t(vm)
        ##R <- Matrix::head(R[order(-rowSums(R**2)),],80)
        R <- Matrix::head(R[order(-rowSums(R**2)),],nfeat)
        ##rownames(R) <- sub(".*:","",rownames(R))        
        W <-  exp(1.5*R/sd(R))
        ##W <-  pmax(R,0)
        W <- W / (1e-5+rowSums(W))
        W <- W * rowSums(R*R,na.rm=TRUE)**0.5
        return(W)        
    }

    .detectSelfLoops <- function(gr,posonly) {                       
        v1 = grep("1:",igraph::V(gr)$name,value=TRUE)
        v3 = grep("3:",igraph::V(gr)$name,value=TRUE)
        
        wt <- (igraph::E(gr)$rho * igraph::E(gr)$importance)
        wt <- wt / max(abs(wt))
        if(posonly) wt <- pmax(wt,0)
        wt <- -log(pmin(1e-8+abs(wt),1))
        
        self.epath <- vector("list",length(v1))
        self.vpath <- vector("list",length(v1))
        self.dist  <- rep(NA,length(v1))
        names(self.vpath) = names(self.epath) = names(self.dist) = v1
        i=1
        for(i in 1:length(v1)) {
            
            d1 = igraph::distances(gr, v1[i], to=v3, mode="out",weights=wt )
            d3 = igraph::distances(gr, v1[i], to=v3, mode="in",weights=wt )
            suppressWarnings( sp1 <- igraph::shortest_paths(
                                  gr, v1[i], to=v3, weights=wt,
                                  mode="out", output="both"))
            suppressWarnings( sp3 <- igraph::shortest_paths(
                                  gr, v1[i], to=v3, weights=wt,
                                  mode="in", output="both"))
            dd = d1 + d3
            j = which.min(dd)
                self.dist[i] <- dd[j]
            self.vpath[[i]] <- NA
            self.epath[[i]] <- NA
            if(!is.infinite(dd[j])) {
                    self.vpath[[i]] <- c(sp1$vpath[[j]]$name,v1[i])
                    self.epath[[i]] <- c(sp1$epath[[j]],sp3$epath[[j]])
            }
        }        
        jj = order(self.dist)
        res <- list( dist=self.dist[jj], vpath=self.vpath[jj],
                        epath=self.epath[jj])
        return(res)
    }

    ##nfeat=50;numedge=100;posonly=FALSE
    .makeTriGraph <- function(res.pls, data, nfeat, numedge, posonly) {
        
        cat("<makeTriSystemGraph:.makeTriGraph> called\n")    
        W.list <- list()
        X.list <- list()        
        k=1
        for(k in 1:3) {
            W.list[[k]] <- .getLoadingImportance(k, res=res.pls, nfeat=nfeat)
            gg <- rownames(W.list[[k]])
            X.list[[k]] <- data[[k]][gg,]
            rownames(X.list[[k]]) <- rownames(W.list[[k]])
        }
        names(W.list) = names(X.list) = names(res.pls$X)
        
        cat("<makeTriSystemGraph:.makeTriGraph> 1\n")
        cat("<makeTriSystemGraph:.makeTriGraph> posonly=",posonly,"\n")
        
        ##numedge=40
        edge.list <- c()    
        k=1
        for(k in 1:3) {
            ## add correlation lines
            p1 <- rownames(W.list[[k]])
            m <- ifelse( k<3, k+1, 1)
            p2 <- rownames(W.list[[m]])
            rho <- stats::cor(res.pls$X[[k]][,p1], res.pls$X[[m]][,p2])  
            ##rho <- cov( res.pls$X[[k]][,p1], res.pls$X[[m]][,p2])  ## COV better?
            if(posonly) rho <- pmax(rho,0)
            q0 = -1
            if(numedge>0) q0 <- Matrix::tail(sort(abs(rho)),numedge)[1]
            idx <- which(abs(rho) > q0, arr.ind=TRUE)                    
            pp <- NULL
            ##if(k<3) 
            pp <- data.frame( from=p1[idx[,1]], to=p2[idx[,2]], rho=rho[idx])   
            if(1) {
                jj <- which(sub("[0-9]:","",pp$from)==sub("[0-9]:","",pp$to))
                if(length(jj)) {
                    pp$rho[jj] <- 0.01  ## no self loops across levels
                }
            }
            self.edges <- data.frame(from=p1, to=p1, rho=0)
            edge.list[[k]] <- rbind(pp, self.edges)
        }
        
        cat("<makeTriSystemGraph:makeGraph> making graph object\n")
        
        ee <- do.call(rbind, edge.list)
        gr <- igraph::graph_from_edgelist(as.matrix(ee[,1:2]), directed=TRUE)
        ww <- lapply(W.list,function(w) rowMeans(w*w)**0.5)  
        names(ww)=NULL; ww=unlist(ww)
        igraph::V(gr)$importance <- ww[igraph::V(gr)$name]  ## node importance  
        
        ## set importance as edge weights
        igraph::E(gr)$rho <- ee$rho    
        igraph::E(gr)$importance <- abs(ee$rho) * sqrt(ww[ee$from]*ww[ee$to])
        
        ##gr <- igraph::simplify(gr)
        gr <- igraph::delete_edges(gr, igraph::E(gr)[which_loop(gr)])
        edges <- data.frame(igraph::get.edgelist(gr), rho=E(gr)$rho,
                            importance=E(gr)$importance )
        colnames(edges)[1:2] <- c("from","to")    
        Matrix::head(edges)
        
        gr
        loops <- .detectSelfLoops(gr,posonly)    
        top.loops <- Matrix::head(loops$epath,10)
        top.loops
        edges$looping <- (igraph::E(gr) %in% unlist(top.loops))
        igraph::E(gr)$looping <- (igraph::E(gr) %in% unlist(top.loops))
        cat("<makeTriSystemGraph:makeGraph> done!\n")
    
        names(W.list) <- names(res.pls$X)
        out <- list(graph=gr, edges=edges, W=W.list)
        return(out)
    }

    ##----------------------------------------------------------------------
    ##----------------------------------------------------------------------
    ##----------------------------------------------------------------------
    
    ## prepend index
    for(i in 1:length(data)) {
        rownames(data[[i]]) <- paste0(i,":",rownames(data[[i]]))
    }

    ## set number of components and features
    NCOMP = 3
    ## NCOMP = length(table(Y))+1  ## number of classes+1
    ##nfeat = 25
    nfeat1 = min(nfeat, nrow(data[[1]]))
    nfeat2 = min(nfeat, nrow(data[[2]]))
    nfeat3 = min(nfeat, nrow(data[[3]]))    
    ##list.keepX = list(rep(40, 3), rep(40,3), rep(40,3))
    list.keepX = list(rep(nfeat1, NCOMP), rep(nfeat2,NCOMP), rep(nfeat3,NCOMP))
    names(list.keepX) <- names(data)
    
    ## set up a full design where every block is connected 
    design = matrix(1, ncol = length(data), nrow = length(data),
                    dimnames = list(names(data), names(data)))
    diag(design) =  0
    design[1,3] = 0.5
    design[3,1] = 0.5
    design      

    
    cat("<makeTriSystemGraph> calling block.splsda!\n")
    res.pls <- mixOmics::block.splsda(
        X = lapply(data,t), Y = Y, ncomp = NCOMP,
        keepX = list.keepX, design = design)    

    cat("<makeTriSystemGraph> calling .makeTriGraph\n")
    ##res.gr <- .makeTriGraph(res.pls, data, nfeat=25, numedge=50)
    res.gr <- .makeTriGraph(
        res.pls, data, nfeat=nfeat, numedge=numedge, posonly=posonly)
    
    res <- c(res.pls, res.gr)
    attr(res,"class") <- attr(res.pls,"class")    
    return(res)
}


##----------------------------------------------------------------------
## Variable importance functions
##----------------------------------------------------------------------

pgx.survivalVariableImportance <-
    function(X, time, status,
             methods=c("glmnet","randomforest","boruta","xgboost","pls"))
{
    ##
    ## multi-class version
    ##
    ##

    imp <- list()
    xnames <- rownames(X)
    sdx <- apply(X,1,sd)    
    if(nrow(X)==1) X <- rbind(X,X)
    
    if(class(status)!="logical" && all(status %in% c(0,1,NA))) {
        stop("status must be logical or 0/1")
    }
    
    y <- survival::Surv(time+0.0001, status)
    
    if("glmnet" %in% methods) {
        
        fam <- "cox"
        
        NFOLD=5
        out0 <- glmnet::cv.glmnet( t(X), y, alpha=0, family=fam,standardize=TRUE, nfolds=NFOLD)
        cf0 <- glmnet::coef.glmnet(out0, s="lambda.min")[,1]
        ##imp[["coxnet.s0"]] <- cf0 * sdx[names(cf0)]
        
        out1 <- glmnet::cv.glmnet( t(X), y, alpha=1, family=fam,standardize=TRUE, nfolds=NFOLD)
        cf1 <- glmnet::coef.glmnet(out1, s="lambda.min")[,1]
        ##imp[["coxnet.s1"]] <- cf1 * sdx[names(cf1)]
        
        out0a <- glmnet::cv.glmnet( t(X), y, alpha=0, family=fam, standardize=FALSE, nfolds=NFOLD)
        cf0a <- glmnet::coef.glmnet(out0a, s="lambda.min")[,1]
        ##imp[["coxnet.a0"]] <- cf0a * sdx[names(cf0a)]
        
        out1a <- glmnet::cv.glmnet( t(X), y, alpha=1, family=fam, standardize=FALSE, nfolds=NFOLD)
        cf1a <- glmnet::coef.glmnet(out1a, s="lambda.min")[,1]
        ##imp[["coxnet.a1"]] <- cf1a * sdx[names(cf1a)]

        imp[["coxnet.a0"]] <- (cf0/max(abs(cf0)) + cf0a/max(abs(cf0a))) * sdx[names(cf0)]
        imp[["coxnet.a1"]] <- (cf1/max(abs(cf1)) + cf1a/max(abs(cf1a))) * sdx[names(cf1)]
        
        
    }
    
    if("randomforest" %in% methods) {
        ##install.packages("caret")
        ##install.packages("randomForest")
        ##install.packages("randomForestSRC")
        ##
        df <- data.frame(time=time, status=status, t(X))
        fit_rf <- randomForestSRC::rfsrc( survival::Surv(time,status) ~ ., data=df)
        vimp <- randomForestSRC::vimp(fit_rf)$importance
        imp[["randomForest"]] <- vimp

        if(0) {
            
            df <- data.frame(time=time+0.01, status=status,t(X))
            tfit = rpart::rpart(survival::Surv(time,status) ~ ., data=df)
            plot(tfit); text(tfit)
            
            ##install.packages("partykit")
            
            (tfit2 <- partykit::as.party(tfit))
            plot(tfit2)

            (tfit3 <- partykit::as.party(rpart::prune(tfit,cp=0.04)))
            plot(tfit3)

            (tfit3 <- partykit::as.party(rpart::prune(tfit,cp=0.05)))
            plot(tfit3)
            
        }

    }

    if("boruta" %in% methods) {
        ##install.packages("Boruta")
        
        imp4 <- rep(0,nrow(X))
        niter=4
        for(k in 1:niter) {
            jj <- sample(ncol(X),ncol(X)*0.9)
            out3 <- Boruta(t(X[,jj,drop=FALSE]), y[jj])
            fd <- factor(out3$finalDecision, levels=c("Rejected","Tentative","Confirmed"))
            fd <- (as.integer(fd)-1)/2
            imp4 <- imp4 + fd/niter
        }
        table(imp4)
        names(imp4) <- rownames(X)
        imp[["Boruta"]] <- imp4
    }
    
    if("xgboost" %in% methods) {
        ##install.packages("xgboost")
        require(xgboost)
        yy <- ifelse( !status, -time, time)
        bst <- xgboost::xgboost(data = t(X), label = yy, booster = "gbtree", 
                                max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
                                verbose=0, objective = "survival:cox")
        ##pred <- predict(bst, t(X))
        ##table(round(pred),y)
        xgmat <- xgboost::xgb.importance(model = bst)
        imp5 <- xgmat$Gain ** 0.33
        imp5 <- imp5[match(rownames(X),xgmat$Feature)]
        names(imp5) <- rownames(X)
        imp5[which(is.na(imp5))] <- 0
        Matrix::tail(sort(imp5))
        imp[["xgboost"]] <- imp5

        ## linear model
        bst2 <- xgboost::xgboost(data = t(X), label = yy, booster = "gblinear", 
                                 max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
                                 verbose=0, objective = "survival:cox")
        xgmat <- xgboost::xgb.importance(model = bst2)
        imp6 <- xgmat$Weight
        names(imp6) <- xgmat$Feature
        imp6 <- imp6[match(rownames(X),names(imp6))]
        Matrix::tail(sort(imp6))
        imp[["xgboost.lin"]] <- imp6
    }

    if("pls" %in% methods) {
        
        res <- pls::plsRcox(t(X),time=time,event=status,nt=5)
        summary(res)
        cf <- res$Coeffs[,1]
        cf[is.na(cf)] <- 0
        cf <- cf * sdx[names(cf)]  ## really?
        imp[["pls.cox"]] <- abs(cf) / max(abs(cf),na.rm=TRUE)
    }

    
    P <- do.call(cbind, imp)
    P <- abs(P)  ## always positive??
    P[is.na(P)] <- 0
    dim(P)
    P <- P[xnames,,drop=FALSE]
    return(P)
}

pgx.multiclassVariableImportance <-
    function(X, y, methods=c("glmnet","randomforest","boruta","xgboost","pls"))
{
    ##
    ## multi-class version
    ##
    ##
    
    imp <- list()
    xnames <- rownames(X)

    dbg("[pgx.multiclassVariableImportance] 0: dim.X = ",dim(X))
    
    if(nrow(X)==1) X <- rbind(X,X)
    sdx <- apply(X,1,sd)    

    ## convert to factor
    y <- factor(y)
    table(y)

    dbg("[pgx.multiclassVariableImportance] 1: dim.X = ",dim(X))
    
    ## resample to minimum size to balance groups
    if(1) {
        NFOLD=5
        ##NSIZE = max(max(table(y)),20)
        NSIZE = 20
        jj <- unlist(tapply(1:length(y), y, function(ii) {
            if(length(ii)<NSIZE) ii <- sample(ii,NSIZE,replace=TRUE)
            return(ii)
        }))
        cat("pgx.multiclassVariableImportance:: augmenting data from",
            length(y),"to",length(jj),"samples\n")
        X <- X[,jj]
        y <- y[jj]
    }


    dbg("[pgx.multiclassVariableImportance] 2: dim.X = ",dim(X))
    
    if("glmnet" %in% methods) {
        require(glmnet)
        fam <- "multinomial"
        ##fam <- "binomial"
        ##if(length(unique(y))>2) fam="multinomial"
        out0 <- glmnet::cv.glmnet( t(X), y, alpha=0, family=fam,standardize=TRUE, nfold=NFOLD)
        cf0 <- Matrix::rowMeans(do.call(cbind, coef(out0, s="lambda.min"))[-1,]**2)
        
        out1 <- glmnet::cv.glmnet( t(X), y, alpha=1, family=fam,standardize=TRUE, nfold=NFOLD)
        cf1 <- Matrix::rowMeans(do.call(cbind, coef(out1, s="lambda.min"))[-1,]**2)        
        
        out0a <- glmnet::cv.glmnet( t(X), y, alpha=0, family=fam, standardize=FALSE, nfold=NFOLD)
        cf0a <- Matrix::rowMeans(do.call(cbind, coef(out0a, s="lambda.min"))[-1,]**2)
        
        out1a <- glmnet::cv.glmnet( t(X), y, alpha=1, family=fam, standardize=FALSE, nfold=NFOLD)
        cf1a <- Matrix::rowMeans(do.call(cbind, coef(out1a, s="lambda.min"))[-1,]**2)

        imp[["glmnet.a0"]] <- (cf0/max(abs(cf0)) + cf0a/max(abs(cf0a))) * sdx[names(cf0)]
        imp[["glmnet.a1"]] <- (cf1/max(abs(cf1)) + cf1a/max(abs(cf1a))) * sdx[names(cf1)]

    }
    
    if("randomforest" %in% methods) {
        ##install.packages("caret")
        ##install.packages("randomForest")
        ##
        require(randomForest)
        fit_rf = randomForest::randomForest(t(X), factor(y))
        imp[["randomForest"]] <- fit_rf$importance[,1]
    }

    if("boruta" %in% methods) {
        ##install.packages("Boruta")
        ##require(Boruta)        
        imp4 <- rep(0,nrow(X))
        niter=4
        for(k in 1:niter) {
            jj <- sample(ncol(X),ncol(X)*0.9)
            out3 <- Boruta::Boruta(t(X[,jj,drop=FALSE]), y[jj])
            fd <- factor(out3$finalDecision,
                         levels=c("Rejected","Tentative","Confirmed"))
            fd <- (as.integer(fd)-1)/2
            imp4 <- imp4 + fd/niter
        }
        table(imp4)
        names(imp4) <- rownames(X)
        imp[["Boruta"]] <- imp4
    }
    
    if("xgboost" %in% methods) {
        ##install.packages("xgboost")
        
        ny <- length(table(y))
        yy <- as.integer(factor(y))-1
        bst <- xgboost::xgboost(data = t(X), label = yy, booster = "gbtree", 
                                max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
                                num_class = ny,
                                verbose=0, objective = "multi:softmax")
        ##pred <- predict(bst, t(X))
        ##table(round(pred),y)
        xgmat <- xgboost::xgb.importance(model = bst)
        imp5 <- xgmat$Gain ** 0.2
        imp5 <- imp5[match(rownames(X),xgmat$Feature)]
        names(imp5) <- rownames(X)
        imp[["xgboost"]] <- imp5

        ## linear model
        bst2 <- xgboost::xgboost(data = t(X), label = yy, booster = "gblinear", 
                                 max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
                                 num_class = ny,
                                 verbose=0, objective = "multi:softmax")
        xgmat <- xgboost::xgb.importance(model = bst2)
        imp6 <- xgmat$Weight
        names(imp6) <- xgmat$Feature
        imp6 <- imp6[match(rownames(X),names(imp6))]
        imp[["xgboost.lin"]] <- imp6
    }

    if("pls" %in% methods) {
        
        n <- min(25,nrow(X))
        colnames(X) <- names(y) <- paste0("sample",1:length(y))
        res <- mixOmics::splsda( t(X), y, keepX = c(n,n))  
        impx <- rowMeans(res$loadings$X[rownames(X),]**2)
        imp[["spls.da"]] <- impx
    }
    
    P <- do.call(cbind, imp)
    P <- abs(P)  ## always positive??
    P[is.na(P)] <- 0
    dim(P)
    P <- P[xnames,,drop=FALSE]
    return(P)
}


##methods=c("glmnet","randomforest","boruta","xgboost")
##methods=c("glmnet","randomforest","xgboost")
pgx.variableImportance <-
    function(X, y, methods=c("glmnet","randomforest","boruta","xgboost","pls"))

{
    if(length(unique(y[!is.na(y)]))!=2) {
        stop("pgx.variableImportance:: only binary outcome variables")
    }
    
    imp <- list()
    xnames <- rownames(X)
    if(nrow(X)==1) X <- rbind(X,X)
    sdx <- apply(X,1,sd)    

    ## convert to factor
    y <- factor(y)
    table(y)

    ## resample to minimum size to balance groups
    if(1) {
        NFOLD=5
        ##NSIZE = max(max(table(y)),20)
        NSIZE = 20
        jj <- unlist(tapply(1:length(y), y, function(ii) {
            if(length(ii)<NSIZE) ii <- sample(ii,NSIZE,replace=TRUE)
            return(ii)
        }))
        cat("pgx.variableImportance:: augmenting data from",
            length(y),"to",length(jj),"\n")
        X <- X[,jj]
        y <- y[jj]
    }    
    
    if("glmnet" %in% methods) {
        
        out0 <- glmnet::cv.glmnet( t(X), y, alpha=0, family="binomial",
                          standardize=TRUE)
        cf0 <- coef(out0, s="lambda.min")[-1,1]
        ##imp[["glmnet.s0"]] <- abs(cf0) * sdx[names(cf0)]
        
        out1 <- glmnet::cv.glmnet( t(X), y, alpha=1, family="binomial",
                          standardize=TRUE)
        cf1 <- coef(out1, s="lambda.min")[-1,1]
        ##imp[["glmnet.s1"]] <- abs(cf1) * sdx[names(cf1)]
        
        out0a <- glmnet::cv.glmnet( t(X), y, alpha=0, family="binomial",
                           standardize=FALSE)
        cf0a <- coef(out0a, s="lambda.min")[-1,1]
        ##imp[["glmnet.a0"]] <- abs(cf0a) * sdx[names(cf0a)]
        
        out1a <- glmnet::cv.glmnet( t(X), y, alpha=1, family="binomial",
                           standardize=FALSE)
        cf1a <- coef(out1a, s="lambda.min")[-1,1]
        ##imp[["glmnet.a1"]] <- abs(cf1a) * sdx[names(cf1a)]

        imp[["glmnet.a0"]] <- (cf0/max(abs(cf0)) + cf0a/max(abs(cf0a))) * sdx[names(cf0)]
        imp[["glmnet.a1"]] <- (cf1/max(abs(cf1)) + cf1a/max(abs(cf1a))) * sdx[names(cf1)]
    }
    
    if("randomforest" %in% methods) {
        ##install.packages("caret")
        ##install.packages("randomForest")
        ##
        
        fit_rf = randomForest::randomForest(t(X), factor(y))
        imp[["randomForest"]] <- fit_rf$importance[,1]
    }

    if("boruta" %in% methods) {
        ##install.packages("Boruta")
        
        imp4 <- rep(0,nrow(X))
        niter=4
        for(k in 1:niter) {
            jj <- sample(ncol(X),ncol(X)*0.9)
            out3 <- Boruta::Boruta(t(X[,jj,drop=FALSE]), y[jj])
            fd <- factor(out3$finalDecision, levels=c("Rejected","Tentative","Confirmed"))
            fd <- (as.integer(fd)-1)/2
            imp4 <- imp4 + fd/niter
        }
        table(imp4)
        names(imp4) <- rownames(X)
        imp[["Boruta"]] <- imp4
    }
    
    if("xgboost" %in% methods) {
        ##install.packages("xgboost")
        
        y1 <- as.integer(factor(y))-1
        bst <- xgboost::xgboost(data = t(X), label = y1, booster = "gbtree", 
                                max_depth = 2, eta = 1, nthread = 2, nrounds = 2, 
                                verbose=0, objective = "binary:logistic")
        ##pred <- predict(bst, t(X))
        ##table(round(pred),y)
        xgmat <- xgboost::xgb.importance(model = bst)
        imp5 <- xgmat$Gain ** 0.2
        imp5 <- imp5[match(rownames(X),xgmat$Feature)]
        names(imp5) <- rownames(X)
        imp[["xgboost"]] <- imp5

        bst2 <- xgboost::xgboost(data = t(X), label = y1, booster = "gblinear", 
                                 eta = 0.3, nthread = 1, nrounds = 20,
                                 verbose=0, objective = "binary:logistic")
        xgmat <- xgboost::xgb.importance(model = bst2)
        imp6 <- xgmat$Weight
        names(imp6) <- xgmat$Feature
        imp6 <- imp6[match(rownames(X),names(imp6))]
        imp[["xgboost.lin"]] <- imp6
    }

    if("pls" %in% methods) {
        ##install.packages("caret")
        ##install.packages("randomForest")
        ##
        
        n <- min(25,nrow(X))
        colnames(X) <- names(y) <- paste0("sample",1:length(y))
        res <- mixOmics::splsda( t(X), y, keepX = c(n,n))  
        ##impx <- res$loadings$X[rownames(X),1]
        impx <- rowMeans(res$loadings$X[rownames(X),]**2)
        imp[["spls.da"]] <- impx
    }
    
    P <- do.call(cbind, imp)
    P <- abs(P)  ## always positive??
    P[is.na(P)] <- 0
    dim(P)
    P <- P[xnames,,drop=FALSE]
    return(P)
}


##----------------------------------------------------------------------
## combine all importance values
##----------------------------------------------------------------------

if(0) {
    
    
    source("../R/gx-heatmap.r")
    source("../R/gx-limma.r")    
    load("../pgx/GSE10846-dlbcl-mRNA-8k.pgx",verbose=1)
    
    y <- factor(ngs$samples$dlbcl.type)
    y1 <- 1*(ngs$samples$dlbcl.type=="GCB")
    X <- as.matrix(ngs$X)
    X <- as.matrix(ngs$gsetX)
    mod1 <- model.matrix(~ y)
    X <- ComBat(X, batch=ngs$samples$Chemotherapy, mod=mod1)
    
    ## pre-select using LIMMA
    res <- gx.limma(X, y, fdr=1, lfc=0)
    jj <- which(res[,"adj.P.Val"] < 0.05 & abs(res[,"logFC"]) > 1)
    length(jj)
    X <- X[rownames(res)[jj],,drop=FALSE]
    ##X <- Matrix::head(X[order(-apply(X,1,sd)),],500)
    sdx <- apply(X,1,sd)
    dim(X)
    
    ##----------------------------------------------------------------------
    ## compute importance values
    ##----------------------------------------------------------------------
    methods=c("glmnet","randomforest","boruta","xgboost")
    methods=c("glmnet","randomforest","xgboost")
    P <- pgx.variableImportance(X, y, methods=methods)
    P <- t( t(P) / apply(P,2,max,na.rm=TRUE))
    P <- pmax(P,0.1)
    P <- P[order(-rowSums(P,na.rm=TRUE)),,drop=FALSE]
    Matrix::head(P)
    
    par(mfrow=c(2,2), mar=c(8,4,2,2))
    frame()
    barplot( t(Matrix::head(P,30)), las=1, horiz=TRUE)
    klr <- grey.colors(ncol(P))
    legend("topright", legend=colnames(P), fill=klr)
    
    R <- (apply(P,2,rank)/nrow(P))**4
    R <- R[order(-rowSums(R)),,drop=FALSE]
    frame()
    barplot( t(Matrix::head(R,30)), las=1, horiz=TRUE)
    legend("topright", legend=colnames(R), fill=klr)
    
    par(mfrow=c(4,4), mar=c(4,4,2,2))
    for(i in 1:16) {
        g <- rownames(R)[i]
        boxplot( X[g,,drop=FALSE] ~ y, col="grey90", main=g)
    }
    Matrix::head(res[rownames(R),])

    ##install.packages("rpart.plot")
    
    
    sel <- Matrix::head(rownames(R),100)
    sel <- Matrix::head(rownames(R),20)
    df <- data.frame( y=y, t(X[sel,,drop=FALSE] ))
    rf <- rpart::rpart( y ~ ., data=df)
    par(mfrow=c(1,1), mar=c(4,4,2,2)*0)
    rpart.plot(rf)

}
