##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

ngs.hiveplot <- function(ngs, pheno, level=1, ntop=400, main="", axis.lab=c("GX","CN","ME"),
                         bkgnd="white", rx=2.2, cex=1, slen=-1)
{
    
    
    res <- omx$zstats[[level]]
    res <- cbind(res, omx$stats[[pheno]][[level]])
    res <- res[!grepl("^PHENO",rownames(res)),]
    rownames(res) <- gsub("^.*:","",rownames(res)) ## strip prefix
    if(slen>0) rownames(res) <- substr(rownames(res),1,slen)
    ann.cex <- cex*0.7

    hpd <- omx.makeHivePlotData_( res, rho.min=0.15, ntop=ntop, rx=rx,
                                 cxi=(0.05+ann.cex/7) )
    write.csv(hpd$nodes.ann, file="/tmp/annode.csv", row.names=FALSE)
    grid::grid.newpage()
    axlab.col <- ifelse(bkgnd=="black", "grey90","grey15")
    plotHive( hpd, np=FALSE, ch=1.4, bkgnd = bkgnd,
             axLabs = c(paste0(main,"\n",axis.lab[1]),axis.lab[2],axis.lab[3]),
             axLab.pos = c(0.4,0.33,0.33),
             axLab.gpar = grid::gpar(col=axlab.col,cex=cex*1.3),
             anNodes = "/tmp/annode.csv",
             anNode.gpar = grid::gpar(cex=0.7*cex,col=axlab.col,lwd=0.50))
}

omx.makeHivePlotData_ <- function( res, rho.min=0.15, cxi=0.11, use.alpha=TRUE,
                                 ntop=1000, ew=4, nv=10, dx.dir=c(-1,1), rx=2.2 )
{
    ##rho.min=0.15;cxi=0.13;use.alpha=TRUE;ew=4;nv=10;ntop=1000;dx.dir=c(-1,1)
    ## omics score
    res <- res[order(-abs(res[,"score"])),]
    if(ntop>0) res <- Matrix::head(res, ntop)

    ## stuff...
    dx <- res[,"gx.rho"]
    dc <- res[,"cn.rho"]
    dm <- res[,"me.rho"]
    n <- nrow(res)
    names(dx) <- paste("x",rownames(res),sep="")
    names(dc) <- paste("c",rownames(res),sep="")
    names(dm) <- paste("m",rownames(res),sep="")
    dx[is.na(dx)] <- 0
    dc[is.na(dc)] <- 0
    dm[is.na(dm)] <- 0

    ## define node properties
    dx0 = dc0 = dm0 = rho.min ## threshold for coloring
    dx0 = dc0 = dm0 = 0 ## threshold for coloring
    kr0 <- c("grey40","green4","red")[1+1*(dx < -dx0)+2*(dx > dx0)]
    n0 <- data.frame( id=1:n,lab=names(dx), axis=1,
                     radius=2*(rank(dx,na.last="keep")/n-0.5), value=dx,
                     size=abs(dx)/max(abs(dx),na.rm=TRUE), color=kr0 )
    kr1 <- c("grey40","blue","orange2")[1+1*(dc < -dc0) + 2*(dc > dc0)]
    n1 <- data.frame( id=(n+1):(2*n), lab=names(dc), axis=2,
                     radius=2*(rank(dc,na.last="keep")/n-0.5), value=dc,
                     size=abs(dc)/max(abs(dc),na.rm=TRUE), color=kr1 )
    kr2 <- c("grey40","purple","gold")[1+1*(dm < -dm0)+2*(dm > dm0)]
    n2 <- data.frame( id=(2*n+1):(3*n), lab=names(dm), axis=3,
                     radius=2*(rank(dm,na.last="keep")/n-0.5), value=dm,
                     size=abs(dm)/max(abs(dm),na.rm=TRUE), color=kr2 )
    nn <- rbind(n0,n1,n2)
    nn$size <- (0.10 + 0.60 * nn$size*2)*2
##  nn$size <- 0.60
    nn$lab   <- as.character(nn$lab)
    nn$axis  <- as.integer(nn$axis)
    nn$color <- as.character(nn$color)

    ## define edge properties
    ee1 <- data.frame(n1=names(dx), n2=names(dc), cor=res[,"cn.tau"])
    ee2 <- data.frame(n1=names(dx), n2=names(dm), cor=res[,"me.tau"])
    ##ee3 <- data.frame(n1=names(dm), n2=names(dc), cor=res[,"tau.cn.me"])
    ee  <- rbind(ee1,ee2)
    ee  <- ee[which(abs(ee$cor) > rho.min),]
    ee$id1 <- nn$id[match(ee$n1,nn$lab)]
    ee$id2 <- nn$id[match(ee$n2,nn$lab)]
    ee$weight <- (abs(ee$cor) - rho.min) / (1 - rho.min)
##  ee$weight <- 4
    if(use.alpha) {
        alpha <- abs(ee$weight * nn[ee$id1,"radius"] * nn[ee$id2,"radius"])
        alpha <- alpha - min(alpha)
        ##alpha <- pmax( (alpha/max(alpha))**2.0, 0.10)
        alpha <- pmax( (alpha/max(alpha))**1.0, 0.20)
    } else {
        alpha <- 1
    }
    ee$weight <- 0.0 + ew*(ee$weight)**1.5  ## add minimum and gamma
    ee$alpha <- alpha
    ##ee$color  <- c("steelblue","greenyellow")[1+1*(ee$cor>0)]
    ee$color  <- c("steelblue","darkorange")[1+1*(ee$cor>0)]
    cc2 <- cbind(t(col2rgb(ee$color)/255), alpha)
    ee$color  <- apply(cc2[,],1,function(x)rgb(x[1],x[2],x[3],x[4]))
    ee <- ee[!duplicated(ee),]
    ee <- ee[ which(as.character(ee$n1)!=as.character(ee$n2)), ]

    ## node annotation matrix
    rr0 <- res[which(dx>0 & abs(res[,"score"])>0),]  ## upregulated
    rr1 <- res[which(dx<0 & abs(res[,"score"])>0),]  ## downregulated
    rr0 <- Matrix::head(rr0, nv)
    rr1 <- Matrix::head(rr1, nv)
    rr0 <- rr0[order(-rr0[,"gx.rho"]),]
    rr1 <- rr1[order(+rr1[,"gx.rho"]),]
    yy0 <- seq(1,1-nv*cxi,-cxi)[1:nrow(rr0)]
    yy1 <- seq(-1,-1+nv*cxi,cxi)[1:nrow(rr1)] + 0.5
    ##rx <- 2.2 ## distance from center
    rownames(rr0) <- paste("x",rownames(rr0),sep="")
    rownames(rr1) <- paste("x",rownames(rr1),sep="")
    yy0 <- yy0 - nn[rownames(rr0),]$radius + 0.2
    yy1 <- yy1 - nn[rownames(rr1),]$radius - 0.2  ## offset from top
    M0 <- data.frame( node.lab=rownames(rr0),
                     node.text=sub("^x","",rownames(rr0)),
                     angle=atan(yy0/rx)/pi*180, radius=sqrt(rx**2+yy0**2),
                     offset=0, hjust=0, vjust=0)
    M1 <- data.frame( node.lab=rownames(rr1),
                     node.text=sub("^x","",rownames(rr1)),
                     angle=180-atan(yy1/rx)/pi*180, radius=sqrt(rx**2+yy1**2),
                     offset=0, hjust=1, vjust=0)
    mm <- c()
    if(  1 %in% dx.dir ) mm <- rbind(mm, M0)
    if( -1 %in% dx.dir ) mm <- rbind(mm, M1)
    rownames(res) <- sub("^x","",rownames(res))

    ## HPD object
    hpd <- c()
    hpd$nodes <- nn
    hpd$edges <- ee
    hpd$type <- "2D"
    hpd$desc <- paste("3 axes --",nrow(nn),"nodes --",nrow(ee),"edges")
    ##    hpd$axis.cols <- c("red","green","blue")
    hpd$axis.cols <- rep("grey40",3)
    hpd$score <- res
    hpd$nodes.ann <- mm
    class(hpd) <- "HivePlotData"
    return(hpd)
}
