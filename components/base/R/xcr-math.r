##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##matlist=list(S1,S2);weights=1
matrix.mean.SAVE <- function( matlist, weights=1) {
    matx <- lapply(matlist,function(x) as.vector(x))
    matx <- as.matrix(do.call(rbind, matx))
    if(length(weights)==1) weights = rep(weights,length(matlist))
    weigths = weights / sum(weights)
    matx = matx * weights
    matx <- colMeans(matx * weights, na.rm=TRUE)
    matrix(matx, nrow=nrow(matlist[[1]]), ncol=ncol(matlist[[1]]))
}

matrix.prod.SAVE <- function( matlist ) {
    matx <- lapply(matlist,function(x) as.vector(x))
    matx <- as.matrix(do.call(cbind, matx))
    nna = rowSums(!is.na(matx))
    matx <- apply(matx, 1, prod, na.rm=TRUE)**(1/nna)
    xp = matrix(matx, nrow=nrow(matlist[[1]]), ncol=ncol(matlist[[1]]))
    return(xp)
}

matrix.mean <- function( ... ) {
    matlist = list(...)
    if(class(matlist[[1]])=="list") matlist = matlist[[1]]
    matnna = lapply(matlist, function(x) 1*!is.na(x))
    p  = Reduce( '+', matnna)
    matlist = lapply(matlist, function(x) {x[is.na(x)]=0;x})
    xsum = Reduce( '+', matlist) / p
    return(xsum)
}

matrix.prod <- function( ... , na.value=1) {
    matlist = list(...)
    if(class(matlist[[1]])=="list") matlist = matlist[[1]]
    p  = Reduce( '+', lapply(matlist, function(x) 1*!is.na(x)))
    i=1
    for(i in 1:length(matlist)) {
        if(any(is.na(matlist[[i]]))) {
            jj <- Matrix::which( is.na(matlist[[i]]), arr.ind=TRUE)
            matlist[[i]][jj] = na.value
        }
    }
    xprod = Reduce( '*', matlist)
    xprod = xprod ** (1/p)
    jj <- Matrix::which(p==0,arr.ind=TRUE)
    if(length(jj)>0) {
        xprod[jj] <- NA
    }
    return(xprod)
}

##na.fill = 1
sparse=NULL
matrix.prodSPARSE <- function( ..., na.fill=1) {
    matlist = list(...)
    if(class(matlist[[1]])=="list") matlist = matlist[[1]]
    if(is.null(sparse)) sparse = is(matlist[[1]], "sparseMatrix")
    p  = Reduce( '+', lapply(matlist, function(x) 1*!is.na(x)))
    i=1
    for(i in 1:length(matlist)) {
        if(any(is.na(matlist[[i]]))) {
            mx = 1
            if(na.fill=="mean") {
                mx = mean(matlist[[i]],na.rm=TRUE)
                ##mx = median(matlist[[i]],na.rm=TRUE)
            }
            jj <- which( is.na(matlist[[i]]), arr.ind=TRUE)
            matlist[[i]][jj] = mx
        }
    }
    xprod = Reduce( '*', matlist)
    if(na.fill!="mean")
        xprod = xprod ** (1/p)
    jj <- which(p==0,arr.ind=TRUE)
    if(sparse) {
        xprod[is.na(xprod)] <- 0
        xprod = Matrix::Matrix(xprod, sparse=TRUE)
    }
    if(length(jj)>0) {
        if(sparse) jj <- jj[which(jj[,2]==1 | jj[,1]==1),]
        xprod[jj] <- NA
    }
    return(xprod)
}


merge.similarities <- function(S.list, sparse=NULL)
{
    length(S.list)
    if(is.null(sparse)) sparse = is(S.list[[1]], "sparseMatrix")
    if( length(S.list)>1 && sparse) {
        ## merge similarities (i.e. multiply) but allow for missing values
        S.list = lapply(S.list, function(m) Matrix::Matrix(m,sparse=TRUE))
        cat("merge similarities...(sparse)\n")
        idx = Matrix::which(S.list[[1]]!=0, arr.ind=TRUE)
        x = S.list[[1]][idx]
        p = rep(1,length(x))
        i=2
        for(i in 2:length(S.list)) {
            not.na <- !is.na(S.list[[i]][,1])
            jj = which(not.na[idx[,1]] & not.na[idx[,2]] )
            x[jj] = x[jj] * S.list[[i]][idx[jj,]]
            p[jj] = p[jj] + 1
        }
        table(p)
        x = x**(1/p)
        table(x!=0)
        jj = which(x!=0)
        n = nrow(S.list[[1]])
        Q = Matrix::sparseMatrix(i=idx[jj,1], j=idx[jj,2], x=x[jj],
                         use.last.ij=TRUE, dims=c(n,n))
        ## set NA??
        jj <- which(Matrix::rowMeans(Q==0)==1)
        if(length(jj)>0) {
            Q[jj,1] <- NA
            Q[1,jj] <- NA
        }

    } else if(length(vars)>1 && !sparse) {
        ## merge similarities (i.e. multiply) but allow for missing values
        cat("merge similarities... (full)\n")
        S.list = lapply(S.list, function(m) Matrix::Matrix(m,sparse=FALSE))
        Q = as.matrix(S.list[[1]])
        p = matrix(1,nrow(Q),ncol(Q))
        p = 1*(!is.na(Q))
        i=2
        for(i in 2:length(S.list))  {
            s1 = as.matrix(S.list[[i]])
            not.na <- !is.na(S.list[[i]][,1])
            j0 = which(!not.na)
            s1[j0,] = 1
            s1[,j0] = 1
            Q = Q * s1
            j1 = which(not.na)
            p[j1,j1] = p[j1,j1] + 1
        }
        Q = Q**(1/p)
    } else {
        Q = S.list[[1]]
    }
    rownames(Q) = colnames(Q) = rownames(S.list[[1]])
    return(Q)
}

tcosine.similarity <- function(X, Y=NULL, method=NULL) {
    if(!is.null(Y)) {
        return( cosine.similarity( t(X), Y=t(Y), method=method))
    }
    return( cosine.similarity( t(X), Y=NULL, method=method))
}

cosine.similarity <- function(X, Y=NULL, method=NULL)
{
    
    ## sparse cosine: A.B / (|A||B|)
    ## handles sparse matrix and missing values
    X <- as(X, "dgCMatrix")
    if(is.null(Y)) {
        ii = which(Matrix::colMeans(is.na(X))!=1)
        zX = X[,ii]
        zX[is.na(zX)] = 0
        s1 <- crossprod( !is.na(X[,ii]), zX**2)
        ss <- (s1 * t(s1))**0.5
        M <- crossprod(zX) / (ss + 1e-20)
        S <- matrix(NA, nrow=ncol(X), ncol=ncol(X))
        S[ii,ii] <- as.matrix(M)
        return(S)
    } else {
        Y <- as(Y, "dgCMatrix")
        ii = which(Matrix::colMeans(is.na(X))!=1)
        jj = which(Matrix::colMeans(is.na(Y))!=1)
        zX = X[,ii]
        zX[is.na(zX)] = 0
        zY = Y[,jj]
        zY[is.na(zY)] = 0
        s1 <- crossprod( !is.na(X[,ii]), zY**2)
        s2 <- crossprod( !is.na(Y[,jj]), zX**2)
        ss <- (s1 * t(s2))**0.5
        M <- crossprod(zX,zY) / (ss + 1e-20)
        S <- matrix(NA, nrow=ncol(X), ncol=ncol(Y))
        S[ii,jj] <- as.matrix(M)
        return(S)
    }
}

##th=0.01;block=100;k=100;ties.method="random"
tcosine.sparse <- function(X, k=100, th=0.01, block=100, ties.method="random",
                           gpu=FALSE ) {
    dim(X)
    X = X / (1e-20 + Matrix::rowSums(X**2)**0.5)
    nblock = ceiling(nrow(X)/block)
    nblock
    idx <- c()
    i=1
    for(i in 1:nblock) {
        j0 <- (i-1)*block + 1
        j1 <- min((i-1)*block + block,nrow(X))
        jj = j0:j1
        ##rho = tcosine.similarity( X[jj,], X[,] )
        if(gpu==TRUE) {
            
            rho <- gpuTcrossprod(X[jj,], X[,])
        } else {
            rho <- tcrossprod( X[jj,], X[,] )
        }
        if(ties.method=="random") {
            idx.j = apply( rho,1,function(x) Matrix::head(order(x,rnorm(length(x)),decreasing=TRUE),k) )
        } else {
            idx.j = apply( rho,1,function(x) Matrix::head(order(x,decreasing=TRUE),k) )
        }
        idx.j = as.vector(idx.j)
        idx.i = as.vector(sapply(1:nrow(rho), rep, k))
        x = rho[cbind(idx.i,idx.j)]
        idx.i = idx.i + j0 - 1
        x1 = cbind( idx.i, idx.j, x=x)
        if(th>0) x1 <- x1[ which(x1[,3] > th),,drop=FALSE]
        if(nrow(x1)>0) {
            ## add symmetric part
            x2 = cbind( x1[,2], x1[,1], x=x1[,3])
            x1 <- rbind( x1, x2)
            idx <- rbind( idx, x1)
        }
    }
    S <- Matrix::sparseMatrix( i=idx[,1], j=idx[,2], x=idx[,3],
                      dims=c(nrow(X),nrow(X)), use.last.ij=TRUE )
    dim(S)
    rownames(S) = colnames(S) = rownames(X)

    sp = round(100*mean(S!=0,na.rm=TRUE),digits=2)
    sp.na = round(100*mean(is.na(S),na.rm=TRUE),digits=2)
    cat("tcosine.sparse: sparsity=",sp,"%\n")
    if(sp.na>0) cat("tcosine.sparse: NA=",sp.na,"%\n")

    return(S)
}


##X=G;blocksize=1000;th=0.05;k=500;mc.cores=10;blocksize=1000
tcosine.sparse.paral <- function(X, th=0.01, k=100, mc.cores=4, blocksize=100,
                                 verbose=1){
    ##
    ## Normalizing X then doing tcrossprod is same as doing
    ## tcosine.similarity on non-normalized X.
    ##
    X = X / (1e-20 + Matrix::rowSums(X**2)**0.5)
    iter = ceiling(nrow(X)/blocksize)
    sim <-c()
    i=1
    for(i in 1:iter){
        down = i*blocksize-blocksize+1
        up = min(i*blocksize,nrow(X))
        if(verbose) print(paste0(down," - ",up))
        j=1
        s <- parallel::mclapply(down:up, function(j) {
            cr = tcrossprod(X[j,], X)[1,];
            jj = as.vector(which(cr>=th))
            jj = Matrix::head(jj[order(cr[jj], decreasing=TRUE)],k)
            if(length(jj)>0) {
                ##data.frame(from=rep(rownames(X)[j],length(cr)),to=names(cr),val=as.numeric(cr))
                data.frame(from=rep(j,length(jj)), to=jj, val=as.numeric(cr[jj]))
            } else {
                data.frame(from=NA,to=NA,val=NA)[0,]
            }
        }, mc.cores = mc.cores)
        sim[down:up] <- s
    }
    ##sim1 <- data.table::rbindlist(sim)
    sim1 <- do.call(rbind, sim)
    dim(sim1)
    S <- Matrix::sparseMatrix( i=sim1[,1], j=sim1[,2], x=sim1[,3])
    if(verbose) {
        sp = round(100*mean(S!=0,na.rm=TRUE),digits=2)
        sp.na = round(100*mean(is.na(S),na.rm=TRUE),digits=2)
        cat("tcosine.sparse: sparsity=",sp,"%\n")
        if(sp.na>0) cat("tcosine.sparse: NA=",sp.na,"%\n")
    }
    return(S)
}

seq.similarity.paral.SAVE <- function(aln, th=0.75, k=500, mc.core=100, blocksize=1000){
    na = which(is.na(aln))
    no.na = which(!is.na(aln))
    y = seq.encode(aln[no.na], fill.na=1, verbose=0)
    ##y = y / (1e-20 + Matrix::rowSums(y,na.rm=TRUE))
    y = y / (1e-20 + Matrix::rowSums(y**2)**0.5)
    print(dim(y))
    rownames(y)=no.na

    iter = ceiling(nrow(y)/blocksize)
    sim <-c()
    for(i in 1:iter){
        down=i*blocksize-blocksize+1
        up=min(i*blocksize,nrow(y))
        print(paste0(down," - ",up))
        s<-mclapply(down:up, function(j){
            cr = tcrossprod(y[j,], y);
            cr = Matrix::head(sort(cr[1,],decreasing = TRUE),k)
            cr = cr[which(cr>=th)];
            if(length(cr)>0) {data.frame(from=rep(rownames(y)[j],length(cr)),to=names(cr),val=as.numeric(cr))} else {data.frame(from="",to="",val="")[0,]}
        }, mc.cores = mc.core)
        sim[down:up] <- s
    }
    if(length(na)>0){
        sim[length(no.na)+1] <-list(data.frame(from=na,to=rep(1,length(na)),val=NA))
    }
    return(sim)
}


cosine.similarity.EXACT <- function(X) {
    vec.cos <- function(a,b) {
        j <- which(!is.na(a) & !is.na(b))
        as.numeric( sum(a[j]*b[j]) / sum(a[j]**2)**0.5 / sum(b[j]**2)**0.5)
    }
    M <- matrix(NA,ncol(X),ncol(X))
    for(i in 1:ncol(X))
        for(j in i:ncol(X))
            M[i,j] = M[j,i] = vec.cos( X[,i], X[,j] )
    ##M[which(is.na(M))] <- NA
    M
}

jaccard_similarity <- function(m)
{
    
    A <- tcrossprod(m)
    im <- which(A > 0, arr.ind=TRUE, useNames = F)
    b <- rowSums(m)
    Aim <- A[im]
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim)
    sp = Matrix::sparseMatrix( i=im[,1], j=im[,2], x=x, dims=dim(A) )
    rownames(sp) <- colnames(sp) <- rownames(m)
    sp
}

length.similarityEXACT <- function(x,p=1)
{
    D = 1.0*abs(outer(x, x, FUN='-'))
    L = 1.0*abs(outer(x, x, FUN='pmax'))
    ##R = (1 - D / L)**p
    R = exp( -p*(D/L))
    return(R)
}

length.similaritySPARSE <- function(x,r=0.1)
{
    X = length.encode(x,r=r)
    S = cosine.similarity(t(X))
    rownames(S) <- colnames(S) <- NULL
    return(S)
}

length.encode <- function(x, r=0.1, a=0.25)
{
    
    x0 = x
    x = x[which(!is.na(x))]
    logx = log(x)
    maxlen = max(logx, na.rm=TRUE)
    minlen = min(logx, na.rm=TRUE)
    c(minlen,maxlen)
    dx = log(1 + a*r)
    ##dx = 0.05  ## 5% intervals
    brks = seq(minlen-dx, maxlen+dx, dx)
    ix = 2 + as.integer(cut(logx, breaks=brks))
    endx = max(ix,na.rm=TRUE)+2
    M0 = model.matrix( ~ 0 + factor(ix, levels=1:endx) )
    M = matrix(NA, nrow=length(x0), ncol=ncol(M0))
    M[which(!is.na(x0)),] <- M0
    M = Matrix::Matrix(M, sparse=TRUE)
    colnames(M) <- NULL
    rownames(M) <- NULL
    if(r==0) {
        S = cosine.similarity(t(M))
        return(S)
    }
    n = ceiling(r/dx)
    n = 2  ## alway 2 for now..
    X = Matrix::Matrix(M)
    i = 1
    for(i in 1:n) {
        dM = cbind( M[,-i:-1], matrix(0,nrow(M),ncol=i))
        X  =  X + dM
    }
    for(i in 1:n) {
        dM = cbind(matrix(0,nrow(M),ncol=i), M[,1:(ncol(M)-i)])
        X  =  X + dM
    }
    return(X)
}

if(0) {
    x[1]=NA
    length.encode(x, r=0.1, a=0.5)
    length.encode(x, r=0.1, a=0.33)
    length.encode(x, r=0.1, a=0.25)

}


##===================================================================================
##===================================================================================
##===================================================================================

clustering.score <- function(xy, labels, kinter=NULL) {
    ## Defined as the ratio between inter-distance between clusters
    ## and intra-distance within the clusters.
    ##
    nlab = length(unique(labels))
    groups = tapply( 1:nrow(xy), labels, list)
    groups

    xy = scale(xy)
    xy[is.na(xy)] = 0
    xy.mean = t(sapply(groups, function(g) colMeans(xy[g,,drop=FALSE])))
    inter.dist <- pos2dist( xy.mean )
    intra.dist <- sapply( groups, function(g) mean(pos2dist(xy[g,,drop=FALSE])))
    intra.dist

    diag(inter.dist) <- NA
    ##kinter=1
    if(!is.null(kinter)) {
        d1 = inter.dist*NA
        for(i in 1:nrow(inter.dist)) {
            jj <- Matrix::head(order(inter.dist[i,]),kinter)
            d1[i,jj] <- d1[jj,i] <- inter.dist[i,jj]
        }
        d1
        inter.dist = d1
    }
    inter.dist
    intra.dist

    intra.dist0 = intra.dist + 0.1 * mean(intra.dist,na.rm=TRUE)
    S = outer(intra.dist0, intra.dist0, FUN='+')**0.5
    inter.t = inter.dist / S
    inter.t
    score = mean(inter.t, na.rm=TRUE)
    return(score)
}


##pred.var="antigen.epitope";K=5;ntest=1000;x.var=c("cdr3.alpha","cdr3.beta")
knn.predict <- function(data, pred.var, x.var=NULL, K, samples=NULL, ntest=0.33 )
{
    
    
    qvars = names(data$Q)
    x.var0 = x.var
    if( is.null(x.var) || "<all others>" %in% x.var ) {
        x.var <- setdiff(qvars,pred.var)
    }
    x.var = intersect(x.var, qvars)
    if(is.null(samples)) samples = 1:nrow(data$Q[[1]])
    X = do.call(cbind, data$Q[x.var])[samples,]
    Y = as.character(data$tcr[samples,pred.var])
    names(Y) <- rownames(data$tcr)[samples]
    dim(X)

    if(class(X)=="matrix") {
        X <- Matrix::Matrix(X, sparse=TRUE)
    }

    ## random split into prediction and training samples
    cat("ntest=",ntest,"\n")
    if(ntest>0 && ntest<1) ntest = round(ntest*nrow(X))
    cat("ntest=",ntest,"\n")
    cat("nrowX=",nrow(X),"\n")
    j1 = sample(nrow(X), ntest) ## prediction
    j0 = setdiff(1:nrow(X),j1)
    length(j1)
    length(j0)

    ## distance calculations for test-set (superfast!)
    pos = dist.xy = NULL
    ##system.time( dist.xy <- 1 - qlcMatrix::corSparse( t(X[j1,]), t(X[j0,]) ) )
    if(sum(is.na(X))>0) {
        system.time( dist.xy <- 1 - cosine.similarity( t(X[j1,]), t(X[j0,]) ) )
    } else {
        system.time( dist.xy <- 1 - qlcMatrix::cosSparse( t(X[j1,]), t(X[j0,]) ) )
    }
    dim(dist.xy)

    ## Get best prediction from K-nearest neighbours
    knn <- apply(dist.xy, 1, function(x) Matrix::head(order(x),K))
    if(NCOL(knn)==1) {
        pred.yy  <- matrix(Y[j0][knn],ncol=1)
    } else {
        pred.yy  <- t(apply(knn, 2, function(i) Y[j0][i]))
    }
    colnames(pred.yy) <- paste0("nnb",1:ncol(pred.yy))
    dim(pred.yy)
    pred.max   <- apply(pred.yy,1,function(x) names(sort(-table(x)))[1])
    names(pred.max) <- names(Y)[j1]

    ## Accuracies
    test.y = Y[j1]
    acc=sens=specf=NA
    acc = mean( pred.max == test.y )
    acc

    res <- list( y.pred=pred.max, y=Y[j1],
                pred.var=pred.var, x.input=x.var0, x.var=x.var,
                train.idx=j0, test.idx=j1,
                accuracy=acc, sensitivity=sens, specificity=specf,
                pos=pos )
    return(res)
}

##align=TRUE
tagged.hamming <- function(aa,bb,align=TRUE)
{
    aligned.dist <- function(seq1, seq2) {
        seq1 <- gsub("[ ]","",seq1)
        seq2 <- gsub("[ ]","",seq2)
        seq.distance(c(seq1,seq2),verbose=0)$distance[1,2]
    }
    hamming.dist <- function(s,t) {
        s1=sapply(s,strsplit,split=" ")[[1]]
        t1=sapply(t,strsplit,split=" ")[[1]]
        n <- max(length(s1),length(t1))
        s1 = c(s1,rep(" ",n))[1:n]
        t1 = c(t1,rep(" ",n))[1:n]
        ##sum(s1!=t1) / n
        sum(s1!=t1)
    }
    tag.hamming0 <- function(a,b) {
        a1 = parse.tags(a)
        b1 = parse.tags(b)
        tags = intersect(names(a1),names(b1))
        if(length(tags)==0) return(NA)
        tags
        dist <- rep(NA,length(tags))
        names(dist) <- tags
        for(i in 1:length(tags)) {
            k = tags[i]
            is.seq = grepl("cdr|seq",k)
            if(is.seq && align) {
                dist[i] = aligned.dist(a1[k],b1[k])
            } else if(is.seq && !align) {
                dist[i] = hamming.dist(a1[k],b1[k])
            } else {
                dist[i] = 1*(a1[k] != b1[k])
            }
        }
        dist
    }
    a = aa[1]
    b = bb[1]
    aa.tags = unique(unlist(lapply(aa,function(a) names(parse.tags(a)))))
    bb.tags = unique(unlist(lapply(bb,function(b) names(parse.tags(b)))))
    all.tags <- sort(intersect(aa.tags, bb.tags))
    all.tags
    if(length(all.tags)==0 || is.null(all.tags)) return(NA)
    if(length(aa)==1 && length(bb)>1) aa <- rep(aa,length(bb))
    if(length(bb)==1 && length(aa)>1) bb <- rep(bb,length(aa))
    D = lapply(1:length(aa), function(i) tag.hamming0(aa[i], bb[i]))
    D <- t(sapply(D, function(e) e[match(all.tags,names(e))] ))
    D
    if(length(all.tags)==1) {
        D <- t(D)
        rownames(D) <- NULL
    }
    colnames(D) = all.tags
    D
}

hamming.match <- function(a, b, d=0) {
    any(rowSums(tagged.hamming(a,b)) < d)
}
hamming.min <- function(a, b) {
    min(rowSums(tagged.hamming(a,b)))
}

uscale <- function(x, symm=FALSE) {
    uscale.func <- function(x) (x-min(x))/(max(x)-min(x))
    if(NCOL(x)==1) {
        y = uscale.func(x)
    } else {
        y =  apply(x,2,uscale.func)
    }
    y[is.na(y)] <- NA
    if(symm) y <- (y-0.5)*2
    return(y)
}

fastcorSAVE <- function(x){
    sx <- scale(x)
    (t(sx) %*% sx) / (nrow(x)-1)
}

fastcor <- function(a) {
    a = Matrix::Matrix(a) ## automatically switch to sparse
    ##a <- as(a, "dgCMatrix")
    ##b <- as(b, "dgCMatrix")
    n <- nrow(a)
    muX <- colMeans(a)
    system.time( covmat <- (as.matrix(crossprod(a)) - n*tcrossprod(muX))/(n-1) )
    sdvec <- sqrt(diag(covmat))
    r <- covmat / tcrossprod(sdvec)
    r
}


mat2sp.DEPRECATED <- function(x) {
    
    idx = which(x!=0, arr.ind=TRUE)
    dn=list( rownames(x), colnames(x))
    spX <- Matrix::sparseMatrix( i=idx[,1], j=idx[,2], x=x[idx],
                                dims=c(nrow(x),ncol(x)), dimnames=dn )
    spX
}

sparse.corDEPRECATED <- function(a,b=NULL){
    ## see also corSparse in package 'qlcMatrix'
    
    if(is.null(b)) b=a
    a=Matrix(a, sparse=TRUE)
    b=Matrix(b, sparse=TRUE) ## make sure it is sparse coded
    n <- nrow(a)
    aMeans <- Matrix::colMeans(a,na.rm=TRUE)
    bMeans <- Matrix::colMeans(b,na.rm=TRUE)
    ##aSD <- apply(a,2,sd,na.rm=TRUE)
    ##bSD <- apply(b,2,sd,na.rm=TRUE)
    aSD <- ((Matrix::colSums(a**2,na.rm=TRUE) - n*aMeans**2)/(n-1))**0.5
    bSD <- ((Matrix::colSums(b**2,na.rm=TRUE) - n*bMeans**2)/(n-1))**0.5
    mm <- tcrossprod(aMeans, bMeans)
    ss <- tcrossprod(aSD, bSD)
    rr <- as.matrix(crossprod(a,b))
    covmat <- (rr - n*mm) / (n-1)
    cormat <- covmat / ss
    return(cormat)
}



## NEED RETHINK!! for non centereds, non-scaled values!!!
##fun="mean";m=1000,n=1000
mat.downsample <- function(mat, m, n=-1, FUN=c("mean","max"),
                           clust=FALSE )
{
    downsample.rows <- function(x, n=1000, FUN, clust=FALSE ) {
        if( NCOL(x) <= n  ) return(x)
        rowSymmMax <- function(a) {
            mx <- cbind( apply(a,1,max,na.rm=TRUE),
                        -apply(a,1,min,na.rm=TRUE))
            mx[is.infinite(mx)] <- 0 ## OK???
            apply(mx,1,max) * c(1,-1)[max.col(mx)]
        }
        if(clust) {
            x <- x[,nclust(x, distance="euclidean", link="complete")$order]
        }
        nn <- ncol(x)
        nd <- round(nn/n)
        jj <- seq(1,nn,nd)  ## start index of block
        col.mx <- colMeans(x,na.rm=TRUE)
        col.sd <- apply(x,2,sd,na.rm=TRUE)
        col.mx0 <- sapply(jj, function(i) mean(col.mx[i:min(nn,i+nd-1)],na.rm=TRUE))
        ##col.sd0 <- sapply(jj, function(i) max(col.sd[i:min(nn,i+nd-1)],na.rm=TRUE))
        dx <- scale(x, scale=FALSE)  ## just mean center
        ii <- lapply(jj, function(i) { i:min(nn,i+nd-1)} )
        sdx <- apply(dx,2,sd,na.rm=TRUE)
        ii <- lapply(ii, function(a) a[order(-sdx[a])] )  ## order on SD
        if(FUN=="max") {
            sx <- sapply(ii, function(i) rowSymmMax(dx[,i,drop=FALSE]))
        } else if(FUN=="mean") {
            sx <- sapply(ii, function(i) rowMeans(dx[,i,drop=FALSE],na.rm=TRUE))
        } else {
            stop("error. unknown function")
        }
        if( class(sx)=="numeric") sx <- matrix(sx,nrow=1)
        sx <- t(t(sx) + col.mx0) ## add back pooled col.mean
        dim(sx)
        ## set rownames
        aa <- sapply(ii, function(i) paste(colnames(x)[i],collapse=" "))
        colnames(sx) <- aa
        rownames(sx) <- rownames(x)
        sx
    }
    if(NCOL(mat)==0 && n>0) return(mat)
    ##fun="max"
    FUN <- FUN[1]
    if(NCOL(mat)==1 && class(mat)!="matrix") {
        dx <- matrix(mat,ncol=1)
        rownames(dx) <- names(mat)
        n <- -1
    } else {
        dx <- mat
    }
    if(nrow(dx)==0 && m>0) return(mat)
    if(n < ncol(dx) && n>0) dx <- downsample.rows(dx, n, FUN=FUN, clust=clust)
    if(m < nrow(dx) && m>0) dx <- t(downsample.rows(t(dx), m, FUN=FUN, clust=clust))
    dim(dx)
    if(NCOL(mat)==1) dx <- dx[,1]
    return(dx)
}

##===================================================================================
##===================================================================================
##===================================================================================
