source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")
source("../R/xcr-graph.r")
source("../R/pgx-graph.R")
library(Rtsne.multicore)
library(Rtsne)
library(qlcMatrix)
library(corpora)
##devtools::install_github("bwlewis/rthreejs")

## all public datasets
PGX.FILES <- dir("../pgx", patter="GSE|argin|rieck|schmiedel")
PGX.FILES <- grep("-LT.pgx$",PGX.FILES,value=TRUE)
pub.id <- sub("-.*","",PGX.FILES)
PGX.FILES <- PGX.FILES[!duplicated(pub.id)]
PGX.FILES

allFC <- list()
pgx=PGX.FILES[2]
for(pgx in PGX.FILES) {
    load(file.path("../pgx",pgx),verbose=1)
    rownames(ngs$X) <- toupper(sub(".*:","",rownames(ngs$X)))
    meta.fx <- sapply(ngs$gx.meta$meta,function(x) x$meta.fx)
    rownames(meta.fx) <- toupper(rownames(ngs$gx.meta$meta[[1]]))
    allFC[[pgx]] <- meta.fx
}

## find common genes
all.gg <- toupper(as.character(unlist(sapply(allFC, rownames))))
gg.tbl <- table(all.gg)
table(gg.tbl)
length(unique(all.gg))

gg <- head(names(sort(-gg.tbl)),8000)
length(gg)
allFC <- lapply(allFC, function(x) {
    x <- x[match(gg,rownames(x)),]
    rownames(x) <- gg
    return(x)
})
str(allFC)

names(allFC) <- sub("-.*","",names(allFC))
for(i in 1:length(allFC)) {
    colnames(allFC[[i]]) <- paste0("[",names(allFC)[i],"] ",colnames(allFC[[i]]))
}

allFC <- do.call(cbind, allFC)
dim(allFC)
summary(colMeans(is.na(allFC)))

save(allFC, file="../files/allFoldChanges-pub-8k.rda")
##load(file="../files/allFoldChanges.rda", verbose=1)


##================================================================================
##================================================================================
##================================================================================
calcAllPairsMA <- function(X, y) {
    pairs <- t(combn(unique(y),2))
    M <- c()
    A <- c()
    Q <- c()
    i=1
    for(i in 1:nrow(pairs)) {
        j0 <- which( y == pairs[i,1])
        j1 <- which( y == pairs[i,2])

        ## compute mean, diff, q-value
        f1 <- Matrix::rowMeans(X[,j1,drop=FALSE]) - Matrix::rowMeans(X[,j0,drop=FALSE])
        m1 <- Matrix::rowMeans(X[,c(j0,j1),drop=FALSE])

        ## add reverse
        f2 <- cbind(f1, -f1)
        m2 <- cbind(m1, m1)
        colnames(f2) <- c( paste(pairs[i,2:1],collapse="_vs_"),
                          paste(pairs[i,1:2],collapse="_vs_") )
        colnames(m2) <- colnames(f2)

        M <- cbind(M, f2)
        A <- cbind(A, m2)
    }
    return(list(M=M, A=A, Q=Q))
}


if(0) {

    allX <- list()
    allY <- list()
    pgx=PGX.FILES[2]
    for(pgx in PGX.FILES) {
        load(file.path("../pgx",pgx),verbose=1)
        rownames(ngs$X) <- toupper(sub(".*:","",rownames(ngs$X)))
        allX[[pgx]] <- ngs$X
        ##allY[[pgx]] <- gsub("[-_.]","",ngs$samples$group)
        allY[[pgx]] <- ngs$samples$group
    }
    names(allX)
    
    ## find common genes
    ##gg <- Reduce(intersect, lapply(allX, rownames))
    gg.tbl <- table(unlist(sapply(allX, rownames)))
    ##gg <- names(which( gg.tbl >= length(allX)))
    table(gg.tbl)
    gg <- head(names(sort(-gg.tbl)),2000)
    length(gg)
    allX <- lapply(allX, function(x) {
        x <- x[match(gg,rownames(x)),]
        rownames(x) <- gg
        return(x)
    })
    
    
## compute all pairs MA
    allM <- c()
    allA <- c()
    i=1
    for(i in 1:length(allX)) {
        X=allX[[i]]
        y=allY[[i]]
        res <- calcAllPairsMA(X, as.character(y))
        bx <- strsplit(PGX.FILES[i], split="[-_]")[[1]][1]
        colnames(res$A) <- paste0("[",bx,"]",colnames(res$A))
        colnames(res$M) <- paste0("[",bx,"]",colnames(res$M))
        allM <- cbind(allM, res$M)
        allA <- cbind(allA, res$A)
    }
    dim(allM)
    dim(allA)
    head(colnames(allA))
    head(colnames(allM))
    
    save(allM, allA, file="../files/allMA-pub.rda")
    load(file="../files/allMA.rda", verbose=1)
}
