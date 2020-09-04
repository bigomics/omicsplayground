##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

martMMUSCULUS <- function() {


}

hg19GeneLengths <- function(symbols=NULL) {
    require(parallel)
    require(GenomicFeatures)
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    require(org.Hs.eg.db)
    if(is.null(symbols))
        symbols = as.character(unlist(as.list(org.Hs.egSYMBOL)))
    exons.db = exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')
    egs  = unlist(mget(symbols[ symbols %in% keys(org.Hs.egSYMBOL2EG) ],org.Hs.egSYMBOL2EG))
    exons.len = mclapply(egs[],function(eg) {
        exons = exons.db[[eg]]
        if(is.null(exons)) return(NA)
        exons = reduce(exons)
        sum( width(exons) )
    })
    exons.len <- unlist(exons.len)
    exons.len <- exons.len[!is.na(exons.len)]
    exons.len <- tapply(exons.len, names(exons.len), mean)
    return(exons.len)
}

mm10GeneLengths <- function(symbols=NULL) {
    ##biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
    require(parallel)
    require(GenomicFeatures)
    require(TxDb.Mmusculus.UCSC.mm10.knownGene)
    require(org.Mm.eg.db)
    if(is.null(symbols))
        symbols = as.character(unlist(as.list(org.Mm.egSYMBOL)))
    head(symbols)
    exons.db = exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by='gene')
    egs  = unlist(mget(symbols[ symbols %in% keys(org.Mm.egSYMBOL2EG) ],org.Mm.egSYMBOL2EG))
    exons.len = mclapply(egs[],function(eg) {
        exons = exons.db[[eg]]
        if(is.null(exons)) return(NA)
        exons = reduce(exons)
        sum( width(exons) )
    })
    exons.len <- unlist(exons.len)
    exons.len <- exons.len[which(!is.na(exons.len))]
    exons.len <- tapply(exons.len, names(exons.len), mean)
    return(exons.len)
}


exons.len <- hg19GeneLengths()
head(exons.len)
write.csv( data.frame(symbol=names(exons.len), exon.length = exons.len),
          file="../files/hg19GeneExonLengths.csv", row.names=FALSE, quote=FALSE)

mm10.len <- mm10GeneLengths()
head(mm10.len,100)
tail(mm10.len,100)
write.csv( data.frame(symbol=names(mm10.len), exon.length = mm10.len),
          file="../files/mm10GeneExonLengths.csv", row.names=FALSE, quote=FALSE)


