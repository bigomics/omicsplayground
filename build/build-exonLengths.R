##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

martMMUSCULUS <- function() {


}

hg19GeneLengths <- function(symbols=NULL) {

    require(org.Hs.eg.db)
    
    if(is.null(symbols)) {
        symbols = as.character(unlist(as.list(org.Hs.egSYMBOL)))
    }
    exons.db = GenomicFeatures::exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')
    egs  = unlist(mget(symbols[ symbols %in% biomaRt::keys(org.Hs.egSYMBOL2EG) ],
                       envir=org.Hs.egSYMBOL2EG))
    exons.len = parallel::mclapply(egs[],function(eg) {
        exons = exons.db[[eg]]
        if(is.null(exons)) return(NA)
        exons = reduce(exons)
        sum( Biostrings::width(exons) )
    })
    exons.len <- unlist(exons.len)
    exons.len <- exons.len[!is.na(exons.len)]
    exons.len <- tapply(exons.len, names(exons.len), mean)
    return(exons.len)
}

mm10GeneLengths <- function(symbols=NULL) {
    ##biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
    require(org.Mm.eg.db)    
    if(is.null(symbols))
        symbols = as.character(unlist(as.list(org.Mm.egSYMBOL)))
    Matrix::head(symbols)
    exons.db = GenomicFeatures::exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by='gene')
    egs  = unlist(mget(symbols[ symbols %in% biomaRt::keys(org.Mm.egSYMBOL2EG) ],
                       envir=org.Mm.egSYMBOL2EG))
    exons.len = parallel::mclapply(egs[],function(eg) {
        exons = exons.db[[eg]]
        if(is.null(exons)) return(NA)
        exons = reduce(exons)
        sum( Biostrings::width(exons) )
    })
    exons.len <- unlist(exons.len)
    exons.len <- exons.len[which(!is.na(exons.len))]
    exons.len <- tapply(exons.len, names(exons.len), mean)
    return(exons.len)
}


if(0) {
    exons.len <- hg19GeneLengths()
    Matrix::head(exons.len)
    write.csv( data.frame(symbol=names(exons.len), exon.length = exons.len),
          file="../files/hg19GeneExonLengths.csv", row.names=FALSE, quote=FALSE)
    
    mm10.len <- mm10GeneLengths()
    Matrix::head(mm10.len,100)
    Matrix::tail(mm10.len,100)
    write.csv( data.frame(symbol=names(mm10.len), exon.length = mm10.len),
              file="../files/mm10GeneExonLengths.csv", row.names=FALSE, quote=FALSE)
}

