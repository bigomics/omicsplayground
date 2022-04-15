##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##===================================================================
##======= Script to build PGX object from GEO dataset ===============
##===================================================================
##
##

OPG="~/Playground/omicsplayground/"
RDIR = file.path(OPG,"R")
FILES = file.path(OPG,"lib")
FILESX = file.path(OPG,"libx")
PGX.DIR = file.path(OPG,"data")
source(file.path(RDIR,"pgx-include.R"))
source(file.path(RDIR,"pgx-getgeo.R"))

archs.h5 = file.path(FILESX,"human_matrix.h5")

geo.ids = c("GSE141499","GSE80629","GSE80631","GSE2432","GSE106391")
##geo.ids = c("GSE2432","GSE106391")
id = "GSE106391"
id = "GSE2432"
id
for(id in geo.ids) {

    ##------------------------------------------------------------
    ## Get GEO data
    ##------------------------------------------------------------
    
    ## load series and platform data from GEO
    geo <- pgx.getGEOseries(id=id, archs.h5=archs.h5)
    names(geo)
    head(geo$counts)
    head(geo$samples)    
    head(geo$meta)
    geo$contrasts
    geo$samples$group
    
    ##------------------------------------------------------------
    ## Set PGX information
    ##------------------------------------------------------------

    counts = geo$counts
    samples = geo$samples
    contrasts = geo$contrasts
    
    pgx = pgx.createPGX(
        counts = counts,
        samples = samples,
        contrasts = contrasts
    )


    ## Add cluster contrasts
    pgx <- pgx.clusterSamples(pgx, perplexity=30)
    Y = pgx$samples[,"cluster",drop=FALSE]
    ct <- makeDirectContrasts(Y, ref="others")
    pgx$contrasts
    if(ncol(pgx$contrasts)==0) {
        pgx$contrasts <- ct$exp.matrix
    } else {
        pgx$contrasts <- cbind(pgx$contrasts, ct$exp.matrix)
    }
    
    gx.methods    = c("trend.limma","edger.qlf","deseq2.wald")
    gset.methods  = c("fisher","gsva","fgsea")
    extra.methods = c("deconv","wordcloud","connectivity")
    extra.methods  = c("meta.go","infer","drugs","deconv","wordcloud","connectivity")
    
    pgx <- pgx.computePGX(
        pgx,
        max.genes = 20000,
        max.genesets = 5000, 
        gx.methods = gx.methods,
        gset.methods = gset.methods,
        extra.methods = extra.methods,
        use.design = TRUE,      ## no.design+prune are combined 
        prune.samples = FALSE,  ##
        do.cluster = TRUE,                
        progress = NULL,
        lib.dir = FILES 
    )
    
    ##extra.methods <- c("connectivity")
    ##pgx <- compute.extra(pgx, extra.methods, lib.dir=FILES) 

    ##-------------------------------------------------------------------
    ## save PGX object
    ##-------------------------------------------------------------------
    
    rda.file = paste0("../data/",id,".pgx")
    ngs = pgx   ## still old...
    ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
    ngs$datatype = geo$info$type

    is.mouse <- (mean(grepl("[a-z]",rownames(pgx$counts))) > 0.5)
    ngs$organism = ifelse(is.mouse, "mouse", "human")
    ngs$description = paste0(geo$info$title,". ",geo$info$summary)
    
    rda.file
    ngs.save(ngs, file=rda.file)

}
