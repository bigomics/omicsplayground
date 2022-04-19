##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2021 BigOmics Analytics Sagl. All rights reserved.
##

##===================================================================
##================ Script to build PGX object =======================
##===================================================================
##
##
##
##


RDIR = "../R"
FILES = "../lib"
PGX.DIR = "../data"
source("../R/pgx-include.R")
source("../R/pgx-tcga.R")
##source("options.R")

library('rhdf5')
h5ls("../libx/tcga_matrix.h5")
X = h5read("../libx/tcga_matrix.h5",'data/expression')
genes = h5read("../libx/tcga_matrix.h5",'meta/genes')
case.id = h5read("../libx/tcga_matrix.h5",'meta/gdc_cases.submitter_id')
colnames(X) = case.id
rownames(X) = genes

dim(X)
gx.hist(log2(pmax(X[,1:40])))
summary(colSums(X))
gx.hist( log2(pmax(X[,1:40]/1000,1)))


surv.file <- file.path(FILES, "rtcga-survival.csv")
surv <- read.csv(surv.file, row.names=1)
head(surv)
surv$months <- round(surv$times/365*12, 2)
surv$status <- surv$patient.vital_status

kk = intersect(rownames(surv),colnames(X))
counts = X[,kk]
samples = surv[kk,]

##------------------------------------------------------------
## Make survival contrasts
##------------------------------------------------------------

contrasts <- samples[,0]
cc=surv$cancer_type[1]
cc.types <- sort(unique(surv$cancer_type))
for(cc in cc.types) {
    sel <- which(surv$cancer_type == cc & surv$status == 1)
    med.survtime <- median(surv$months[sel], na.rm=TRUE)
    sel1 <- sel[surv$months[sel] < med.survtime]
    sel0 <- sel[surv$months[sel] >= med.survtime]
    sel0 <- rownames(surv)[sel0]
    sel1 <- rownames(surv)[sel1]
    ct <- -1*(rownames(samples) %in% sel0) + 1*(rownames(samples) %in% sel1)
    contrasts <- cbind(contrasts, ct)        
}
colnames(contrasts) <- cc.types
dim(contrasts)
head(contrasts)
colSums(contrasts!=0)

sel <- (colSums(contrasts<0) >= 3 & colSums(contrasts>0) >= 3)
table(sel)
contrasts <- contrasts[,sel]

##------------------------------------------------------------
## Set PGX information
##------------------------------------------------------------

pgx = pgx.createPGX(
    counts = counts,
    samples = samples,
    contrasts = contrasts
)

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

rda.file = "../data/tcga-survival.pgx"
ngs = pgx   ## still old...
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = geo$info$type

is.mouse <- (mean(grepl("[a-z]",rownames(pgx$counts))) > 0.5)
ngs$organism = ifelse(is.mouse, "mouse", "human")
ngs$description = paste0(geo$info$title,". ",geo$info$summary)

rda.file
ngs.save(ngs, file=rda.file)
