## library(limma)
## source("../R/gx-heatmap.r")
## source("../R/gx-limma.r")
## source("../R/gx-util.r")
## source("../R/gset-fisher.r")
## source("../R/gset-gsea.r")
## source("../R/gset-meta.r")

SAVE.PARAMS <- ls()

##-----------------------------------------------------------
## Load gene sets
##-----------------------------------------------------------

pp <- rownames(ngs$counts)
is.mouse = (mean(grepl("[a-z]",sub(".*:","",pp))) > 0.8)
is.mouse
if(is.mouse) {
    cat("Loading mouse gene sets...\n")
    load(file=file.path(FILES,"gmt-all-mouse.rda"))
} else {
    cat("Loading human gene sets...\n")
    load(file=file.path(FILES,"gmt-all.rda"))
}
table(sub(":.*","",names(gmt.all)))
summary(sapply(gmt.all,length))

##-----------------------------------------------------------
## Filter gene sets
##-----------------------------------------------------------
cat("Filtering gene sets...\n")

require(Matrix)
require(org.Hs.eg.db)
##GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
##genes = head(as.character(unlist(as.list(org.Hs.egSYMBOL))),1000)
genes = unique(as.character(ngs$genes$gene_name))
gmt.all = mclapply(gmt.all, function(gs) intersect(gs, genes))
table(sub(":.*","",names(gmt.all)))

## filter gene sets on size
gmt.size = sapply(gmt.all,length)
summary(gmt.size)
gmt.all = gmt.all[which(gmt.size >= 15 & gmt.size <= 1e3)]
gmt.size = sapply(gmt.all,length)
gmt.all = gmt.all[order(-gmt.size)]
summary(gmt.size)
gmt.all = gmt.all[!duplicated(names(gmt.all))]
length(gmt.all)
table(sub(":.*","",names(gmt.all)))

##-----------------------------------------------------------
## create the full GENE matrix (always collapsed by gene)
##-----------------------------------------------------------

genes = as.character(ngs$genes$gene_name)
X <- apply(ngs$counts, 2, function(x) tapply(x, genes, sum))
dim(X)
X <- edgeR::cpm(X, log=TRUE )
X <- X[which(apply(X,1,sd)>0),,drop=FALSE]
X <- X[!(rownames(X) %in% c(NA,""," ") ),,drop=FALSE]
X <- limma::normalizeQuantiles(X)
dim(X)

##-----------------------------------------------------------
## create the GENESETxGENE matrix
##-----------------------------------------------------------
cat("Building gene set matrix...\n")

##GMT = sapply( gmt.all, function(s) 1*(rownames(X) %in% s))
##GMT = Matrix(GMT, sparse=TRUE)
GMT <- gmt2mat.nocheck(gmt.all[], bg=rownames(X))  ## in gset-gsea.r
dim(GMT)
dim(X)
table(rownames(X) %in% rownames(GMT))
table(names(gmt.all) %in% colnames(GMT))
GMT[1:4,1:4]
GMT <- GMT[rownames(X),names(gmt.all)]

summary(Matrix::colSums(GMT))
dim(GMT)

##-----------------------------------------------------------
## Prioritize gene sets by fast rank-correlation
##-----------------------------------------------------------

if(!exists("MAX.GENES")) MAX.GENES <- 20000
if(MAX.GENES < 0) MAX.GENES <- 20000
MAX.GENES

if(MAX.GENES > 0) {
    require(limma)
    ## Reduce gene sets by selecting top varying genesets. We use the
    ## very fast sparse rank-correlation for approximate single sample
    ## geneset activation.
    gsetX = qlcMatrix::corSparse( GMT[,], apply(X[,],2,rank) )
    gsetX = limma::normalizeQuantiles(gsetX) ##???
    grp <- ngs$samples$group
    gsetX.bygroup <- t(apply(gsetX,1,function(x) tapply(x,grp,mean)))
    sdx <- apply(gsetX.bygroup,1,sd)
    names(sdx) <- colnames(GMT)
    jj = head(order(-sdx), MAX.GENES) 
    must.include <- "hallmark|kegg|^go|^celltype"
    jj = unique( c(jj, grep(must.include,colnames(GMT),ignore.case=TRUE)))
    jj = jj[order(colnames(GMT)[jj])]
    length(jj)
    GMT = GMT[,jj]
    gmt.all = gmt.all[colnames(GMT)]
}
dim(GMT)

sum(duplicated(rownames(GMT)))
dim(X)
dim(GMT)
ngs$GMT <- GMT
ngs$gmt.all <- gmt.all

##-----------------------------------------------------------
## get design and contrast matrix
##-----------------------------------------------------------
design = ngs$model.parameters$design
contr.matrix = ngs$model.parameters$contr.matrix
## exp.matrix = (design %*% contr.matrix)

test.methods=c("fisher","ssgsea","gsva", "spearman", "camera", "fry",
               "fgsea","gsea.permPH","gsea.permGS","gseaPR")
test.methods = c("fisher","fgsea")
test.methods = c("fisher","gsva","ssgsea","spearman","camera","fry","fgsea") ## no GSEA
if(!is.null(USER.GENESETTEST.METHODS)) test.methods = USER.GENESETTEST.METHODS

##-----------------------------------------------------------
## Run methods
##-----------------------------------------------------------
cat(">>> Testing gene sets with methods:",test.methods,"\n")
test.methods

Y <- ngs$samples
##gmt=ngs$gmt.all[1:100]
gmt=ngs$gmt.all[]
gc()

gset.meta = gset.fitContrastsWithAllMethods(
    gmt = ngs$gmt.all, X = X, Y = Y, design=design, ## genes=GENES,
    contr.matrix=contr.matrix, methods=test.methods,
    mc.threads=1, mc.cores=NULL, batch.correct=TRUE )

print(gset.meta$timings)

rownames(gset.meta$timings) <- paste0("[testgenesets]",rownames(gset.meta$timings))
ngs$timings <- rbind(ngs$timings, gset.meta$timings)
ngs$gset.meta <- gset.meta
ngs$gsetX = ngs$gset.meta$matrices[["meta"]]  ## META??!

##-----------------------------------------------------------------------
##------------------------ clean up -------------------------------------
##-----------------------------------------------------------------------

## remove large outputs... (uncomment if needed)
ngs$gset.meta$outputs <- NULL

remove(X)
remove(Y)
rm(list=setdiff(ls(),SAVE.PARAMS))

