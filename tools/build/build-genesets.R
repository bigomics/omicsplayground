##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

source("../R/gx-util.r")
source("../R/gset-gsea.r")
source("../R/gset-fisher.r")

##----------------------------------------------------------------------
##--------------- cell type signatures from ImSig ----------------------
##----------------------------------------------------------------------
S <- read.csv("../opt/imsig-signatures-genes.csv", skip=2)
imsig.gmt <- tapply( S$Gene.Symbol, S$ImSig.cell.type, list)
imsig.gmt <- lapply(imsig.gmt, as.character)
write.gmt(imsig.gmt, file="../files/gmt/celltype_imsig.gmt")
write.gmt(imsig.gmt, file="../files/gmt0/celltype_imsig.gmt")


##----------------------------------------------------------------------
##----------- cell type signatures from xCELL (collapsed) --------------
##----------------------------------------------------------------------

xcell.data <- read.gmt("../opt/xCell_celltype_signatures.txt")
cell.type <- sub("_.*","",names(xcell.data))
table(cell.type)
xcell.gmt <- tapply( xcell.data, cell.type, function(s) sort(unique(unlist(s))))
write.gmt(xcell.gmt, file="../files/gmt/celltype_xcell.gmt")
write.gmt(xcell.gmt, file="../files/gmt0/celltype_xcell.gmt")


##----------------------------------------------------------------------
##----------- build the GMT-all object (all gene sets) -----------------
##----------------------------------------------------------------------

##tt = read.gmt("../../Data/Enrichr/Tissue_Jensen.txt")
##head(sapply(tt,length),100)

require(parallel)
gmt.files2 = dir("../files/gmt0", pattern=".gmt$|.txt$", full.names=TRUE)
gmt.all = mclapply(gmt.files2, read.gmt)
names(gmt.all) = gmt.files2
names(gmt.all) = gsub(".*/|.txt$|.gmt$", "", names(gmt.all))
gmt.db = gsub("[_.-].*|.txt$|.gmt$", "", names(gmt.all))
gmt.db = toupper(gmt.db)

table(gmt.db)
for(i in 1:length(gmt.all)) {
    names(gmt.all[[i]]) <- sub("\\(GO:","(GO_",names(gmt.all[[i]]))
    names(gmt.all[[i]]) <- sub(":","",names(gmt.all[[i]]))
    ## names(gmt.all[[i]]) <- tolower(names(gmt.all[[i]]))
    names(gmt.all[[i]]) <- paste0(toupper(gmt.db[i]),":",names(gmt.all[[i]]))
}
j0 = grep("_up", names(gmt.all))
j1 = grep("_down", names(gmt.all))
for(i in j0) {
    names(gmt.all[[i]]) <- paste0(names(gmt.all[[i]])," (up)")
}
for(i in j1) {
    names(gmt.all[[i]]) <- paste0(names(gmt.all[[i]])," (down)")
}
names(gmt.all) <- NULL
gmt.all <- unlist(gmt.all,recursive=FALSE, use.names=TRUE)
length(gmt.all)

## get rid of trailing numeric values
gmt.all <-  mclapply(gmt.all, function(x) gsub("[,].*","",x), mc.cores=4)

## order by length and take out duplicated sets (only by name)
gmt.all <- gmt.all[order(-sapply(gmt.all,length))]
gmt.all <- gmt.all[!duplicated(names(gmt.all))]


## save
gmt.all <- gmt.all[order(names(gmt.all))]
table(sub(":.*","",names(gmt.all)))
save(gmt.all, file="../files/gmt-all.rda")
saveRDS(gmt.all, file="../files/gmt-all.rds")


## NEED RETHINK!!!! this should be improved using BioMart...
cat("Converting gsets to mouse ID...\n")
library(org.Mm.eg.db)
mouse.genes = as.character(unlist(as.list(org.Mm.egSYMBOL)))
names(mouse.genes) = toupper(mouse.genes)
gmt.all <- mclapply(gmt.all[], function(s) setdiff(as.character(mouse.genes[s]),NA), mc.cores=4)
save(gmt.all, file="../files/gmt-all-mouse.rda")
saveRDS(gmt.all, file="../files/gmt-all-mouse.rds")
remove(gmt.all)




##======================================================================
##======================================================================
##======================================================================

