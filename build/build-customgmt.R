##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")
source("../R/gset-fisher.r")
source("../R/gset-gsea.r")
source("../R/gset-meta.r")

## HUGO gene families
G = read.csv("hgnc-genefamilies.txt",sep="\t")
fam = tapply( G$Approved.Symbol, G$Gene.family.description, list )
##names(fam) = toupper(gsub("[ -/]","_",names(fam)))
##names(fam) = gsub("_$","",names(fam))
##names(fam) = gsub("[,.]","",names(fam))
fam.size = sapply(fam, length)
summary(fam.size)
fam = fam[which(fam.size >= 50 & fam.size <=1000) ]
fam.size = sapply(fam, length)
summary(fam.size)
length(fam)
write.gmt(fam, file="gmt/hgnc_genefamilies.gmt", source="www.genenames.org")


## MSIGDB gene families
G = read.csv("msigdb-genefamilies.csv")
fam = apply( G, 2, function(x) setdiff(x,""))
##names(fam) = toupper(gsub("[ -/]","_",names(fam)))
##names(fam) = gsub("_$","",names(fam))
##names(fam) = gsub("[,.]","",names(fam))
fam.size = sapply(fam, length)
summary(fam.size)

fam = fam[which(fam.size >= 50 & fam.size <=9999) ]
fam.size = sapply(fam, length)
summary(fam.size)
length(fam)
write.gmt(fam, file="gmt/msigdb_genefamilies.gmt", source="MSigDB")
write.gmt(fam, file="gmt0/msigdb_genefamilies.gmt", source="MSigDB")

## ABC/GCB
gmt <- read.gmt("~/IRB/Data/gmt/IOR-custom-genesets_1704.gmt",dir=NULL)
abcgcb.genes <- unique(unlist(gmt[grep("ABC_GCB",names(gmt))]))
write(abcgcb.genes, file="abcgcb-genes.txt")

## DSigDB
dsig <- read.csv("~/IRB/Data/gmt/src/DSigDBv1.0_all.txt",header=FALSE)
dsig <- as.character(dsig[,1])
head(dsig)
jj = c(grep(" : ",dsig),length(dsig)+1)
head(jj)
gmt <- list()
i=1
for(i in 1:(length(jj)-1)) {
    aa = dsig[ jj[i]:(jj[i+1]-1)]
    gmt[[ aa[1] ]] <- aa[2:length(aa)]
}
length(gmt)
gmt1 = gmt[grep("Compound\\(D1\\)",names(gmt))]
length(gmt1)
gmt2 = gmt[grep("\\(D2",names(gmt))]
length(gmt2)
gmt3 = gmt[grep("Compound\\(D3\\)",names(gmt))]
length(gmt3)
names(gmt1) = sub("Compound.*: ","",names(gmt1))
names(gmt2) = sub("Compound.*: ","",names(gmt2))
names(gmt2) = sub("D2 ","",names(gmt2))
names(gmt3) = sub("Compound.*: ","",names(gmt3))

write.gmt(gmt1, file="gmt/dsigdb_FDAapprovedD1.gmt", source="DSigDB")
write.gmt(gmt2, file="gmt/dsigdb_KinaseInhibitorsD2.gmt", source="DSigDB")
write.gmt(gmt3, file="gmt/dsigdb_signaturesD3.gmt", source="DSigDB")

