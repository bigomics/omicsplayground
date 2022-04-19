##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

require(org.Hs.eg.db)
source("../R/pgx-functions.R")

cat("************************************************************************\n")
cat("*********************** BUILD ANNOTATIONS START ************************\n")
cat("************************************************************************\n")

cat("===================== BUILDING GENE ANNOTATIONS ======================\n")

eg="1"
getGeneInfo <- function(eg) {    
    a1 <- getMyGeneInfo(eg)
    a2 <- getHSGeneInfo(eg, as.link=FALSE)
    aa <- c(a1, a2[setdiff(names(a2),names(a1))])
    aa <- sapply(aa, paste, collapse=" /// ")
    aa['entrez'] <- eg
    kk <- c("entrez","symbol","name","alias","map_location","summary",
      "OMIM","KEGG","GO")
    aa <- aa[match(kk,names(aa))]
    names(aa) <- kk
    aa
}

all.eg <- keys(org.Hs.eg.db)
length(all.eg)
info = lapply( all.eg[1:2], getGeneInfo)
info = do.call(rbind,info)
dim(info)

cat("writing to gene-info.csv")
write.csv(info, file="../lib/gene-info.csv", row.names=FALSE)

##--------------------------------------------------------
## Human to mouse translation
##--------------------------------------------------------

##install.packages("homologene")
require("homologene")

hs.symbols <- sort(unique(unlist(as.list(org.Hs.egSYMBOL))))
length(hs.symbols)
mm <- homologene::human2mouse(hs.symbols)
head(mm)




cat("************************************************************************\n")
cat("************************ BUILD ANNOTATIONS DONE ************************\n")
cat("************************************************************************\n")
