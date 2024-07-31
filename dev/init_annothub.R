##
## Initialize (prepopulate) AnnotationHub cache with main species
##
##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

if(1) {
  if(!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
  if(!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
  if(!require("org.Rn.eg.db")) BiocManager::install("org.Rn.eg.db")
}

if(0) {
  ##dir("~/.cache/R/AnnotationHub")
  ah <- AnnotationHub::AnnotationHub()  
  species <- c(
    "Homo sapiens",
    "Mus musculus",
    "Rattus norvegicus"
    ##  "Saccharomyces cerevisiae",
    ##  "Drosophila melanogaster",
    ##  "Caenorhabditis elegans"
  )
  s <- species[1]
  for(s in species) {
    cat("Retrieving OrgDb for",s,"...\n")
    res <- AnnotationHub::query(ah, pattern = c(s, "OrgDb"))  
    suppressMessages(
      db <- try(res[[length(res)]])
    )
  }
}

