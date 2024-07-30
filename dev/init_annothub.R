##
## Initialize (prepopulate) AnnotationHub cache with main species
##
##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

##dir("~/.cache/R/AnnotationHub")
ah <- AnnotationHub::AnnotationHub()

species <- c(
  "Homo sapiens",
  "Mus musculus",
  "Rattus norvegicus"
#  "Saccharomyces cerevisiae",
#  "Drosophila melanogaster",
#  "Caenorhabditis elegans"
)

s <- species[1]
for(s in species) {
  cat("Retrieving OrgDb for",s,"...\n")
  res <- AnnotationHub::query(ah, pattern = c(s, "OrgDb"))  
  suppressMessages(
    db <- try(res[[length(res)]])
  )
}

