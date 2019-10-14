## Build all examples data sets
##
##


## all scripts
all.scripts <- dir(".", pattern="pgx-.*R$")

## These are the default example scripts (uncomment if you really want to do all)
all.scripts <- grep("GSE10846|GSE114716|GSE22886|GSE28492|GSE32591|GSE53784|GSE88808|GSE72056|GSE92332|GSE98638|geiger2016|rieckmann2017|tcga-brca|tcga-prad", all.scripts, ignore.case=TRUE, value=TRUE)
all.scripts

script=all.scripts[1]
for(script in all.scripts) {

    ## skip if already done
    pgx.file <- gsub("pgx-|[.]R$","",script)
    if(any(grepl(pgx.file, dir("../data")))) next

    ## run script
    cat(">>>>>>>>>>>>>>>>>> processing",script,"<<<<<<<<<<<<<<<\n")
    source(script, local=FALSE)

    ## clean up
    rm(list=setdiff(ls(),c("script","all.scripts")))
}

## scan and update datasets info
source("update-datasets-info.R")

