##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Build all examples data sets
##
##

## all scripts
all.scripts <- dir("../scripts", pattern="pgx-.*R$")

## These are some default example scripts (uncomment if you really want to do all)
all.scripts <- c("pgx-geiger2016-arginine.R","pgx-GSE72056-scmelanoma.R")
all.scripts
script=all.scripts[1]

for(script in all.scripts) {

    ## skip if already done
    pgx.file <- gsub("pgx-|[.]R$","",script)
    ##if(any(grepl(pgx.file, dir("../data")))) next

    ## run script
    cat(">>>>>>>>>>>>>>>>>> processing",script,"<<<<<<<<<<<<<<<\n")
    script1 <- paste0("../scripts/",script)
    try.err <- try(source(script1, local=FALSE))
    if(class(try.err)=="try-error") {
        cat("WARNING:: Error in source",script,"\n")
    }
    
    ## clean up
    rm(list=setdiff(ls(),c("script","all.scripts")))
}

## scan and update datasets info
## source("update-datasets-info.R")

