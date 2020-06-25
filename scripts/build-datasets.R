## Build all examples data sets
##
##


## all scripts
all.scripts <- dir(".", pattern="pgx-.*R$")

## These are some default example scripts (uncomment if you really want to do all)
all.scripts <- c("pgx-geiger2016-arginine.R","pgx-GSE72056-scmelanoma.R",
                 "pgx-GSE22886-immune.R","pgx-tcga-brca.R","pgx-GSE102908-ibet.R")
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
## source("update-datasets-info.R")

