

all.scripts <- dir(".", pattern="pgx-.*R$")
##all.scripts <- grep("dlbcl|ipi|melano|rieck",all.scripts,value=TRUE)
all.scripts

script=all.scripts[1]
for(script in all.scripts) {

    ## skip if already done
    pgx.file <- gsub("pgx-|[.]R$","",script)
    if(any(grepl(pgx.file, dir("../pgx")))) next

    ## run script
    cat(">>>>>>>>>>>>>>>>>> processing",script,"<<<<<<<<<<<<<<<\n")
    source(script, local=FALSE)

    ## clean up
    rm(list=setdiff(ls(),c("script","all.scripts")))
}
