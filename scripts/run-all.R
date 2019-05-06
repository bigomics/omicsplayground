


all.scripts <- dir(".", pattern="pgx-gei.*R$")
all.scripts <- dir(".", pattern="pgx-.*R$")
all.scripts <- grep("dlbcl|ipi|melano|rieck",all.scripts,value=TRUE)
all.scripts

script=all.scripts[1]

for(script in all.scripts) {

    cat(">>>>>>>>>>>>>>>>>>>>>>>",script,"<<<<<<<<<<<<<<<<<<<<<<<\n")
    source(script, local=FALSE)
    
    ## clean up
    rm(list=setdiff(ls(),c("script","all.scripts")))    
}
