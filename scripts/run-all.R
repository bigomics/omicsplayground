if(0) {
    system("cp /home/kwee/Projects/R/gx/gx-limma.r ../R/")
    system("cp /home/kwee/Projects/R/gx/gx-heatmap.r ../R/")
    system("cp /home/kwee/Projects/R/gx/gx-volcano.r ../R/")
    system("cp /home/kwee/Projects/R/gx/gx-util.r ../R/")
    system("cp /home/kwee/Projects/R/gset/gset-fisher.r ../R/")
    system("cp /home/kwee/Projects/R/gset/gset-gsea.r ../R/")
    system("cp /home/kwee/Projects/R/gset/gset-meta.r ../R/")
    system("cp /home/kwee/Projects/R/ngs/ngs-salmon.r ../R/")
    system("cp /home/kwee/Projects/R/ngs/ngs-cook.r ../R/")
    system("cp /home/kwee/Projects/R/ngs/ngs-fit.r ../R/")
    system("cp /home/kwee/Projects/R/pgx/pgx-*.R ../R/")
    ##system("cp /home/kwee/Projects/R/pgx/module-*.Rmd ../R/")
}


rm(list=ls())
all.scripts <- dir(".", pattern="pgx-gei.*R$")
all.scripts <- dir(".", pattern="pgx-.*R$")
all.scripts

FAST=TRUE
FAST=FALSE
SMALL=4000
EXT="4x"
script=all.scripts[13]
run.param = c("run.param",ls())

for(script in all.scripts[11]) {
    cat(">>> running",script,"...\n")
    source(script, local=FALSE)
    rm(list=setdiff(ls(),run.param))
}
