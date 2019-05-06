source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")
source("../R/pgx-functions.R")

FILES = "../files/"
RDIR = "../R/"
PGX.FILES <- dir("../pgx", pattern=".pgx$",full.names=TRUE)
PGX.FILES = grep("-LT.pgx", PGX.FILES, value=TRUE, invert=TRUE)  ## only heavy
##PGX.FILES = grep("GSE|rieck|schmiedel", PGX.FILES, value=TRUE)  ## only public
##PGX.FILES = grep("sallusto|geiger|guarda", PGX.FILES, value=TRUE)
PGX.FILES

rda.file = PGX.FILES[7]
rda.file

for(rda.file in PGX.FILES) {

    ## load
    cat(">>> loading PGX file: ",rda.file,"\n")
    load(file=rda.file,verbose=1)
    dim(ngs$X)
    names(ngs)

    ## re-execute analysis parts
    source("../R/compute-extra2.R")

    ngs.save(ngs, file=rda.file, update.date=FALSE)


}

