
RDIR="../R"
FILES="../lib"
PGX.DIR="."
dir.exists(PGX.DIR)

source(file.path(RDIR,"pgx-functions.R"), local=TRUE)  ## pass local vars
source(file.path(RDIR,"pgx-api.R"), local=TRUE)  ## pass local vars
source(file.path(RDIR,"pgx-files.R"), local=TRUE)  ## pass local vars

pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)
pgx.initDatasetFolder(PGX.DIR, force=FALSE, verbose=1)

##---------------------------------------------------------------
## create local signatures H5 file
##---------------------------------------------------------------
if(0) {
    source(file.path(RDIR,"pgx-include.R"), local=TRUE)  ## pass local vars    
    source(file.path(RDIR,"pgx-signature.R"), local=TRUE)  ## pass local vars

    h5.file = "./sigdb.h5"
    pgx.files <- dir(".", pattern="pgx$")
    pgx.files
    ## pgx.files <- head(pgx.files,10)  ## small test
    pgx.createSignatureDatabaseH5( pgx.files, h5.file, update.only=FALSE)
    h5ls(h5.file)
    pgx.addEnrichmentSignaturesH5( h5.file, mc.cores=4, lib.dir=FILES, methods=c("gsea"))
}
