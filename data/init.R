
RDIR="../R"
FILES="../lib"
PGX.DIR="."
dir.exists(PGX.DIR)

source(file.path(RDIR,"pgx-functions.R"), local=TRUE)  ## pass local vars
source(file.path(RDIR,"pgx-files.R"), local=TRUE)  ## pass local vars
pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=TRUE)
pgx.initDatasetFolder(PGX.DIR, force=FALSE, verbose=TRUE)

if(0) {

    info <- pgx.scanInfoFile(PGX.DIR, verbose=TRUE)
    head(info)[,1:2]

}
