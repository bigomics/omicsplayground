
RDIR="../R"
FILES="../lib"
PGX.DIR="."
dir.exists(PGX.DIR)

source(file.path(RDIR,"pgx-functions.R"), local=TRUE)  ## pass local vars
source(file.path(RDIR,"pgx-api.R"), local=TRUE)  ## pass local vars
source(file.path(RDIR,"pgx-files.R"), local=TRUE)  ## pass local vars
pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)
##pgx.initDatasetFolder(PGX.DIR, force=FALSE, verbose=1)
