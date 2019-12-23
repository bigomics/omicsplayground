
RDIR="../R"
FILES="../lib"
PGX.DIR="."
dir.exists(PGX.DIR)

source("../R/pgx-functions.R", local=TRUE)  ## pass local vars
source("../R/pgx-files.R", local=TRUE)  ## pass local vars
pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)
##pgx.initDatasetFolder(PGX.DIR, force=FALSE, verbose=1)
