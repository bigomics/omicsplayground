
## scan and update datasets info
source("../R/pgx-functions.R")
cat("scanning PGX file info...\n")
pgxinfo <- pgx.scanInfo(pgx.dir="../data", inc.progress=FALSE)
cat("dim.pgxinfo=",dim(pgxinfo),"\n")
info.file <- "../data/datasets-info.csv"
write.csv(pgxinfo, file=info.file)
Sys.chmod(info.file, mode="0666")
