## Scans files in PGX.DIR (default="../pgx") and updates
##
##

PGX.DIR="../pgx"

cat("scanning PGX file info...\n")
PGXINFO <- pgx.scanInfo(pgx.dir=PGX.DIR, inc.progress=FALSE)
cat("dim.PGXINFO=",dim(PGXINFO),"\n")
INFO.FILE <- "./pgx-datasets-info.csv"
##Sys.chmod("../shiny-dev", mode="0777")
Sys.chmod(INFO.FILE, mode="0666")
write.csv(PGXINFO, file=INFO.FILE)
