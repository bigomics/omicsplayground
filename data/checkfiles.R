##
##
##

pgx.files <- dir(".", pattern=".pgx$")

fields <- c('samples','counts','genes','tsne2d','tsne3d','X','model.parameters','gx.meta','gset.meta','gsetX','GMT')

pn.list <- c()
ok.list <- c()
p=pgx.files[1]
for(p in pgx.files[]) {

    try.error <- try(load(p, verbose=0))
    ok=TRUE
    if(class(try.error)=="try-error") ok=FALSE
    ok <- ok && all(fields %in% names(ngs))
    ok
    if(!ok) {
        cat("loading ",p,": ERROR\n")
        ##confirm <- askYesNo("Do you want to disable this file?")
        ##if(confirm)
        file.rename(p, paste0(p,"_BROKEN"))
    } else {
        cat("loading ",p,": OK\n")
    }
    ok.list[p] <- ok
    pn.list[p] <- is.POSvsNEG(ngs)
}

if(0) {

    table(ok.list)
    table(pn.list)
    cbind(ok.list)

}
