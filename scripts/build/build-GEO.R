##
## This builds GEO datasets
##
##

source("../../R/pgx-geoquery.R")
source("../../R/pgx-functions.R")
source("../../R/ngs-functions.R")

platforms = c(
    "GPL570", ## Affymetrix HG-U133plus2
    "GPL6244", ## Affymetrix HuGene 1.0st
    "GPL10558", ## Illumina HT-12 v4.0
    "GPL11154", ## Illumina HiSeq 2000
    "GPL16791", ## Illumina HiSeq 2500
    "GPL20401" ## Illumina HiSeq 4000
)
GEODIR="geo"

## BiocManager::install("GEOmetadb")    
library(GEOmetadb)
if(!file.exists(file.path(GEODIR,'GEOmetadb.sqlite'))) getSQLiteFile(destdir=GEODIR)

con <- dbConnect(SQLite(),file.patch(GEODIR,'GEOmetadb.sqlite'))
##dbDisconnect(con)
geo_tables <- dbListTables(con)
geo_tables

rs <- dbGetQuery(con,paste("select gse,title from gse where",
                           "contributor like '%Ivo%Kwee%'",sep=" "))
rs

platforms2 <- paste0("('",paste(platforms,collapse="','"),"')")
sql <- paste("SELECT DISTINCT gse.title,gse.gse",
             "FROM",
             "  gsm JOIN gse_gsm ON gsm.gsm=gse_gsm.gsm",
             "  JOIN gse ON gse_gsm.gse=gse.gse",
             "  JOIN gse_gpl ON gse_gpl.gse=gse.gse",
             "  JOIN gpl ON gse_gpl.gpl=gpl.gpl",
             "WHERE",
             "  gsm.molecule_ch1 like '%total RNA%' AND",
             "  gsm.characteristics_ch1 is not null AND",             
             ## "  gse.title LIKE '%breast cancer%' AND",
             "  gpl.gpl IN ",platforms2," AND",             
             "  gpl.organism LIKE '%Homo sapiens%'",sep=" ")
rs <- dbGetQuery(con,sql)
dim(rs)
sel1 <- rs$gse  ## datasets with phenotype information

## Select GSE datasets
library(dplyr)
db = src_sqlite(file.path(GEODIR,'GEOmetadb.sqlite'))
gse = data.frame(tbl(db,'gse'))
gpl = data.frame(tbl(db,'gse_gpl'))
gsm = data.frame(tbl(db,'gse_gsm'))
gsm0 = tbl(db,'gsm')
nsamples = table(gsm$gse)
nsamples <- nsamples[match(gpl$gse,names(nsamples))]

s1 <- ( gpl$gpl %in% platforms )
s2 <- ( nsamples >= 20 & nsamples <= 200 )
sel2 <- gpl$gse[ which(s1 & s2) ]  ## dataset with enought samples and in selected GPL
length(sel2)    


gse1 <- gse[ which(gse$gse %in% sel1 & gse$gse %in% sel2),]
dim(gse1)
head(gse1)[,1:3]
tail(gse1)[,1:3]
grep("[ ]aging",gse1$title,ignore.case=TRUE,value=TRUE)
grep("[ ]prostate",gse1$title,ignore.case=TRUE,value=TRUE)

library(GEOquery)
id = gse1$gse[1]
id = "GSE21653"  ## BRCA
id = "GSE138126" ## 
id = "GSE10846"  ## DLBCL
id = "GSE108467"
##id = "GDS507"
id = "GSE53784"

## examples
sel.ids = tail(gse1$gse,10) ## newer first..
id = sel.ids[1]
for(id in sel.ids) {

    cat(">>>>>>>>>>>>> Quering GEO dataset",id,"<<<<<<<<<<<<\n")
    geo <- NULL
    fn1 = file.path(GEODIR,paste0("geo-",id,".error"))
    fn2 = file.path(GEODIR,paste0("geo-",id,".rds"))
    if( file.exists(fn1) || file.exists(fn2)) {
        ## cat(">>> ALREADY EXISTS:: ",id,"\n")
        next
    }
    geo <- pgx.GEOquery(id)
    is.null(geo$error)
    if(!is.null(geo$error)) {
        cat(">>> ",geo$error,"\n")
        write(geo$error, file=fn1)
    } else {
        lapply(geo, dim)
        saveRDS(geo, file=fn2)
    }
}


source("../../R/ngs-functions.R")

rds.list <- dir(GEODIR,pattern=".rds$",full.names=TRUE)
rds <- rds.list[1]
for(rds in rds.list) {
    cat(">>>>>>> dataset",rds,"\n")
    geo <- readRDS(rds)
    lapply(geo, dim)
    cat("CONTRASTS: ",paste(colnames(geo$contrasts),collapse=", "),"\n")


    
    

}




