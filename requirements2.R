##source("http://bioconductor.org/biocLite.R")

if(!require(devtools)) install.packages("devtools")
if(!require(BiocManager)) install.packages("BiocManager")
install.pkg <- function(pkg) {
    if(!pkg %in% installed.packages()) {
        require(BiocManager)
        BiocManager::install(pkg, dependencies=TRUE,
                             ask=FALSE, update=FALSE)
    }
}
install.pkgs <- function(pkgs) {
    for(pkg in pkgs) install.pkg(pkg)
}
remove.pkg <- function(pkg) {
    if(pkg %in% installed.packages()) remove.packages(pkg)
}


##---------------------------------------------------------------------
## CRAN packages
##---------------------------------------------------------------------

install.pkg("fpc")

##---------------------------------------------------------------------
## from GITHUB
##---------------------------------------------------------------------
##devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)

##---------------------------------------------------------------------
## from local folder
##---------------------------------------------------------------------
##install.packages("ext/nclust1_1.9.4.tar.gz",repos=NULL,type="source")

##---------------------------------------------------------------------
## remove unneccessary Big Shit...
##---------------------------------------------------------------------
remove.pkg("reactome.db")  ## >2GB!!
remove.pkg("BH")  ## boost header files
remove.pkg("RNAseqData.HNRNPC.bam.chr14")
remove.pkg("EnsDb.Hsapiens.v86")
remove.pkg("org.Mm.eg.db")
remove.pkg("tximportData")
remove.pkg("TxDb.Hsapiens.UCSC.hg19.knownGene")


