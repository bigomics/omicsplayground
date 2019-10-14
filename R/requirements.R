##
##
## NOTE: This file is supposed to run in the folder .../R/
##

if(!require(devtools)) install.packages("devtools")
if(!require(BiocManager)) install.packages("BiocManager")
require(devtools)
require(BiocManager)

LOCAL.PKGS <- sub("_.*","",dir("../ext/packages"))
LOCAL.PKGS

install.pkg <- function(pkg, force=FALSE) {
    if(force || !pkg %in% installed.packages()) {

        if(pkg %in% LOCAL.PKGS) {
            ## if available locally, we install local version
            cat("installing",pkg,"from local folder...\n")
            pkg1 = dir("../ext/packages",pattern=paste0(pkg,"_"),full.name=TRUE)
            try( install.packages(pkg1,repos=NULL,type="source") )
        } else {
            cat("installing",pkg,"from CRAN/BioConductor...\n")
            try( BiocManager::install(pkg, dependencies=NA,
                                      ask=FALSE, update=FALSE))
        }
    } else {
        cat("package",pkg,"already installed\n")
    }
}
install.pkgs <- function(pkgs, force=FALSE) {
    for(pkg in pkgs) install.pkg(pkg, force=force)
}
remove.pkg <- function(pkg) {
    if(pkg %in% installed.packages()) remove.packages(pkg)
}
remove.pkgs <- function(pkgs, force=FALSE) {
    for(pkg in pkgs) remove.pkg(pkg)
}

BIG.NOTUSED <- c(
    "reactome.db", ## >2GB!!
    "BH",
    "RNAseqData.HNRNPC.bam.chr14",
    "org.Mm.eg.db",
    "tximportData",
    "EnsDb.Hsapiens.v86",
    "EnsDb.Mmusculus.v79",
    "TxDb.Hsapiens.UCSC.hg19.knownGene",
    "TxDb.Mmusculus.UCSC.mm10.knownGene"
)    

PKG.MANUAL <- c("gputools","Seurat","EPIC","PCSF","NNLM","iTALK",
                "fpc","grid","gridGraphics","Rgraphviz",
                "fastcluster","mygene","diptest","infercnv")

##---------------------------------------------------------------------
## Install base packages
##---------------------------------------------------------------------

base.pkg = c("shiny","flexdashboard","shinydashboard")
install.pkgs(base.pkg)

##---------------------------------------------------------------------
## Automatically scan all used packages and install
##---------------------------------------------------------------------

pkg.used <- system("grep 'library(\\|require(' *R *r ../shiny/*Rmd", intern=TRUE)
pkg.used <- gsub(".*require\\(|.*library\\(","",pkg.used)
pkg.used <- gsub("\"|\\).*","",pkg.used)
pkg.used <- grep("[ ]|quietly",pkg.used,value=TRUE,invert=TRUE)
pkg.used <- sort(unique(pkg.used))

install.pkgs( setdiff(pkg.used, c(PKG.MANUAL,BIG.NOTUSED)) )


##---------------------------------------------------------------------
## reinstall problematics ones
##---------------------------------------------------------------------

require(devtools)
##install.packages("gridGraphics")
##install.pkg("grid", force=TRUE)
install_version("gridGraphics", version="0.3-0", repos="http://cran.us.r-project.org")
install.pkg("Rgraphviz", force=TRUE)
install.pkg("fastcluster", force=TRUE)
install.pkg("mygene", force=TRUE)
install.pkg("diptest", force=TRUE)

remove.pkg("fpc")
install.pkgs(c('mclust', 'flexmix', 'prabclus', 'diptest', 'mvtnorm', 'robustbase', 'kernlab', 'trimcluster'))
##install.packages("../ext/packges/fpc_2.1-10.tar.gz",repos=NULL,type="source")
install_version("fpc", version="2.1-10", repos="http://cran.us.r-project.org")

##---------------------------------------------------------------------
## Install latest from GITHUB
##---------------------------------------------------------------------
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
devtools::install_github("IOR-Bioinformatics/PCSF",
                         ## repos=BiocInstaller::biocinstallRepos(),
                         dependencies=TRUE, type="source", force=TRUE)
devtools::install_github('linxihui/NNLM')
devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
## devtools::install_github("broadinstitute/infercnv", ref="RELEASE_3_9")
devtools::install_github('adymimos/rWordCloud')

##---------------------------------------------------------------------
## Seurat needs to be downgraded and dependencies to installed...
##---------------------------------------------------------------------
install.pkgs(c("ROCR", "mixtools", "lars", "ica", "tsne", "ape", "dtw", "SDMTools", "ggridges", "fitdistrplus", "doSNOW","diffusionMap","fpc","hdf5r"))
install.pkgs(c('cowplot', 'Rtsne', 'pbapply', 'RANN', 'dplyr', 'irlba', 'plotly', 'Hmisc', 'tidyr', 'metap', 'lmtest', 'png', 'reticulate', 'RcppEigen', 'RcppProgress'))
install.packages("../ext/packages/Seurat_v2.3.3.tar.gz",repos=NULL,type="source")  ## old version

##---------------------------------------------------------------------
## make sure local ones are preferred
##---------------------------------------------------------------------
install.packages("../ext/packages/nclust1_1.9.4.tar.gz",repos=NULL,type="source")
install.packages("../ext/packages/nclust_2.1.1.tar.gz",repos=NULL,type="source")

install.pkgs(c("HiddenMarkov","coin","rjags","future","argsparse"))
install.packages("../ext/packages/infercnv_1.1.3mod.tar.gz",repos=NULL,type="source")  ## old version

install.pkg("KEGGREST", force=TRUE)
install.packages("../ext/packages/pathview_1.16.7.tar.gz",repos=NULL,type="source")  ## old version

##---------------------------------------------------------------------
## remove unneccessary big packages...
##---------------------------------------------------------------------
remove.pkgs(BIG.NOTUSED)
