##
##
## NOTE: This file is supposed to run in the folder .../R/
##

##install.packages("devtools")
##install.packages("BiocManager")
if(!require(devtools)) install.packages("devtools")
if(!require(BiocManager)) install.packages("BiocManager")

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
    "BH","PCSF",
    "DeMixT", ## purify
    "RNAseqData.HNRNPC.bam.chr14",
    "org.Mm.eg.db",
    "tximportData",
    "EnsDb.Hsapiens.v86",
    "EnsDb.Mmusculus.v79",
    "TxDb.Hsapiens.UCSC.hg19.knownGene",
    "TxDb.Mmusculus.UCSC.mm10.knownGene"
)    

PKG.MANUAL <- c(
    "gputools","Seurat","EPIC","PCSF","NNLM","iTALK",
    "fpc","grid","gridGraphics","Rgraphviz","rWordCloud",
    "shinyparticles","FastGGM","monocle3","proteus",
    "fastcluster","mygene","diptest","infercnv","pathview")

##---------------------------------------------------------------------
## Install base packages
##---------------------------------------------------------------------

base.pkg = c("shiny","flexdashboard","shinydashboard")
install.pkgs(base.pkg)

## priority
install.pkgs("randomcoloR")

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
install.pkg("randomForest", force=TRUE)

remove.pkg("fpc")
install.pkgs(c('mclust', 'flexmix', 'prabclus', 'diptest', 'mvtnorm', 'robustbase', 'kernlab', 'trimcluster'))
##install.packages("../ext/packges/fpc_2.1-10.tar.gz",repos=NULL,type="source")
install_version("fpc", version="2.1-10", repos="http://cran.us.r-project.org")

install.pkg("webshot")
webshot::install_phantomjs()

##---------------------------------------------------------------------
## Install latest from GITHUB
##---------------------------------------------------------------------
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
##devtools::install_github("IOR-Bioinformatics/PCSF", dependencies=TRUE, type="source")
devtools::install_github('linxihui/NNLM')
devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
devtools::install_github('adymimos/rWordCloud', force=TRUE)
remotes::install_github("dreamRs/shinyparticles")

##---------------------------------------------------------------------
## ONLY DEV.VERSION
##---------------------------------------------------------------------
if(0) {

    remotes::install_github("trevorld/r-argparse")
    devtools::install_github("wwylab/DeMixT")
    
    ## ---- monocle3 (only DEV!!! many install problems in R 3.5.2!!!)
    install.pkg("uwot", force=TRUE)
    BiocManager::install("batchelor")
    install.pkgs(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                   'S4Vectors', 'SingleCellExperiment','SummarizedExperiment'))
    devtools::install_github('cole-trapnell-lab/leidenbase')
    devtools::install_github('cole-trapnell-lab/monocle3')
}

## ----- proteus
## devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=FALSE)

##---------------------------------------------------------------------
## Seurat needs to be downgraded and dependencies to installed...
##---------------------------------------------------------------------
install.pkgs(c("ROCR", "mixtools", "lars", "ica", "tsne", "ape", "dtw", "SDMTools", "ggridges", "fitdistrplus", "doSNOW","diffusionMap","fpc","hdf5r"))
install.pkgs(c('cowplot', 'Rtsne', 'pbapply', 'RANN', 'dplyr', 'irlba', 'plotly', 'Hmisc', 'tidyr', 'metap', 'lmtest', 'png', 'reticulate', 'RcppEigen', 'RcppProgress'))
install.pkg("RcppEigen",force=TRUE)
install.packages("../ext/packages/Seurat_v2.3.3.tar.gz",repos=NULL,type="source")  ## old version

##---------------------------------------------------------------------
## make sure local ones are preferred
##---------------------------------------------------------------------
install.packages("../ext/packages/nclust1_1.9.4.tar.gz",repos=NULL,type="source")
install.packages("../ext/packages/nclust_2.1.1.tar.gz",repos=NULL,type="source")

install.pkgs(c("RcppParallel"))
install.packages("../ext/packages/FastGGM.tar.gz", repos = NULL, type = "source")

install.pkgs(c("HiddenMarkov","coin","rjags","future","argparse"))
install.pkg("coin", force=TRUE)
install.packages("../ext/packages/infercnv_1.1.3mod.tar.gz",repos=NULL,type="source")  ## old version

install.pkg("KEGGREST", force=TRUE)
install.packages("../ext/packages/pathview_1.16.7.tar.gz",repos=NULL,type="source")  ## old version

##---------------------------------------------------------------------
## remove unneccessary big packages...
##---------------------------------------------------------------------

remove.pkgs(BIG.NOTUSED)
