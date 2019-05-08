##source("http://bioconductor.org/biocLite.R")

if(!require(devtools)) install.packages("devtools")
if(!require(BiocManager)) install.packages("BiocManager")
install.pkg <- function(pkg, force=FALSE) {
    if(force && (pkg %in% installed.packages())) {
        ##remove.packages(pkg)
    }
    if(force || !pkg %in% installed.packages()) {
        require(BiocManager)
        cat("installing",pkg,"\n")
        BiocManager::install(pkg, dependencies=NA,
                             ask=FALSE, update=FALSE)
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

##---------------------------------------------------------------------
## CRAN packages
##---------------------------------------------------------------------
install.pkg("Rgraphviz", force=TRUE)

##---------------------------------------------------------------------
## Bioconductor packages
##---------------------------------------------------------------------
install.pkg("Rgraphviz")

##---------------------------------------------------------------------
## from GITHUB
##---------------------------------------------------------------------
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
devtools::install_github("IOR-Bioinformatics/PCSF",
                         ## repos=BiocInstaller::biocinstallRepos(),
                         dependencies=TRUE, type="source", force=TRUE)
devtools::install_github('linxihui/NNLM')
devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)

##---------------------------------------------------------------------
## from local folder
##---------------------------------------------------------------------

remove.pkg("fpc")
install.pkgs(c('mclust', 'flexmix', 'prabclus', 'diptest', 'mvtnorm', 'robustbase', 'kernlab', 'trimcluster'))
install.packages("ext/fpc_2.1-10.tar.gz",repos=NULL,type="source")
install.packages("ext/nclust1_1.9.4.tar.gz",repos=NULL,type="source")
install.packages("ext/nclust_2.1.1.tar.gz",repos=NULL,type="source")
install.packages("ext/pathview_1.16.7.tar.gz",repos=NULL,type="source")  ## old version
install.packages("ext/FARDEEP_1.0.1.tar.gz",repos=NULL,type="source")  ## old version
##install.packages("ext/gputools_1.1.tar.gz",repos=NULL,type="source")  ## old version

## Seurat needs to be downgraded and dependencies to installed...
install.pkgs(c("ROCR", "mixtools", "lars", "ica", "tsne", "ape", "dtw", "SDMTools", "ggridges", "fitdistrplus", "doSNOW","diffusionMap","fpc","hdf5r"))
install.pkgs(c('cowplot', 'Rtsne', 'pbapply', 'RANN', 'dplyr', 'irlba', 'plotly', 'Hmisc', 'tidyr', 'metap', 'lmtest', 'png', 'reticulate', 'RcppEigen', 'RcppProgress'))
install.packages("ext/Seurat_v2.3.3.tar.gz",repos=NULL,type="source")  ## old version

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


