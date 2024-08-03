# This file is supposed to run from the root Playground folder
if (basename(getwd()) != "omicsplayground") {
    stop("Please run from the OmicsPlayground root folder")
}

if(!require("renv")) install.packages("renv")
if(!require("BiocManager")) install.packages("BiocManager")
if(!require("remotes")) install.packages("remotes")
if(!require("devtools")) install.packages("devtools")
if(!require("reticulate")) install.packages("reticulate")

## ---------------------------------------------------------------------
## Automatically scan all used packages and install
## ---------------------------------------------------------------------
## We use renv to detect dependencies. Renv is looking for library and
## require statements in the r/R source files.
renv.out <- renv::dependencies(path = "components", root = getwd(), errors = "ignored")
pkg.used <- unique(renv.out$Package)

## Define remote locations or versions
github_url <- function(repo) paste0("github::",repo)
github_url <- function(repo) {
  if(grepl("@",repo)) {
    branch <- sub(".*@","",repo)
    repo <- sub("@.*","",repo)    
    paste0("url::https://github.com/",repo,"/archive/refs/heads/",branch,".zip")
  } else {
    paste0("url::https://github.com/",repo,"/archive/HEAD.zip")
  }
}

remotes.url <- c()
add_github <- function(repo) {
  pkg.name <- gsub(".*[/]|@.*","",repo)
  remotes.url[pkg.name] <<- github_url(repo)
}
add_github("bigomics/PCSF")
add_github("bigomics/playdata")
add_github("bigomics/playbase")
add_github("bigomics/bigdash")
add_github("bigomics/bigLoaders")
add_github("bigomics/fgsea")
add_github("bigomics/wizardR")
add_github("GfellerLab/EPIC")
add_github("broadinstitute/infercnv")
add_github("GfellerLab/SuperCell")
add_github("linxihui/NNLM")
add_github("Coolgenome/iTALK")
add_github("wt2015-github/FastGGM")
add_github("satijalab/azimuth")
add_github("JohnCoene/waiter")
add_github("JohnCoene/firebase@omics")
add_github("JohnCoene/bsutils")
add_github("ropensci/iheatmapr")
##add_github("rstudio/bslib@v0.6.1")
add_github("rstudio/htmltools")
add_github("bigomics/biomaRt")
add_github("Bioconductor/BiocFileCache")
add_github("cysouw/qlcMatrix")
add_github("cole-trapnell-lab/leidenbase")
add_github('cole-trapnell-lab/monocle3')
add_github('bartongroup/Proteus')
add_github('cran/riverplot')
add_github('Ironholds/rgeolocate')

#remotes::install_version("Matrix", version = "1.6.0")
#install.packages("plsRcox", dependencies = FALSE) # pscRcox requires Matrix >= 1.6.0

pkg.remotes <- sort(remotes.url[names(remotes.url) %in% pkg.used])
pkg.imports <- sort(setdiff(pkg.used, names(pkg.remotes)))

## determine packages in playbase
playbase.imports <- strsplit(packageDescription("playbase")$Imports,split='\\n')[[1]]
playbase.imports <- gsub(",$","",playbase.imports)
playbase.remotes <- strsplit(packageDescription("playbase")$Remotes,split='\\n')[[1]]
playbase.remotes <- gsub(".*github.com/|/archive.*","",playbase.remotes)
playbase.remotes <- gsub(".*contrib/|_[0-9].*$","",playbase.remotes)
playbase.remotes <- sub(".*/","",playbase.remotes)                    
playbase.pkg <- sort(unique(c(playbase.imports, playbase.remotes)))

bigomics.pkg <- c("playbase","playdata","bigdash","bigLoaders","wizardR")

## we only report packages that are not in playbase
table(pkg.imports %in% playbase.pkg)
table(names(pkg.remotes) %in% playbase.pkg)
pkg.imports <- setdiff(pkg.imports, playbase.pkg)

sel1 <- !names(pkg.remotes) %in% playbase.pkg
sel2 <- names(pkg.remotes) %in% bigomics.pkg
pkg.remotes <- sort(pkg.remotes[sel1 | sel2])

pkg.imports
names(pkg.remotes)

if(file.exists("DESCRIPTION")) {
  message("WARNING: overwriting existing DESCRIPTION file.")
}
desc.header <- readLines("dev/description.app")
desc.file <- "DESCRIPTION"
write(desc.header, file=desc.file)
write("Imports:", file=desc.file, append=TRUE)
write(paste0("    ",pkg.imports,","), file=desc.file, append=TRUE)
write("Remotes:", file=desc.file, append=TRUE)
write(paste0("    ",pkg.remotes,","), file=desc.file, append=TRUE)

## ---------------------------------------------------------------------
## Installing packages
## ---------------------------------------------------------------------

## some are perhaps already installed by playbase or other packages
installed.pkg <- sort(as.character(installed.packages()[,"Package"]))
missing.imports <- setdiff( pkg.imports, installed.pkg )
missing.remotes <- pkg.remotes[!names(pkg.remotes) %in% installed.pkg]

print(">>> Installing missing CRAN/BioConductor packages...")
message( ">>> Installing ",length(missing.imports),
        " missing CRAN/BioConductor packages:\n",
        paste(missing.imports,collapse="\n"))
if(length(missing.imports>0)) {
  BiocManager::install( missing.imports, ask=FALSE, force=FALSE )
}

message(">>> Installing ",length(missing.remotes)," missing remote packages")
if(length(missing.remotes>0)) {
  missing.remotes <- sub("^url::","",missing.remotes)
  for(url in missing.remotes) {
    message("Installing remote: ", url)
    remotes::install_url( url, ask=FALSE, force=FALSE )
  }
}

## ---------------------------------------------------------------------
## Install Kaleido for plotly
## ---------------------------------------------------------------------

print(">>> Installing Kaleido/plotly packages...")
## Install a clean reticulate and miniconda
# install.packages('reticulate', force=TRUE) # remove reticulate install since its already done.. and we get checksum error for some reason at this step
unlink("~/.local/share/r-miniconda", recursive = TRUE)
reticulate::install_miniconda()
reticulate::conda_install("r-reticulate", "python-kaleido")
reticulate::conda_install("r-reticulate", "plotly", channel = "plotly")
reticulate::use_miniconda("r-reticulate")

## ---------------------------------------------------------------------
## remove unneccessary big packages??
## ---------------------------------------------------------------------
print(">>> removing not needed packages...")

BIG.NOTUSED <- c(
    "reactome.db", ## >2GB!!!
    "terra",
    ## "DeMixT", ## purify
    "RNAseqData.HNRNPC.bam.chr14",
    ## "RSpectra",  ## ???
    "tximportData"
    ## "EnsDb.Hsapiens.v86",
    ## "EnsDb.Mmusculus.v79",
    ## "TxDb.Hsapiens.UCSC.hg19.knownGene",  ## need for import
    ## "TxDb.Mmusculus.UCSC.mm10.knownGene"  ## need for import
)
## remove.pkgs(BIG.NOTUSED)


message("**** FINISHED INSTALLING REQUIREMENTS *****")
