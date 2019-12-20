##
## https://github.com/denalitherapeutics/archs4
##
##

##BiocManager::install("denalitherapeutics/archs4")
library("archs4")

archs4dir <- "~/.archs4data"
archs4dir <- "~/bigomics/data/archs4data"
archs4_local_data_dir_create(archs4dir)

cwd = getwd()
setwd(archs4dir)
system("wget https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5")
##system("wget https://s3.amazonaws.com/mssm-seq-matrix/human_hiseq_transcript_v2.h5")
system("wget https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix.h5")
##system("wget https://s3.amazonaws.com/mssm-seq-matrix/mouse_hiseq_transcript_v2.h5")
##system("wget ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz")
##system("wget ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz")

setwd(cwd)

## create_augmented_feature_info(archs4dir)
archs4_local_data_dir_validate()



library(archs4)
a4 <- Archs4Repository()
ids <- c('GSE89189', 'GSE29943', "GSM1095128", "GSM1095129", "GSM1095130")
sample.info <- sample_info(a4, ids)
head(sample.info)

yg <- as.DGEList(a4, "GSE89189", feature_type = "gene")
