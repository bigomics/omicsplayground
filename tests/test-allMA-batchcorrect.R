source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")
source("../R/xcr-graph.r")
source("../R/pgx-graph.R")
source("../R/pgx-correct.R")

library(Rtsne.multicore)
library(Rtsne)
library(qlcMatrix)
library(corpora)
##devtools::install_github("bwlewis/rthreejs")

load(file="../files/allMA.rda", verbose=1)

dim(allM)
head(allM)[,1:4]

batch <- sub("\\].*","]",colnames(allM))
table(batch)

pgx.performBatchCorrection


rho <- cor(allM, use="pairwise")
dim(rho)
sum(is.na(rho))
rho[is.na(rho)] <- 0

dist <- (1 - rho)**4
pos <- Rtsne( dist, is_distance=TRUE, num_thread=2)$Y




