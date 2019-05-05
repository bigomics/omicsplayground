source("../R/gx-heatmap.r")
library(Rtsne.multicore)
library(Rtsne)
library(ccRemover)
require(scran)
library(Seurat)

## Read in the expression matrix
X <- read.table(file = "../opt/seurat/nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)
X <- log2(1+X)
data("HScc_genes")
##data(t.cell_data)

X <- as.matrix(X)
dim(X)
X <- head(X[order(-apply(X,1,sd)),],8000)
dim(X)
head(X)


##------------------------------------------------------------------------------------
##------------------------------------ SEURAT ----------------------------------------
##------------------------------------------------------------------------------------

## Also read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../opt/seurat/regev_lab_cell_cycle_genes.txt")
s_genes <- cc.genes[1:43]
g2m_genes <- cc.genes[44:97]

# Create our Seurat object and complete the initalization steps
dim(X)
marrow <- CreateSeuratObject(raw.data = 2**X)
marrow <- NormalizeData(object = marrow)
##marrow <- FindVariableGenes(object = marrow, do.plot = FALSE, display.progress = FALSE)
##marrow <- ScaleData(object = marrow, display.progress = FALSE)
marrow <- CellCycleScoring(object = marrow,
                           s.genes = s_genes, g2m.genes = g2m_genes,
                           set.ident = TRUE)

## view cell cycle scores and phase assignments
head(x = marrow@meta.data)
table(marrow@meta.data$Phase)

##------------------------------------------------------------------------------------
##-------------------------------- SCRAN ---------------------------------------------
##------------------------------------------------------------------------------------

## Getting pre-trained marker sets
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
str(hs.pairs)

##
library(EnsDb.Hsapiens.v86)
tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_id", "gene_name"))
hgnc <- tx$gene_name
names(hgnc) <- tx$gene_id
hgnc.pairs <- lapply(hs.pairs, function(x)
    data.frame(first = hgnc[x[[1]]], second = hgnc[x[[2]]]))

# Classifying (note: test.data!=training.data in real cases)
assignments <- cyclone(X, hgnc.pairs)
head(assignments$scores)
table(assignments$phases)

# Visualizing
col1 <- factor(assignments$phases)
plot(assignments$score$G1, assignments$score$G2M, col=col1, pch=16)

##------------------------------------------------------------------------------------
##--------------------- ccRemover (Barron and Li) ------------------------------------
##------------------------------------------------------------------------------------
##biocLite("ccRemover")
library(ccRemover)
data("HScc_genes")
head(human_cell_cycle_genes)
head(mouse_cell_cycle_genes)

set.seed(10)
data(t.cell_data)
head(t.cell_data[,1:5])
## Center data and select small sample for example
X_cen <- t(scale(t(X[,]), center=TRUE, scale=FALSE))

## Extract gene names
gene_names <- rownames(X)
if_cc <- (gene_names %in% human_cell_cycle_genes$HGNC.symbol)
table(if_cc)
dat <- list(x=X_cen[,1:20], if_cc=if_cc)

## Run ccRemover (slow...)
xhat <- ccRemover(dat, cutoff = 3, max_it = 4, nboot = 200, ntop = 10)

##------------------------------------------------------------------------------------
##-------------------------------- using GSEA ----------------------------------------
##------------------------------------------------------------------------------------

require(GSVA)

load("../files/gmt-all.rda",verbose=1)
##gmt.cc <- gmt.all[grep("cell.cycle",names(gmt.all),ignore.case=TRUE)]
##gmt.cc <- gmt.all[grep("cell.cycle|g1.s|g2.m",names(gmt.all),ignore.case=TRUE)]
gmt.cc <- gmt.all[grep("cell.cycle|g1.s|g2.m|[-_ ]g1[-_ /]|[-_ ]g2[-_ /]|s.phase",names(gmt.all),ignore.case=TRUE)]
head(names(gmt.cc))
names(gmt.cc)

gsetX <- gsva( X[,], gmt.cc[], parallel.sz=8)
dim(gsetX)

g1s <- table(unlist(gmt.all[grep("g1.s",names(gmt.all),ignore.case=TRUE)]))
g2m <- table(unlist(gmt.all[grep("g2.m",names(gmt.all),ignore.case=TRUE)]))
g1s.genes <- head(names(sort(g1s,decreasing=TRUE)),200)
g2m.genes <- head(names(sort(g2m,decreasing=TRUE)),200)
g0.genes <- setdiff( unlist(gmt.cc), c(g1s.genes, g2m.genes))

g1s.g2m.common <- intersect(g1s.genes,g2m.genes)
g1s.genes <- setdiff(g1s.genes, g1s.g2m.common)
g2m.genes <- setdiff(g2m.genes, g1s.g2m.common)
length(g1s.genes)
length(g2m.genes)

gsm.markers <- list("GSM.MARKERS:G0/G1"=g0.genes,
                    "GSM.MARKERS:G2/M"=g2m.genes,
                    "GSM.MARKERS:G1/S"=g1s.genes)
gsm.phase <- gsva( X[,], gsm.markers, parallel.sz=8)
dim(gsm.phase)

gsm.idx <- max.col(t(gsm.phase))
seurat.phase <- marrow@meta.data$Phase
scran.phase  <- assignments$phases
table(seurat.phase, scran.phase)
table(seurat.phase, gsm.idx)

annot <- cbind( seurat=seurat.phase,
               cyclone=scran.phase,
               gsm=gsm.idx)

tt <- sort(table(unlist(gmt.cc)),decreasing=TRUE)
gmt.genes <- head(intersect(names(tt),rownames(X)),100)
xc <- X[gmt.genes,]
xc <- head(xc[order(-apply(xc,1,sd)),],40)
xc <- 1*(xc > apply(xc,1,median,na.rm=TRUE))
annot1 <- cbind(annot, t(xc))
dim(annot1)

jj <- unique(c(1,order(-apply(gsetX,1,sd))))
gs1 <- head(gsetX[jj,],80)
dim(gs1)

par(mfrow=c(1,1), oma=c(4,4,4,25))
gx.heatmap(gs1, mar=c(4,2), cexRow=0.8,
           col.annot=annot1, annot.ht=2.5,
           keysize=0.8, key=FALSE
)


pdf("test-cellcycle.pdf",w=12,h=15)
par(mfrow=c(1,1), oma=c(2,4,4,30))
gs2 <- gsetX[grep("g1.s|g2.m|[-_ ]g1[-_ /]|[-_ ]g2[-_ /]|s.phase",rownames(gsetX),ignore.case=TRUE),]
dim(gs2)
gx.heatmap(gs2, mar=c(4,2), cexRow=0.8,
           col.annot=annot1, annot.ht=1.5,
           keysize=0.5, key=FALSE
)
dev.off()


rho <- cor( t(gs1), t(xc))
par(mfrow=c(1,1), oma=c(2,4,4,40))
gx.heatmap(rho, mar=c(2,2), cexRow=0.8, ## scale="none",
           ##col.annot=annot1, annot.ht=2.5,
           keysize=0.8, key=FALSE
)




##------------------------------------------------------------------------------------
##-------------------------------- compare -------------------------------------------
##------------------------------------------------------------------------------------



seurat.phase <- marrow@meta.data$Phase
scran.phase  <- assignments$phases
table(seurat.phase, scran.phase)

s1 <- cbind(marrow@meta.data$S.Score,assignments$scores$S)
s2 <- cbind(marrow@meta.data$G2M.Score,assignments$scores$G2M)
meta.s   <- rowMeans(apply(s1,2,rank))
meta.g2m <- rowMeans(apply(s2,2,rank))

require(irlba)
jj  <- head(order(-apply(X,1,sd)),1000)
pos <- irlba(X[jj,],nv=2)$v

col1 <- factor(seurat.phase)
col2 <- factor(scran.phase)
plot(pos, col=col1)
plot(pos, col=col2)

plot(meta.s, meta.g2m, col=col1)
plot(meta.s, meta.g2m, col=col2)


devtools::install_github("varemo/piano")
library(piano)
data("gsa_input")
gsc <- loadGSC(gsa_input$gsc)
gsares <- runGSA(geneLevelStats=gsa_input$pvals , directions=gsa_input$directions,
                 gsc=gsc, nPerm=500)
networkPlot(gsares,class="non",significance=0.01)
nw <- networkPlot(gsares,class="non",significance=0.01,lay=5)
networkPlot(gsares,class="distinct",direction="both",lay=nw$layout,geneSets=nw$geneSets)
