source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")

## https://satijalab.org/seurat/pbmc3k_tutorial.html
##
##

load(file="../pgx/sallusto2019b-tenx-s1000-4k.pgx",verbose=1)
load(file="../pgx/sallusto2019b-tenx-s2000-8k.pgx",verbose=1)
load(file="../pgx/tenx-pbmc1k-8k.pgx",verbose=1)
##load(file="../pgx/GSE98638-liver-scRNA-8k.pgx",verbose=1)

library(Seurat)
##library(scran)
##library(scater)

CANONICAL.MARKERS = list(
    "CD4 T cells"=c("IL17R"),
    "CD14+ Monocytes"=c("CD14","LYZ"),
    "B cells"=c("MS4A1"),
    "CD8 T cells"="CD8A",    
    "FCGR3A+ Monocytes"=c("FCGR3A","MS4A7"),
    "NK cells"=c("GNLY","NKG7"),
    "Dendritic cells"=c("FCER1A","CST3"),
    "Megakaryocytes"=c("PPBP")
)

##-------------- create Seurat object from counts -----------------------

obj <- CreateSeuratObject(raw.data = ngs$counts, min.cells = 3,
                           min.genes = 200, project = "mydata_scRNAseq")
obj
head(ngs$genes)
sum(is.na(ngs$counts))
head(ngs$samples)

##obj <- AddMetaData(object = obj, metadata = ngs$samples$cell.type, col.name = "cell.type")
obj <- AddMetaData(object = obj, metadata = ngs$samples$.cell_cycle, col.name = ".cell_cycle")
obj <- AddMetaData(object = obj, metadata = ngs$samples$treatment, col.name = "treatment")

##sel <- head(intersect(rownames(ngs$genes), unlist(CANONICAL.MARKERS)),3)
VlnPlot(object = obj, features.plot = c("GZMB","CCL5","CD4"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = obj, gene1 = "GZMB", gene2 = "CCL5")
GenePlot(object = obj, gene1 = "GZMB", gene2 = "nGene")

##-------------- QC, mito genes -----------------------

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = obj@data), value = TRUE)
percent.mito <- Matrix::colSums(obj@raw.data[mito.genes, ])/Matrix::colSums(obj@raw.data)
obj <- AddMetaData(object = obj, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = obj, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1,2))
GenePlot(object = obj, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = obj, gene1 = "nUMI", gene2 = "nGene")

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold.
obj <- FilterCells(object = obj, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(3500, 0.05))

##-------------- normalization, filtering, correction -----------------------
mean(colSums(obj@raw.data))
obj <- NormalizeData(object = obj, normalization.method = "LogNormalize",
                      scale.factor = 10000)

##-------------- detection of variable genes -----------------------
obj <- FindVariableGenes(object = obj, mean.function = ExpMean,
                          dispersion.function = LogVMR, x.low.cutoff = 0.0125,
                          x.high.cutoff = 3, y.cutoff = 0.5)
length(x = obj@var.genes)

##-------------- scale and simulataneous batch correction
##obj <- ScaleData(object = obj)
obj <- ScaleData(object = obj, vars.to.regress = c("nUMI", "percent.mito"))

##-------------------------- PCA --------------------------------------
## Examine and visualize PCA results a few different ways
obj <- RunPCA(object = obj, pcs.compute = 100,
              pc.genes = obj@var.genes, do.print = TRUE,
              pcs.print = 1:5, genes.print = 5)

##PrintPCA(object = obj, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
##VizPCA(object = obj, pcs.use = 1:2)
PCAPlot(object = obj, dim.1 = 1, dim.2 = 2)

## show the PC vectors as heatmaps
obj <- ProjectPCA(object = obj, do.print = FALSE)
PCHeatmap(object = obj, pc.use = 1:6, cells.use = 500, do.balanced = TRUE,
    label.columns = FALSE, use.full = FALSE)

##-------------------------- clustering and tSNE --------------------------------------
# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
PCElbowPlot(object = obj)
obj <- FindClusters(object = obj, reduction.type = "pca", dims.use = 1:20,
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = obj)
obj <- RunTSNE(object = obj, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = obj)

saveRDS(obj, file = "seurat-object.rds")

##-------------------------------------------------------------------------------

# find *all* markers compared to all remaining cells
obj.markers <- FindAllMarkers(object = obj, only.pos = TRUE,
                               min.pct = 0.25, thresh.use = 0.25)
toptable <- obj.markers %>% group_by(cluster) %>% top_n(2,avg_logFC)
toptable

## Violin plot
top.markers <- toptable$gene
VlnPlot(object = obj, features.plot = top.markers, point.size.use=0.3)
VlnPlot(object = obj, features.plot = top.markers, point.size.use=0.3,
        use.raw=TRUE, y.log=TRUE)

##-------------------------------------------------------------------------------

## Feature plot
markers <- unlist(CANONICAL.MARKERS)
markers <- intersect(markers, rownames(ngs$counts))
markers
FeaturePlot(object = obj, features.plot = markers,
            cols.use = c("grey", "blue"), reduction.use = "tsne")

## Heatmap plot
require(tidyverse)
top10 <- obj.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top10
DoHeatmap(object = obj, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)


##=====================================================================================
##================================ CELL CYCLE =========================================
##=====================================================================================
library(Seurat)

load(file="../pgx/GSE72056-melanoma-scRNA-vsphenotype-TC8k-LT.pgx",verbose=1)
exp.mat <- round(2**ngs$X)

# Read in the expression matrix The first row is a header row, the first
# column is rownames
exp.mat <- read.table(file = "../opt/seurat/nestorawa_forcellcycle_expressionMatrix.txt",
    header = TRUE, as.is = TRUE, row.names = 1)

dim(exp.mat)
head(exp.mat)[,1:4]

# Also read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../opt/seurat/regev_lab_cell_cycle_genes.txt")

# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

# Create our Seurat object and complete the initalization steps
dim(exp.mat)
marrow <- CreateSeuratObject(raw.data = exp.mat)
marrow <- NormalizeData(object = marrow)
marrow <- FindVariableGenes(object = marrow, do.plot = FALSE, display.progress = FALSE)
marrow <- ScaleData(object = marrow, display.progress = FALSE)

marrow <- RunPCA(object = marrow, pc.genes = marrow@var.genes, pcs.print = 1:4,
                 genes.print = 10)

##PCHeatmap(object = marrow, pc.use = 4, do.balanced = TRUE, label.columns = FALSE, remove.key = FALSE)
PCHeatmap(object = marrow, pc.use = 4, do.balanced = TRUE, label.columns = FALSE, remove.key = TRUE)

marrow <- CellCycleScoring(object = marrow, s.genes = s.genes, g2m.genes = g2m.genes,
    set.ident = TRUE)

## view cell cycle scores and phase assignments
head(x = marrow@meta.data)
table(marrow@meta.data$Phase)

## Visualize the distribution of cell cycle markers across
##RidgePlot(object = marrow, features.plot = c("PCNA", "TOP2A", "MCM6", "MKI67"),  nCol = 2)
marrow <- RunPCA(object = marrow, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
PCAPlot(object = marrow)

## global PCA
marrow <- RunPCA(object = marrow, pc.genes = NULL, do.print = FALSE)
PCAPlot(object = marrow)

y <- log2(1 + as.matrix(exp.mat)["P2RX7",])
phase <- marrow@meta.data$Phase
boxplot( y ~ phase)
table(y>4)

plot( y, marrow@meta.data$S.Score, ylab="P2RX7")
plot( y, marrow@meta.data$G2M.Score, ylab="P2RX7")

plot( marrow@meta.data$S.Score, marrow@meta.data$G2M.Score,
     xlab="S.Score", ylab="G2M.Score", col=marrow@meta.data$Phase)
abline(h=0, v=0, lty=2)



##----------------------------------------------------------------------
## For each gene, Seurat models the relationship between gene
## expression and the S and G2M cell cycle scores. The scaled
## residuals of this model represent a ‘corrected’ expression matrix,
## that can be used downstream for dimensional reduction.
marrow <- ScaleData(object = marrow, vars.to.regress = c("S.Score", "G2M.Score"),
    display.progress = FALSE)

##marrow <- RunPCA(object = marrow, pc.genes = marrow@var.genes, genes.print = 10)
marrow <- RunPCA(object = marrow, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
PCAPlot(object = marrow)


## Alternate Workflow. As an alternative, we suggest regressing out
## the difference between the G2M and S phase scores.
marrow@meta.data$CC.Difference <- marrow@meta.data$S.Score - marrow@meta.data$G2M.Score
marrow <- ScaleData(object = marrow, vars.to.regress = "CC.Difference", display.progress = FALSE)
##marrow <- RunPCA(object = marrow, pc.genes = marrow@var.genes, genes.print = 10)
marrow <- RunPCA(object = marrow, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
PCAPlot(object = marrow)


##=====================================================================================
##=====================================================================================
##=====================================================================================
