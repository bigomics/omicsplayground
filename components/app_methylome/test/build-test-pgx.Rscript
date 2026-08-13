## Build a small methylomics pgx from the GSE43976 stratified subset.
## Phenotype is derived from the data (Horvath DNAm age split at the median),
## because the GEO matrix ships anonymised sample names and no metadata.
suppressMessages({library(playbase); library(wateRmelon)})

## Local env has playdata 0.1.0, which predates playdata::PHOSPHOSITE that
## playbase 4.1.2 dereferences during probe-type detection. Nothing to do with
## methylation - stub the detector so the build can proceed.
try(assignInNamespace(
  "annotate_phospho_residue",
  function(features, detect.only = FALSE) if (detect.only) FALSE else NULL,
  ns = "playbase"
), silent = TRUE)

setwd("/tmp/claude-1000/-home-massagno-bigomics-GitHub-omicsplayground-methyl-analysis/4d9e655f-36af-48da-8177-14e37d553518/scratchpad/meth")
d <- readRDS("computed.rds")
X <- d$X
A <- d$age

age <- A[, "horvath.horvath.age"]
grp <- ifelse(age >= median(age), "older", "younger")

samples <- data.frame(
  age_group    = grp,
  dnam_age     = round(age, 1),
  hannum_age   = round(A[, "hannum.hannum.age"], 1),
  phenoage     = round(A[, "phenoage.phenoage.age"], 1),
  skinblood    = round(A[, "skinblood.skinblood.age"], 1),
  row.names    = colnames(X),
  stringsAsFactors = FALSE
)

contrasts <- data.frame(
  older_vs_younger = ifelse(samples$age_group == "older", "older", "younger"),
  row.names = rownames(samples),
  stringsAsFactors = FALSE
)

message("[build] X: ", nrow(X), " x ", ncol(X),
        " | groups: ", paste(names(table(grp)), table(grp), collapse = " "))

pgx <- playbase::pgx.createPGX(
  counts       = X,
  samples      = samples,
  contrasts    = contrasts,
  organism     = "Human",
  datatype     = "methylomics",
  meth_type    = "450K array",
  dma          = "Differentially methylated positions",
  name         = "GSE43976-methylation-mini",
  creator      = "methylome-profiler demo",
  description  = "GSE43976 whole-blood 450K subset (22,080 probes x 95 samples). Stratified to retain all Horvath/Hannum clock probes and imprinted DMRs. Phenotype derived from DNAm age.",
  X            = X,
  is.logx      = FALSE,
  norm_method  = "none",
  filter.genes = FALSE,
  only.known   = FALSE,
  only.proteincoding = FALSE,
  max.genesets = 200
)

message("[build] createPGX done. slots: ", paste(names(pgx), collapse = ", "))

pgx <- playbase::pgx.computePGX(
  pgx,
  max.genes     = 25000,
  gx.methods    = c("trend.limma"),
  gset.methods  = c("fisher"),
  extra.methods = c(),
  do.cluster    = TRUE
)

## Derived per-sample metrics live on the sample sheet, following the .cell_cycle precedent
pgx$samples$.dnam_age  <- round(age, 1)
pgx$samples$.age_accel <- round(as.numeric(residuals(lm(age ~ seq_along(age)))), 2)

saveRDS(pgx, "GSE43976-methyl-mini.pgx")
message("[build] saved. X dim: ", paste(dim(pgx$X), collapse = " x "),
        " | genes cols: ", paste(head(colnames(pgx$genes), 12), collapse = ","))
