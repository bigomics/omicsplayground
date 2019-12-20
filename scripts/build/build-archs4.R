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
cwd
setwd(archs4dir)
system("wget https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5")
system("wget https://s3.amazonaws.com/mssm-seq-matrix/human_hiseq_transcript_v2.h5")
system("wget https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix.h5")
system("wget https://s3.amazonaws.com/mssm-seq-matrix/mouse_hiseq_transcript_v2.h5")
system("wget ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz")
system("wget ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz")

create_augmented_feature_info(datadir=archs4dir)
archs4_local_data_dir_validate(datadir=archs4dir)
setwd(cwd)
cwd

## EXAMPLES
##

library(archs4)
a4 <- Archs4Repository(archs4dir)
ids <- c('GSE89189', 'GSE29943', "GSM1095128", "GSM1095129", "GSM1095130")
sample.info <- sample_info(a4, ids)
head(sample.info)

library(dplyr)
samples <- sample_info(a4, "GSE52564", columns = "Sample_title")
select(samples, sample_id, Sample_title)
sample_covariates(a4) %>% head()

yg <- my.as.DGEList(a4, "GSE89189", feature_type="gene")
yg <- my.as.DGEList(a4, "GSE53986", feature_type="gene")

dim(yg)
names(yg)

dim(yg$counts)
dim(yg$samples)
dim(yg$genes)
table(yg$samples[,"Sample_source_name_ch1"])



x = a4
id = "GSE89189"
id = "GSE53986"

features = NULL; sample_columns = c("Sample_title", "Sample_source_name_ch1"); feature_type = "gene"; row_id = "symbol";check_missing_samples = TRUE




my.as.DGEList <- function (x, id, features = NULL, sample_columns = c("Sample_title", 
    "Sample_source_name_ch1"), feature_type = c("gene", "transcript"), 
    row_id = c("ensembl", "symbol"), check_missing_samples = TRUE) 
{
    ##assert_class(x, "Archs4Repository")
    ##feature_type <- match.arg(feature_type)
    feature_type <- feature_type[1]
    ##row_id <- match.arg(row_id)
    row_id <- row_id[1]
    if (!is.null(features)) {
        if (is.character(features)) {
            features <- feature_lookup(x, features, type = type)
            if (any(is.na(features$a4name))) {
                warning("Removing 'not found' features from query")
                features <- filter(features, !is.na(a4name))
            }
        }
        assert_data_frame(features, min.rows = 1L)
        assert_subset(c("ensembl_id", "a4name"), colnames(features))
        features <- distinct(features, ensembl_id, .keep_all = TRUE)
    }
    si <- sample_info(x, id, columns = sample_columns, check_missing_samples = check_missing_samples)
    si <- as.data.frame(si, stringsAsFactors = FALSE)
    si <- distinct(si, sample_id, .keep_all = TRUE)
    rownames(si) <- si$sample_id
    not.found <- filter(si, is.na(organism))
    if (nrow(not.found)) {
        stop("The following samples could not be found: ", paste(sprintf("%s::%s", 
            not.found$series_id, not.found$sample_id), collapse = "; "))
    }
    org <- unique(si$organism)
    if (length(org) != 1L) {
        stop("You are querying across species")
    }
    finfo <- feature_info(x, feature_type, org, augmented = TRUE)
    finfo <- as.data.frame(finfo, stringsAsFactors = FALSE)
    if (feature_type == "gene") {
        if (row_id == "ensembl") {
            finfo <- filter(finfo, !is.na(ensembl_id))
            rownames(finfo) <- finfo[["ensembl_id"]]
        }
        else {
            rownames(finfo) <- finfo[["a4name"]]
        }
    }
    else {
        dup.ensid <- duplicated(finfo[["ensembl_id"]])
        if (any(dup.ensid)) {
            warning("Duplicated ensembl identifiers when version is removed, ", 
                "rownames maintain their versioned id")
            rownames(finfo) <- finfo[["ensembl_id_full"]]
        }
        else {
            rownames(finfo) <- finfo[["ensembl_id"]]
        }
    }
    if (!is.null(features)) {
        finfo <- subset(finfo, a4name %in% features[["a4name"]])
    }
    h5col <- paste0("sample_h5idx_", feature_type)
    isna <- is.na(si[[h5col]])
    if (any(isna)) {
        msg <- paste0("The following samples do not have ", feature_type, 
            "-level quantitation and will not be included in the expression ", 
            "container:\n  ", paste(sprintf("%s::%s", si$series_id[isna], 
                si$sample_id[isna]), collapse = "; "))
        warning(msg, immediate. = TRUE)
        si <- si[!isna, , drop = FALSE]
    }
    if (nrow(si) == 0L) {
        stop("No samples left to assemble expression data")
    }
    

    libinfo <- sample_table(x) %>% select(series_id, sample_id, 
                                          ## lib.size = libsize, norm.factors = normfactor,
                                          a4libsize) %>% 
        distinct(sample_id, .keep_all = TRUE)
    libinfo <- select(si, series_id, sample_id) %>% left_join(libinfo, 
        by = c("series_id", "sample_id"))
    na.overflow <- is.na(libinfo$lib.size) | is.na(libinfo$norm.factors)
    if (any(na.overflow)) {
        warning("Removing ", sum(na.overflow), " samples due to libsize NA overflow issues")
        libinfo <- filter(libinfo, !na.overflow)
        si <- subset(si, sample_id %in% libinfo$sample_id)
    }
    rownames(si) <- si[["sample_id"]]
    counts <- local({
        key <- paste(org, feature_type, sep = "_")
        h5.fn <- file_path(x, key)
        index <- list(finfo$h5idx, si[[h5col]])
        cnts <- rhdf5::h5read(h5.fn, "data/expression", index = index)
        colnames(cnts) <- rownames(si)
        rownames(cnts) <- rownames(finfo)
        cnts
    })
    out <- suppressWarnings(edgeR::DGEList(counts, genes = finfo, 
        samples = si))
    xref <- match(colnames(out), libinfo$sample_id)
    if (any(is.na(xref))) {
        stop("Problem matching sample_id to libinfo data.frame")
    }
    if (!all(colnames(out) == libinfo$sample_id[xref])) {
        stop("Mismatch in outgoing DGEList to libinfo data.frame")
    }
    out$samples$lib.size <- libinfo$lib.size[xref]
    out$samples$norm.factors <- libinfo$norm.factors[xref]
    out
}
