#!/usr/bin/env Rscript
##
## Path B: build a params.RData from the fixture CSVs + a scenario the way a
## plain script would, i.e. raw counts + `preprocess` settings, no precomputed X.
##
##   R_LIBS_USER=libs/pb-edgy Rscript run_script.R <scenario-id>
##
## Writes runs/<scenario>/B/params.RData. Executing it is run_pgx.R's job.
## See docs/script-pipeline.md for the params contract and the counts trap.
##

suppressMessages(library(playbase))

HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 1) stop("usage: run_script.R <scenario-id>")
scenario_id <- argv[1]

cfg <- jsonlite::fromJSON(file.path(HERE, "scenarios.json"), simplifyDataFrame = FALSE)
raw <- Filter(function(s) identical(s$id, scenario_id), cfg$scenarios)
if (length(raw) != 1) stop("unknown scenario: ", scenario_id)
raw <- raw[[1]]

## defaults <- scenario, one level deep (same merge drive_app.mjs does)
sc <- utils::modifyList(cfg$defaults, raw)
sc$qc <- utils::modifyList(cfg$defaults$qc, if (is.null(raw$qc)) list() else raw$qc)
sc$compute <- utils::modifyList(cfg$defaults$compute, if (is.null(raw$compute)) list() else raw$compute)
fx <- cfg$fixtures[[sc$fixture]]
if (is.null(fx)) stop("unknown fixture: ", sc$fixture)

qc <- sc$qc
cmp <- sc$compute
out_dir <- file.path(HERE, "runs", scenario_id, "B")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("[run_script] scenario=", scenario_id, " fixture=", sc$fixture, " (", fx$datatype, ")")

## ------------------------------------------------------------------ inputs
## Everything up to the app's checked_samples_counts()$COUNTS boundary.
counts <- playbase::read_counts(file.path(fx$dir, "counts.csv"))
samples <- playbase::read_samples(file.path(fx$dir, "samples.csv"))
contrasts <- playbase::read_contrasts(file.path(fx$dir, "contrasts.csv"))

## The app does NOT pass the raw contrasts through. The upload flow converts them
## to a LABEL matrix ("notact" / "act12h") before they reach params, whereas an
## old-style contrasts.csv on disk holds +1/-1. Verified against a captured
## params.RData: A held labels, an unconverted B held "-1". A script that skips
## this feeds a different contrast encoding into createPGX than the app does.
contrasts <- playbase::contrasts.convertToLabelMatrix(contrasts, samples)
annot_table <- playbase::read_annot(file.path(fx$dir, "counts.csv"))

## checkINPUT + the e29 log-scale back-conversion, exactly as
## opg-exampledata/check.R:128-134 does it.
res <- playbase::pgx.checkINPUT(counts, "COUNTS")
if ("e29" %in% names(res$checks)) {
  message("[run_script] e29: data is log-scale, converting back to linear")
  res$df <- 2**res$df
  if (min(res$df, na.rm = TRUE) >= 1) res$df <- res$df - 1
  counts <- res$df
}
message("[run_script] counts: ", nrow(counts), " x ", ncol(counts),
  " (", sum(is.na(counts)), " NA, ", sum(counts == 0, na.rm = TRUE), " zeros)")

## Probe type: detected AND reduced exactly the way the app does
## (detection at upload_server.R:41-49, reduction at 1160-1187). The app keys the
## result by organism and collapses with "+", so mirror both steps -- storing the
## raw list here would show up as a spurious gate-1 mismatch.
detected <- tryCatch(
  playbase::check_species_probetype(
    probes = rownames(counts),
    datatype = fx$datatype,
    test_species = unique(c(fx$organism, c("Human", "Mouse", "Rat"))),
    annot.cols = NULL
  ),
  error = function(e) { message("[run_script] probetype detection failed: ", e$message); NULL }
)
org <- fx$organism
task_failed <- length(detected) == 0 || !(org %in% names(detected)) ||
  is.null(detected[[org]]) || all(is.na(detected[[org]]))
probe_type <- if (task_failed) "error" else paste(detected[[org]], collapse = "+")
message("[run_script] probe_type = ", probe_type)

## ------------------------------------------------------- preprocess options
## Transcribed from the edgy app's `preprocess` reactive
## (upload_module_normalization.R:1504-1520). filter_threshold stays CHARACTER.
preprocess <- list(
  datatype = fx$datatype,
  is_npx = isTRUE(fx$is_npx),
  zero_as_na = isTRUE(qc$zero_as_na),
  normalize = isTRUE(qc$normalize),
  norm_method = qc$normalization_method,
  ref_gene = if (identical(qc$normalization_method, "reference")) qc$ref_gene else NULL,
  filter_missing = isTRUE(qc$filtermissing),
  filter_threshold = qc$filterthreshold,
  impute = isTRUE(qc$impute),
  impute_method = qc$impute_method,
  remove_outliers = isTRUE(qc$remove_outliers),
  outlier_threshold = qc$outlier_threshold,
  meth_type = fx$meth_type
)

## The app's exported norm_method is the ENCODED one: with normalization switched
## off it reports "skip_normalization", not the (still-selected) method
## (upload_module_normalization.R:1339-1343). Note this is only the reported/stored
## value -- `preprocess$norm_method` above keeps the real method, because
## pgx.preprocess decides via the `normalize` flag and still needs the method name
## for the CPM|TMM log2-prior branch.
norm_method_reported <- if (isTRUE(qc$normalize)) qc$normalization_method else "skip_normalization"

## ------------------------------------------------------------ batch correct
## Mirrors the app's bc_method reactive (upload_module_normalization.R:1351-1360).
if (isTRUE(qc$batchcorrect)) {
  batch.correct.method <- qc$bec_method
  batch.pars <- if (is.null(qc$bec_param)) "<autodetect>" else unlist(qc$bec_param)
} else {
  batch.correct.method <- "no_batch_correct"
  batch.pars <- "<autodetect>"
}

## ------------------------------------------------------------ filter flags
## upload_module_computepgx.R:1085-1100
flt <- unlist(cmp$filter_methods)
has <- function(x) x %in% flt
append.symbol <- has("append.symbol")

## The app always folds meta.go + infer in, even when do_extra is FALSE.
extra.methods <- if (isTRUE(cmp$do_extra)) unlist(cmp$extra_methods) else character(0)
extra.methods <- unique(c("meta.go", "infer", extra.methods))

## Same reader global.R:198 uses, so the limits match the app's exactly.
opts <- playbase::pgx.readOptions(file = file.path(HERE, "..", "..", "etc", "OPTIONS"))
max.genes <- as.integer(opts$MAX_GENES)
max.genesets <- as.integer(opts$MAX_GENESETS)

## The app never sends NULL for these three; mirroring it keeps gate 1 focused on
## real divergence instead of harness noise. All three confirmed against a captured
## params.RData from a live browser session.

## upload_module_computepgx.R initialises this even when no GMT is uploaded.
custom_geneset <- list(gmt = NULL, info = NULL)

## sc_compute_settings.PARS is built unconditionally (computepgx.R:1139-1152), so a
## bulk run still carries all 8 flags as FALSE rather than a NULL.
sc_compute_settings <- list(
  nfeature_threshold = FALSE, mt_threshold = FALSE, hb_threshold = FALSE,
  compute_supercells = FALSE, regress_mt = FALSE, regress_hb = FALSE,
  regress_ribo = FALSE, regress_ccs = FALSE
)

## computepgx.R:1152-1161: proteomics resolves to Olink NPX / Nulisa NPQ / MS from
## the is.olink()/is.nulisa() reactives; every other datatype stays NULL.
datatype_subtype <- if (identical(fx$datatype, "proteomics")) {
  if (isTRUE(fx$is_npx)) fx$datatype_subtype else "MS"
} else {
  fx$datatype_subtype
}

## ----------------------------------------------------------------- params
params <- list(
  organism = fx$organism,
  samples = samples,
  counts = counts,      ## RAW -- pgx.preprocess re-derives the reduction
  countsX = NULL,       ## MUST be NULL or createPGX ignores `preprocess`
  preprocess = preprocess,
  azimuth_ref = NULL,
  contrasts = as.matrix(contrasts),
  probe_type = probe_type,
  annot_table = annot_table,
  custom.geneset = custom_geneset,
  custom_fc = NULL,
  norm_method = norm_method_reported,
  settings = list(
    ## The app stores a LIST here, not a bare method string
    ## (upload_module_normalization.R:1331-1337).
    imputation_method = list(
      zero_as_na = isTRUE(qc$zero_as_na),
      imputation = if (isTRUE(qc$impute)) qc$impute_method else "no_imputation"
    ),
    bc_method = if (isTRUE(qc$batchcorrect)) list(method = qc$bec_method, param = batch.pars) else "no_batch_correct",
    remove_outliers = if (isTRUE(qc$remove_outliers)) qc$outlier_threshold else "no_outlier_removal",
    norm_method = norm_method_reported,
    custom_fc = NULL
  ),
  sc_compute_settings = sc_compute_settings,
  ## edgy-only. Drives optional AI report generation inside computePGX, not the
  ## computed data. The app builds it from live session state (llm_model +
  ## a credentials closure), which a script has no equivalent of -- NULL means
  ## "no AI reports", which is what a headless run wants.
  ai_features = NULL,
  prune.samples = TRUE,
  filter.genes = has("remove.notexpressed"),
  exclude.genes = NULL,
  remove.xy.probes = has("remove.xy.probes"),
  meth_type = fx$meth_type,
  only.known = has("remove.unknown"),
  average.duplicated = has("average.duplicated"),
  only.proteincoding = FALSE, ## DEPRECATED in the app too
  only.hugo = append.symbol,
  convert.hugo = append.symbol,
  batch.correct.method = batch.correct.method,
  batch.pars = batch.pars,
  covariates = NULL,
  dma = NULL,
  do.cluster = TRUE,
  cluster.contrasts = FALSE,
  max.genes = max.genes,
  max.genesets = max.genesets,
  gx.methods = unlist(cmp$gene_methods),
  dotimeseries = isTRUE(cmp$dotimeseries),
  gset.methods = unlist(cmp$gset_methods),
  extra.methods = extra.methods,
  use.design = FALSE,
  libx.dir = Sys.getenv("LIBX_DIR", "/home/massagno/bigomics/GitHub/libx"),
  name = paste0("cmp_", scenario_id),
  datatype = fx$datatype,
  datatype_subtype = datatype_subtype,
  description = sc$description,
  metadata = NULL,
  creator = "compute-comparer",
  date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  pgx.save.folder = out_dir
)

saveRDS(params, file = file.path(out_dir, "params.RData"))
message("[run_script] wrote ", file.path(out_dir, "params.RData"))
