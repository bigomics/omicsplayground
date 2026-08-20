##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## The differential-methylation model, fitted in-app so covariates can be
## chosen at analysis time.
##
## Why this exists: the result stored on the pgx at upload is an unadjusted
## two-group comparison. In bulk tissue that is usually a test of cell
## composition rather than of the phenotype - most published blood CpG-age
## associations lose significance once composition is adjusted for (Jaffe &
## Irizarry 2014, Genome Biol 15:R31). Anything that shifts the neutrophil to
## lymphocyte ratio (infection, obesity, smoking, age) produces a large,
## highly significant and meaningless hit list. So the board refits limma
## itself with whatever covariates the user selects.

## --------------------------------------------------------------- masking ---

## Probes whose signal is genotype rather than methylation: a common SNP in
## the probe body, at the target CpG, or at the single-base extension site.
## Roughly 10% of the 450K array. Read straight from the manifest we already
## load, so this costs nothing.
mp_snp_probes <- function(pgx, maf = 0.05) {
  ann <- mp_manifest(pgx)
  if (is.null(ann)) return(character(0))
  cols <- intersect(c("Probe_maf", "CpG_maf", "SBE_maf"), colnames(ann))
  if (!length(cols)) return(character(0))
  hit <- Reduce(`|`, lapply(cols, function(k) {
    v <- suppressWarnings(as.numeric(ann[[k]]))
    !is.na(v) & v > maf
  }))
  out <- rownames(ann)[hit]
  ## ~10% of the 450K carries a common SNP; zero means the manifest did not
  ## load and the user would otherwise get an unmasked run believing it masked.
  if (!length(out)) warning("[methylome] SNP mask matched no probes - check the manifest")
  out
}

## Cross-reactive / multi-mapping probes, from the bundled published lists.
##
## Anchored on OPG rather than on the working directory: the app runs with its
## cwd in components/app/R, and a candidate list of relative paths alone once
## resolved to nothing there, leaving this mask a silent no-op - 53,498 probes
## the checkbox claimed to exclude were being tested, while the summary line
## still reported the SNP mask's count and so looked like both had applied.
## pgx.system.file() is not usable here: it only ever resolves
## components/board.<name>, and this is an app, not a board.
##
## The relative fallbacks are for the headless test harness, which runs from
## the test/ directory without the app's globals.
MP_XREACTIVE_CACHE <- new.env(parent = emptyenv())

mp_xreactive_probes <- function() {
  if (!is.null(MP_XREACTIVE_CACHE$p)) return(MP_XREACTIVE_CACHE$p)
  rel <- file.path("masking", "cross_reactive_probes.csv")
  cand <- character(0)
  if (exists("OPG")) cand <- file.path(OPG, "components", "app_methylome", "inst", rel)
  cand <- c(cand, file.path("../inst", rel),
            file.path("components/app_methylome/inst", rel))
  f <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(f)) {
    ## Refuse rather than warn. A warning goes to a log nobody reads while the
    ## checkbox stays ticked and excludes nothing - which is exactly how 53,498
    ## cross-reactive probes were silently tested once already. The file is
    ## 2.8 MB of published blacklist and .gitignore's blanket *.csv rule kept
    ## it out of the repo for a while; if it is missing again, say so loudly.
    MP_XREACTIVE_CACHE$p <- character(0)
    shiny::validate(shiny::need(FALSE, paste(
      "The cross-reactive probe list is missing from this installation, so that",
      "mask cannot be applied. Untick 'Cross-reactive' to run without it -",
      "leaving it ticked would report a mask that did not happen.")))
  }
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  MP_XREACTIVE_CACHE$p <- unique(as.character(d$probe))
  MP_XREACTIVE_CACHE$p
}

## getAnnotation() dispatches through the annotation object, which must be
## ATTACHED, not merely loaded - with only requireNamespace() it fails with
## "no item called package:... on the search list" and the mask silently
## becomes a no-op. Same trap as agep() and data("age_coefficients").
## Array type, inferred from the probe IDs rather than from metadata.
##
## playbase takes meth_type as an argument to pgx.createPGX and uses it when
## annotating, but does not store it on the pgx - so pgx$meth_type is NULL and
## anything keyed on it silently falls through to 450K. Inferring from the
## probes themselves is both more robust and works on subsets, where a probe
## count would tell you nothing: EPIC keeps roughly 453k of the 450K probes
## and adds ~413k of its own, so any real EPIC dataset carries a large
## fraction of probes that do not exist on the 450K at all.
MP_ARRAY_CACHE <- new.env(parent = emptyenv())

MP_ANNO_PKG <- c(`450K` = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                 EPIC = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

## Attach-and-read a manifest. getAnnotation() dispatches through the
## annotation object, which must be ATTACHED, not merely loaded - with only
## requireNamespace() it fails with "no item called package:... on the search
## list" and callers silently get nothing.
mp_read_manifest <- function(pkg) {
  key <- paste0("anno_", pkg)
  if (!is.null(MP_ARRAY_CACHE[[key]])) return(MP_ARRAY_CACHE[[key]])
  if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)
  if (!paste0("package:", pkg) %in% search()) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  ann <- tryCatch(as.data.frame(minfi::getAnnotation(getExportedValue(pkg, pkg))),
                  error = function(e) {
                    warning("[methylome] could not read ", pkg, ": ",
                            conditionMessage(e)); NULL
                  })
  MP_ARRAY_CACHE[[key]] <- ann
  ann
}

## EPICv2 probe ids carry a replicate suffix - cg00000029_TC21 - so they match
## no entry in either manifest we ship. Left undetected the EPIC-only fraction
## comes out at zero and the array is reported as 450K, which is the worst
## outcome: annotation, SNP masking and the gometh background all degrade with
## nothing said. Detect it by the suffix and refuse.
MP_EPICV2_RE <- "^cg[0-9]+_[A-Z]{2}[0-9]{2}$"

#' Infer the array from the probes present. Returns "450K", "EPIC" or "EPICv2".
mp_array_type <- function(pgx, min_frac = 0.01) {
  ids <- rownames(pgx$X)
  key <- paste0("type_", length(ids), "_", digest_ids(ids))
  if (!is.null(MP_ARRAY_CACHE[[key]])) return(MP_ARRAY_CACHE[[key]])

  if (mean(grepl(MP_EPICV2_RE, ids)) > min_frac) {
    MP_ARRAY_CACHE[[key]] <- "EPICv2"
    return("EPICv2")
  }

  epic <- mp_read_manifest(MP_ANNO_PKG[["EPIC"]])
  k450 <- mp_read_manifest(MP_ANNO_PKG[["450K"]])
  type <- "450K"
  if (is.null(epic)) {
    warning("[methylome] the EPIC manifest is not installed; assuming 450K")
  } else if (is.null(k450)) {
    type <- "EPIC"
  } else {
    epic_only <- setdiff(rownames(epic), rownames(k450))
    frac <- mean(ids %in% epic_only)
    if (frac > min_frac) type <- "EPIC"
  }
  MP_ARRAY_CACHE[[key]] <- type
  type
}

## Cheap stable key for a probe set - avoids re-inferring on every call
## without pulling in a hashing dependency.
digest_ids <- function(ids) {
  n <- length(ids)
  paste0(ids[1], "_", ids[max(1, n %/% 2)], "_", ids[n])
}

MP_EPICV2_MSG <- paste(
  "This looks like an EPIC v2 array. The EPICv2 annotation package",
  "(IlluminaHumanMethylationEPICv2anno.20a1.hg38) is not installed, so probe",
  "annotation, SNP masking and gene-set testing are not available for it -",
  "and its probe ids carry replicate suffixes that match neither the 450K nor",
  "the EPIC manifest. Everything that needs annotation is refused rather than",
  "silently computed against the wrong array."
)

## Refuse an array we cannot annotate. Called from the manifest and from the
## model, so an EWAS with masking switched off still stops here rather than
## running to completion with no annotation behind it.
mp_require_array <- function(pgx) {
  shiny::validate(shiny::need(mp_array_type(pgx) != "EPICv2", MP_EPICV2_MSG))
}

mp_manifest <- function(pgx) {
  type <- mp_array_type(pgx)
  mp_require_array(pgx)
  ann <- mp_read_manifest(MP_ANNO_PKG[[type]])
  if (is.null(ann)) warning("[methylome] no manifest available; SNP masking disabled")
  ann
}

## X and Y probes. Read from the manifest rather than from pgx$genes, whose
## chr column is a cytoband ("Xp22.2") in every playbase-built pgx and would
## need parsing; the manifest already says "chrX".
##
## Why this is a mask at all: in a mixed-sex cohort every X probe separates
## the sexes by dosage compensation alone, so any phenotype even mildly
## correlated with sex picks up a block of spurious X hits. Adjusting for sex
## does not fix it - the effect is not additive on the beta scale.
mp_sexchr_probes <- function(pgx) {
  ann <- mp_manifest(pgx)
  if (is.null(ann) || !"chr" %in% colnames(ann)) return(character(0))
  rownames(ann)[ann$chr %in% c("chrX", "chrY")]
}

## Which probes survive the selected masks.
##
## Note the default: SNP and cross-reactive are on because they are
## unconditionally wrong probes - they measure genotype or the wrong locus
## whatever the cohort. The sex-chromosome mask is conditionally right: it is
## near-universal in mixed-sex population EWAS, and plainly wrong in a
## single-sex cohort or in any study whose question is sex itself, where it
## would silently delete the signal. So it is offered, unticked.
mp_masked_probes <- function(pgx, mask = c("snp", "xreactive")) {
  out <- character(0)
  if ("snp" %in% mask) out <- c(out, mp_snp_probes(pgx))
  if ("xreactive" %in% mask) out <- c(out, mp_xreactive_probes())
  if ("sexchr" %in% mask) out <- c(out, mp_sexchr_probes(pgx))
  unique(out)
}

## ----------------------------------------------------------------- model ---

## Sample columns usable as covariates: not the contrast itself, not free text.
mp_model_vars <- function(pgx) {
  s <- pgx$samples
  if (is.null(s)) return(character(0))
  keep <- vapply(colnames(s), function(k) {
    v <- s[[k]]
    if (is.numeric(v)) return(stats::sd(v, na.rm = TRUE) > 0)
    u <- unique(as.character(v[!is.na(v) & v != ""]))
    length(u) >= 2 && length(u) <= 12
  }, logical(1))
  ## dot-prefixed columns are app-derived metrics, still valid covariates
  sort(colnames(s)[keep])
}

## Numeric sample columns usable as a continuous EWAS outcome. A large share
## of published EWAS regress on an exposure rather than compare two groups -
## age, BMI, gestational age, pack-years, dose - and limma fits those with the
## same machinery, so refusing them would be an arbitrary restriction.
## Two-valued numerics are excluded: those are a two-group contrast wearing
## numbers, and belong on the contrast path where delta-beta is meaningful.
mp_continuous_vars <- function(pgx) {
  s <- pgx$samples
  if (is.null(s)) return(character(0))
  keep <- vapply(colnames(s), function(k) {
    v <- suppressWarnings(as.numeric(as.character(s[[k]])))
    sum(!is.na(v)) >= 6 && stats::sd(v, na.rm = TRUE) > 0 &&
      length(unique(v[!is.na(v)])) > 2
  }, logical(1))
  sort(colnames(s)[keep])
}

## Categorical sample columns usable as an EWAS outcome directly.
##
## playbase builds every contrast as a two-group comparison, so a three-level
## phenotype (disease stage, treatment arm, tissue) is not merely refused by
## the model - it is unreachable from the picker. Offering the column itself is
## the only way into the F-test path.
##
## Two-level columns are included even though a contrast usually exists for
## them: it keeps one code path (the k = 2 branch is what a two-level column
## takes), and a sample column with no matching contrast would otherwise be
## untestable. At least three samples per level, or the group mean the F-test
## compares is noise.
mp_categorical_vars <- function(pgx, max_levels = 8, min_per_level = 3) {
  s <- pgx$samples
  if (is.null(s)) return(character(0))
  cont <- mp_continuous_vars(pgx)
  keep <- vapply(colnames(s), function(k) {
    if (k %in% cont) return(FALSE)
    v <- as.character(s[[k]])
    tb <- table(v[!is.na(v) & v != ""])
    length(tb) >= 2 && length(tb) <= max_levels && all(tb >= min_per_level)
  }, logical(1))
  sort(colnames(s)[keep])
}

## The outcome picker offers contrasts, continuous and categorical variables in
## one list, so which path to take is decided by what was picked. A contrast
## wins on a name clash - it is the more specific object.
mp_is_continuous <- function(pgx, sel) {
  if (is.null(sel) || !nzchar(sel)) return(FALSE)
  if (sel %in% colnames(pgx$contrasts)) return(FALSE)
  sel %in% mp_continuous_vars(pgx)
}

mp_is_categorical <- function(pgx, sel) {
  if (is.null(sel) || !nzchar(sel)) return(FALSE)
  if (sel %in% colnames(pgx$contrasts)) return(FALSE)
  sel %in% mp_categorical_vars(pgx)
}

#' Fit the differential-methylation model.
#'
#' @param contrast a contrast name, or the name of a continuous or categorical
#'   sample column.
#' @param n_sv surrogate variables to append to the design: 0 off, -1 estimate
#'   with sva::num.sv, a positive number to force that many.
#' @return list(table, formula, n, groups, strata, dropped_covars, masked, ...)
## ---------------------------------------------------------- chunked fitting ---

## Why fit in pieces at all: the peak of an EWAS is not limma, it is holding the
## M-value matrix. On an EPIC array with a few hundred samples that is gigabytes,
## and betaToM() allocates a second copy of it. Fitting a chromosome at a time
## keeps only one chromosome's M alive, and gives the progress bar something
## true to say instead of a frozen 40%.
##
## This is exact, not an approximation. lmFit is per-feature least squares, so
## the five per-feature slots below are simply the rows stacked back together;
## eBayes then runs ONCE over the combined fit, so variance moderation still
## borrows across every feature and BH still adjusts over the whole array.
## Moderating per chromosome would have quietly changed every p-value.
mp_rbind_fits <- function(fits) {
  fits <- Filter(function(f) !is.null(f) && length(f$sigma) > 0, fits)
  out <- fits[[1]]
  out$coefficients   <- do.call(rbind, lapply(fits, function(f) f$coefficients))
  out$stdev.unscaled <- do.call(rbind, lapply(fits, function(f) f$stdev.unscaled))
  out$sigma          <- unlist(lapply(fits, function(f) f$sigma), use.names = TRUE)
  out$df.residual    <- unlist(lapply(fits, function(f) f$df.residual), use.names = FALSE)
  out$Amean          <- unlist(lapply(fits, function(f) f$Amean), use.names = TRUE)
  out
}

## Features grouped by chromosome, in karyotype order, with everything the
## annotation could not place last. Chunk labels are what the progress bar says.
mp_fit_chunks <- function(features, chr) {
  chr <- sub("(p|q|cen).*", "", sub("^chr", "", as.character(chr)))
  chr[is.na(chr) | !nzchar(chr)] <- "unplaced"
  lv <- c(as.character(1:22), "X", "Y", "unplaced")
  lv <- lv[lv %in% unique(chr)]
  lv <- c(lv, setdiff(unique(chr), lv))
  split(features, factor(chr, levels = lv))
}

## Progress that no-ops outside a Shiny progress context, so the test scripts
## can call the fit directly.
mp_tick <- function(amount, detail) {
  try(shiny::incProgress(amount, detail = detail), silent = TRUE)
  invisible(NULL)
}

## lmFit over the chunks. `M` lets the caller hand in an already-materialised
## M-value matrix (SVA needs one, so there is no point building it twice);
## otherwise each chunk is converted and dropped as it goes.
mp_lmfit_chunked <- function(B, design, chunks, M = NULL, transform = TRUE,
                             label = "") {
  n <- length(chunks)
  fits <- vector("list", n)
  for (i in seq_len(n)) {
    ii <- chunks[[i]]
    if (!length(ii)) next
    Xc <- if (!is.null(M)) M[ii, , drop = FALSE] else {
      if (transform) playbase.epigenetics::betaToM(B[ii, , drop = FALSE])
      else B[ii, , drop = FALSE]
    }
    fi <- limma::lmFit(Xc, design)
    ## lmFit drops the row names when handed a single-row matrix, and a
    ## chromosome really can hold exactly one feature - chrX did, on the demo.
    ## Unnamed rows then fail to match on the way back and eBayes dies in
    ## fitFDist with "NA covariate values not allowed". Put them back.
    rownames(fi$coefficients) <- rownames(Xc)
    rownames(fi$stdev.unscaled) <- rownames(Xc)
    names(fi$sigma) <- rownames(Xc)
    names(fi$Amean) <- rownames(Xc)
    fits[[i]] <- fi
    rm(Xc, fi)
    mp_tick(0.55 / n, sprintf("%schr %s (%s features)", label, names(chunks)[i],
                              format(length(ii), big.mark = ",")))
  }
  fit <- mp_rbind_fits(fits)
  ## Chunking reorders; hand the rows back in the order they came in.
  fit[match(rownames(B), rownames(fit$coefficients)), ]
}

## ------------------------------------------------------ feature collapsing ---

## Testing units coarser than a probe. This is a screening view, not a
## substitute for a probe-level EWAS, and the literature is specific about why:
##
##  - Averaging a WHOLE gene is not defensible. Promoter methylation silences
##    while gene-body methylation correlates positively with expression (Yang
##    2014, Cancer Cell 26:577), so the two halves cancel - Ruiz-Arenas &
##    Gonzalez 2017 (BMC Bioinformatics 18:553): "averaging the effect of the
##    entire region may provide a signal close to zero because the effects are
##    compensated". Neighbouring CpGs only co-vary within about a kilobase
##    anyway (Mansell 2019, BMC Genomics 20:366), and gene bodies are tens of
##    kb. So there is no whole-gene option here, deliberately.
##  - Splitting each gene into promoter and body is the established partition:
##    IMA (Wang 2012, Bioinformatics 28:729) indexes exactly these RefGene
##    groups, and mCSEA (Martorell-Marugan 2019, Bioinformatics 35:3257) merges
##    TSS1500/TSS200/5'UTR/1stExon into one promoter rather than testing four
##    overlapping sub-groups as if they were independent.
##  - Testing an average is statistically fine as long as the grouping comes
##    from annotation and not from the data - Smyth, Bioconductor #126448 - and
##    a median summary holds type I error at 0.0532 against a nominal 0.05
##    (Gomez 2019, Nucleic Acids Res 47:e98). Hence median, and hence collapsing before
##    anything looks at the outcome.
MP_COLLAPSE <- c(
  "Gene region (promoter / body)" = "gene_region",
  "Probe (full EWAS)" = "probe"
)

## mCSEA's partition. 3'UTR sits with the body: it is downstream sequence, not
## regulatory promoter, and folding it in keeps the feature set to two per gene.
MP_PROMOTER_GROUPS <- c("TSS1500", "TSS200", "5'UTR", "1stExon")
MP_BODY_GROUPS     <- c("Body", "3'UTR")

## Smallest number of probes a feature may be built from. A two-probe median is
## not a summary, and single-probe "features" would just be a probe-level EWAS
## with the multiple-testing burden quietly reduced.
MP_MIN_PROBES <- 3L

## Feature key per probe. NA where the probe cannot be assigned, which drops it.
##
## The manifest annotates a probe against every transcript it falls in, as
## parallel semicolon lists, so one probe can read "1stExon;3'UTR". Taking the
## first entry keeps features disjoint - each probe contributes to exactly one
## feature, so the tests are not sharing data - at the cost of preferring
## whichever isoform the manifest happened to list first.
mp_collapse_key <- function(genes, by) {
  sym <- sub(";.*", "", as.character(genes$symbol))
  sym[sym %in% c("", "NA", "-")] <- NA
  loc <- sub(";.*", "", as.character(genes$genomic_location))
  region <- ifelse(loc %in% MP_PROMOTER_GROUPS, "promoter",
            ifelse(loc %in% MP_BODY_GROUPS, "body", NA_character_))
  ifelse(is.na(sym) | is.na(region), NA_character_, paste0(sym, ":", region))
}

## The beta matrix rows for whatever the fit actually tested. At probe level
## that is a plain row subset; on a collapsed fit the features are not rows of
## pgx$X at all, so the ones being drawn are rebuilt from their member probes.
## Every panel that plots per-sample beta for a hit goes through here - indexing
## pgx$X by feature id directly is a subscript-out-of-bounds waiting to happen.
##
## Uses the same median over complete probes the fit used, so a stripchart shows
## the quantity that was tested rather than a differently-summarised cousin.
mp_feature_beta <- function(pgx, meta, features) {
  X <- mp_beta(pgx)
  by <- if (is.null(meta$collapse)) "probe" else meta$collapse
  if (identical(by, "probe")) {
    keep <- intersect(features, rownames(X))
    shiny::validate(shiny::need(length(keep) > 0,
      "None of these features is in the methylation matrix."))
    return(X[keep, , drop = FALSE])
  }
  key <- mp_collapse_key(pgx$genes[rownames(X), , drop = FALSE], by)
  keep <- !is.na(key) & key %in% features & rowSums(is.na(X)) == 0
  shiny::validate(shiny::need(any(keep),
    "None of the top features could be traced back to its probes."))
  Xs <- X[keep, , drop = FALSE]
  f <- factor(key[keep])
  out <- t(vapply(split(seq_len(nrow(Xs)), f),
                  function(i) matrixStats::colMedians(Xs[i, , drop = FALSE]),
                  numeric(ncol(Xs))))
  colnames(out) <- colnames(Xs)
  out[intersect(features, rownames(out)), , drop = FALSE]
}

## Collapse a beta matrix to features, and build the annotation the downstream
## panels read.
##
## Median, not mean: a median survives one blown probe, it is IMA's default,
## and coMethDMR measured the two as equivalent on type I error (0.0532 vs
## 0.0545 against a nominal 0.05). The cost is that Delta-beta becomes a
## difference of medians - still a methylation difference, still on the beta
## scale, but not the same quantity as a probe-level Delta-beta.
##
## Probes with any missing value across the fitted samples are dropped before
## the median rather than skipped inside it. A summary taken over a different
## probe subset per sample moves with the missingness pattern, and because
## probes within a feature sit at very different baseline beta, that shift can
## exceed the effect being looked for.
mp_collapse_beta <- function(B, genes, by, min_probes = MP_MIN_PROBES) {
  g <- genes[rownames(B), , drop = FALSE]
  key <- mp_collapse_key(g, by)
  n_in <- nrow(B)
  n_unassigned <- sum(is.na(key))
  complete <- rowSums(is.na(B)) == 0
  n_incomplete <- sum(!is.na(key) & !complete)
  ok <- !is.na(key) & complete
  shiny::validate(shiny::need(sum(ok) > 100,
    "Almost no probe carries both a gene region and a complete set of values, so this dataset cannot be collapsed. Test at probe level."))
  B <- B[ok, , drop = FALSE]; key <- key[ok]; g <- g[ok, , drop = FALSE]

  ## Features below the floor are dropped whole, not silently thinned.
  cnt <- table(key)
  big <- names(cnt)[cnt >= min_probes]
  n_small <- sum(cnt < min_probes)
  keep <- key %in% big
  shiny::validate(shiny::need(sum(keep) > 100,
    sprintf("No gene region reaches %d probes in this dataset. Test at probe level.",
            min_probes)))
  B <- B[keep, , drop = FALSE]; key <- key[keep]; g <- g[keep, , drop = FALSE]

  f <- factor(key)
  idx <- split(seq_len(nrow(B)), f)
  X <- t(vapply(idx, function(i) matrixStats::colMedians(B[i, , drop = FALSE]),
                numeric(ncol(B))))
  colnames(X) <- colnames(B)

  ## One annotation row per feature, in the same order as X.
  ng <- nlevels(f); gi <- as.integer(f)
  ## Commonest non-empty value per feature, with the count that backs it.
  ## The count is not decoration: a gene body routinely spans island, shore,
  ## shelf and open sea, so a bare majority label asserts a context most of the
  ## feature's probes may disagree with.
  pick2 <- function(v) {
    v <- as.character(v)
    keepv <- !is.na(v) & nzchar(v)
    out <- list(v = rep(NA_character_, ng), n = rep(NA_integer_, ng))
    if (!any(keepv)) return(out)
    d <- data.table::data.table(g = gi[keepv], v = v[keepv])
    d <- d[, list(n = .N), by = c("g", "v")]
    data.table::setorderv(d, c("g", "n", "v"), c(1L, -1L, 1L))
    d <- d[!duplicated(d$g)]
    out$v[d$g] <- d$v
    out$n[d$g] <- as.integer(d$n)
    out
  }
  pick <- function(v) pick2(v)$v
  ppos <- split(suppressWarnings(as.numeric(sub(";.*", "", g$pos))), f)
  ## Median for anything that needs one number on an axis - the Manhattan draws
  ## features at a position. Start and end because a feature is a span, not a
  ## point: on this array the 90th percentile spans ~10kb and the worst 1.6Mb,
  ## so a lone midpoint claims a precision the feature does not have.
  pos <- vapply(ppos, stats::median, numeric(1), na.rm = TRUE)
  pos_start <- vapply(ppos, function(z) suppressWarnings(min(z, na.rm = TRUE)), numeric(1))
  pos_end   <- vapply(ppos, function(z) suppressWarnings(max(z, na.rm = TRUE)), numeric(1))
  pos_start[!is.finite(pos_start)] <- NA_real_
  pos_end[!is.finite(pos_end)] <- NA_real_
  isl <- pick2(sub(";.*", "", g$Relation_to_Island))
  annot <- data.frame(
    symbol = sub(":(promoter|body)$", "", levels(f)),
    region = sub("^.*:", "", levels(f)),
    chr = pick(sub(";.*", "", g$chr)),
    pos = as.numeric(pos),
    pos_start = as.numeric(pos_start),
    pos_end = as.numeric(pos_end),
    genomic_location = pick(sub(";.*", "", g$genomic_location)),
    Relation_to_Island = isl$v,
    island_n = isl$n,
    n_probes = as.integer(table(f)),
    row.names = levels(f), stringsAsFactors = FALSE
  )
  list(X = X, genes = annot, n_probes_in = n_in,
       n_unassigned = n_unassigned, n_incomplete = n_incomplete,
       n_small_features = n_small, min_probes = min_probes)
}

mp_fit_ewas <- function(pgx, contrast, covars = character(0),
                        cellfracs = NULL, mask = character(0), n_sv = 0,
                        ## Probe level, i.e. a real EWAS, is the honest default
                        ## for a programmatic caller. The app's picker defaults
                        ## to gene instead, and always passes it explicitly.
                        collapse = "probe") {
  mp_require_methylomics(pgx)
  mp_require_array(pgx)
  shiny::validate(shiny::need(!is.null(contrast) && nzchar(contrast),
                              "No outcome selected."))

  beta <- mp_beta(pgx)
  continuous <- mp_is_continuous(pgx, contrast)

  if (continuous) {
    x <- suppressWarnings(as.numeric(as.character(pgx$samples[[contrast]])))
    names(x) <- rownames(pgx$samples)
    x <- x[!is.na(x)]
    ss <- intersect(names(x), colnames(beta))
    shiny::validate(shiny::need(length(ss) >= 6,
      "Too few samples carry a value for this variable."))
    dat <- data.frame(.x = as.numeric(x[ss]), row.names = ss)
  } else {
    grp <- if (mp_is_categorical(pgx, contrast)) {
      g <- as.character(pgx$samples[[contrast]])
      names(g) <- rownames(pgx$samples)
      g[!is.na(g) & g != ""]
    } else {
      mp_contrast_groups(pgx, contrast)
    }
    shiny::validate(shiny::need(!is.null(grp), "This contrast has no group labels."))
    ss <- intersect(names(grp), colnames(beta))
    grp <- droplevels(factor(grp[ss]))
    shiny::validate(shiny::need(nlevels(grp) >= 2,
      "This outcome has only one group. Pick a contrast or a variable that varies."))
    dat <- data.frame(.group = grp, row.names = ss, stringsAsFactors = FALSE)
  }

  ## Covariate frame, aligned to the outcome's samples. The outcome itself is
  ## never also a covariate.
  used <- character(0); dropped <- character(0)
  for (k in setdiff(covars, contrast)) {
    v <- pgx$samples[ss, k]
    if (is.numeric(v)) {
      if (stats::sd(v, na.rm = TRUE) > 0) { dat[[k]] <- v; used <- c(used, k) }
      else dropped <- c(dropped, paste0(k, " (constant)"))
    } else {
      f <- droplevels(factor(as.character(v)))
      if (nlevels(f) >= 2) { dat[[k]] <- f; used <- c(used, k) }
      else dropped <- c(dropped, paste0(k, " (single level)"))
    }
  }
  if (!is.null(cellfracs)) {
    cf <- cellfracs[ss, , drop = FALSE]
    ## One fraction is redundant with the rest (they sum to 1); drop the
    ## largest so the design stays full rank.
    cf <- cf[, order(-colMeans(cf, na.rm = TRUE)), drop = FALSE][, -1, drop = FALSE]
    for (k in colnames(cf)) { dat[[k]] <- as.numeric(cf[, k]); used <- c(used, k) }
  }

  ## Drop samples missing any model variable, then re-check the contrast.
  keep <- stats::complete.cases(dat)
  dat <- dat[keep, , drop = FALSE]
  shiny::validate(shiny::need(nrow(dat) >= 6, "Too few complete samples for this model."))
  if (continuous) {
    shiny::validate(shiny::need(stats::sd(dat$.x) > 0,
      "This variable is constant once incomplete covariates are dropped."))
  } else {
    dat$.group <- droplevels(dat$.group)
    shiny::validate(shiny::need(nlevels(dat$.group) >= 2,
      "Only one group has samples left once incomplete covariates are dropped."))
  }
  ## Three or more groups is a moderated F-test, not a difference of two means:
  ## no signed effect and no single t exists, which several panels care about.
  anova_fit <- !continuous && nlevels(dat$.group) >= 3

  ## A continuous outcome keeps the intercept and tests its own slope; a
  ## contrast drops it and tests the difference of the two group means.
  rhs <- c(if (continuous) ".x" else "0 + .group", used)
  design <- stats::model.matrix(
    stats::as.formula(paste("~", paste(rhs, collapse = " + "))), data = dat)
  colnames(design) <- make.names(colnames(design))

  ## Refuse a rank-deficient design rather than returning silent nonsense.
  ne <- limma::nonEstimable(design)
  shiny::validate(shiny::need(is.null(ne),
    paste0("These terms are collinear with the contrast and cannot be fitted: ",
           paste(ne, collapse = ", "), ". Remove one of them.")))
  shiny::validate(shiny::need(nrow(dat) - ncol(design) >= 2,
    "Not enough residual degrees of freedom - drop a covariate."))

  probes <- rownames(beta)
  masked <- intersect(mp_masked_probes(pgx, mask), probes)
  probes <- setdiff(probes, masked)
  shiny::validate(shiny::need(length(probes) > 100, "Masking removed nearly every probe."))

  ## Fit on M-values, report on beta. M is the right scale for testing
  ## (Du 2010) but is not interpretable as a change in methylation.
  B <- beta[probes, rownames(dat), drop = FALSE]
  ## Collapse after masking, not before: the SNP and cross-reactive lists are
  ## probe ids, so a masked probe has to be gone before its beta is averaged
  ## into a feature. Collapsing a working copy here - pgx is never touched, so
  ## switching resolution is a refit, not a re-upload.
  collapsed <- NULL
  if (!identical(collapse, "probe")) {
    collapsed <- mp_collapse_beta(B, pgx$genes, collapse)
    B <- collapsed$X
    shiny::validate(shiny::need(nrow(B) > 50,
      "Collapsing left too few features to test."))
  }
  ## Chromosome is the chunk key. On a collapsed fit it comes off the collapsed
  ## annotation, since those features are not rows of pgx$genes.
  chr_src <- if (is.null(collapsed)) pgx$genes[rownames(B), , drop = FALSE]
             else collapsed$genes[rownames(B), , drop = FALSE]
  cc <- grep("^chr$|chrom", tolower(colnames(chr_src)))[1]
  chunks <- mp_fit_chunks(rownames(B),
                          if (!is.na(cc)) chr_src[[cc]] else rep(NA, nrow(B)))

  ## SVA estimates genome-wide latent structure, so it needs the whole matrix -
  ## there is no per-chromosome version of that question. When it is on we pay
  ## for one full M and then reuse it for the fit rather than rebuilding chunks.
  M <- if (n_sv != 0) playbase.epigenetics::betaToM(B) else NULL

  ## ------------------------------------------------------- latent factors --
  ## Surrogate variables: the only confounding control left when no
  ## deconvolution reference exists for the tissue. The outcome is protected by
  ## sitting in mod and not in mod0, so sva estimates structure orthogonal to
  ## the question rather than removing the answer.
  sv_note <- NULL
  if (n_sv != 0) {
    shiny::validate(shiny::need(requireNamespace("sva", quietly = TRUE),
      "The sva package is not installed, so latent-factor adjustment is unavailable."))
    outcome_cols <- if (continuous) ".x" else grep("^\\.group", colnames(design), value = TRUE)
    mod0 <- design[, setdiff(colnames(design), outcome_cols), drop = FALSE]
    ## The contrast design is fitted without an intercept (~ 0 + .group), so
    ## dropping the group columns can leave a null model with no constant term.
    if (!any(apply(mod0, 2, function(z) length(unique(z)) == 1))) {
      mod0 <- cbind(Intercept = 1, mod0)
    }
    sv <- tryCatch({
      k <- if (n_sv < 0) sva::num.sv(M, design, method = "be") else as.integer(n_sv)
      ## Cap: ten SVs is already more latent structure than a typical cohort
      ## supports, and every SV costs a residual degree of freedom.
      k <- min(k, 10L, nrow(dat) - ncol(design) - 2L)
      if (k < 1) NULL else {
        ## sva narrates its IRW iterations on stdout; the progress bar already
        ## says what is happening.
        invisible(utils::capture.output(
          obj <- sva::sva(M, mod = design, mod0 = mod0, n.sv = k)))
        obj$sv
      }
    }, error = function(e) e)
    shiny::validate(shiny::need(!inherits(sv, "error"),
      paste0("Surrogate-variable estimation failed on this design: ",
             if (inherits(sv, "error")) conditionMessage(sv) else "",
             ". Untick the latent-factor adjustment, or drop a covariate.")))
    if (is.null(sv)) {
      ## Degrade to the unadjusted model, but say so - a silent fall-back is
      ## exactly the fail-open this board exists to prevent.
      sv_note <- "no latent factors detected"
    } else {
      colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
      design <- cbind(design, sv)
      used <- c(used, colnames(sv))
      shiny::validate(shiny::need(nrow(dat) - ncol(design) >= 2,
        "Not enough residual degrees of freedom once the surrogate variables are added - drop a covariate."))
    }
  }

  if (continuous) {
    ## One eBayes over the combined fit, never per chunk: the moderation has to
    ## borrow across every feature or the p-values are not the ones a
    ## whole-array fit would give.
    fit2 <- limma::eBayes(mp_lmfit_chunked(B, design, chunks, M), trend = TRUE)
    tt <- limma::topTable(fit2, coef = ".x", number = Inf, sort.by = "none")
    ## Effect size on the beta scale is the slope of beta on the outcome, so
    ## the same design is refitted on B. Reported per unit of the outcome -
    ## a delta-beta of 0.004 for age means 0.4% methylation per year.
    dbeta <- mp_lmfit_chunked(B, design, chunks, transform = FALSE,
                              label = "effect size, ")$coefficients[, ".x"]
    head <- paste0("~ ", contrast, " (continuous)")
    ## Tertiles, so the per-CpG and per-region panels have one code path
    ## whether the outcome is a contrast or a continuous exposure.
    br <- unique(c(-Inf, stats::quantile(dat$.x, c(1 / 3, 2 / 3), na.rm = TRUE), Inf))
    strata <- factor(cut(dat$.x, br,
                         labels = if (length(br) == 4) c("low", "mid", "high") else NULL))
    groups <- NULL
    desc <- sprintf("%s from %.4g to %.4g", contrast, min(dat$.x), max(dat$.x))
  } else {
    lv <- levels(dat$.group)
    gcol <- make.names(paste0(".group", lv))
    fit <- mp_lmfit_chunked(B, design, chunks, M)
    ## Every level against the first. With one contrast limma returns the
    ## moderated t; with k-1 it returns the moderated F of "any difference
    ## among the groups", which is the multi-level test. Testing the k group
    ## coefficients directly would instead ask whether every group mean is
    ## zero - true of no M-value anywhere.
    cm <- limma::makeContrasts(contrasts = paste0(gcol[-1], "-", gcol[1]),
                               levels = design)
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm), trend = TRUE)
    tt <- limma::topTable(fit2, number = Inf, sort.by = "none")
    if (anova_fit) {
      ## An F-test has no direction, so the effect size is the spread of the
      ## group means on the beta scale: unsigned, still a methylation
      ## difference, and still meaningful to the |delta-beta| filter.
      mu <- lapply(lv, function(l)
        rowMeans(B[, dat$.group == l, drop = FALSE], na.rm = TRUE))
      dbeta <- do.call(pmax, mu) - do.call(pmin, mu)
      head <- sprintf("~ %s (F-test across %d groups)", contrast, length(lv))
    } else {
      dbeta <- rowMeans(B[, dat$.group == lv[2], drop = FALSE], na.rm = TRUE) -
               rowMeans(B[, dat$.group == lv[1], drop = FALSE], na.rm = TRUE)
      head <- paste0("~ ", lv[2], " vs ", lv[1])
    }
    strata <- dat$.group
    groups <- table(dat$.group)
    desc <- NULL
  }
  names(strata) <- rownames(dat)

  tab <- data.frame(
    probe = rownames(tt),
    logFC_M = if (anova_fit) NA_real_ else tt$logFC,
    dbeta = round(dbeta[rownames(tt)], 4),
    p = tt$P.Value, q = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
  ## Only on the anova branch: the F is what was tested there, and there is no
  ## single log fold change across k-1 contrasts to report in its place.
  if (anova_fit) tab$Fstat <- tt$F

  list(
    table = tab,
    formula = paste0(head,
                     if (length(used)) paste0(" + ", paste(used, collapse = " + ")) else
                       "   (unadjusted)",
                     if (!is.null(sv_note)) paste0("   [", sv_note, "]")),
    ## dmrff needs the standard error, which limma does not return directly;
    ## se = |logFC / t|, so keep the moderated t alongside. NULL on the anova
    ## branch so nothing downstream mistakes an F for a signed t.
    t = if (anova_fit) NULL else tt$t,
    anova = anova_fit,
    n_sv = if (is.null(sv_note)) sum(grepl("^SV[0-9]+$", used)) else 0L,
    n = nrow(dat), groups = groups, desc = desc,
    continuous = continuous,
    ## Which unit was tested. Every panel keyed on CpG ids reads this and
    ## refuses rather than showing an empty answer as if it were a real one.
    collapse = collapse,
    collapse_annot = if (is.null(collapsed)) NULL else collapsed$genes,
    collapse_note = if (is.null(collapsed)) NULL else sprintf(
      "%s gene regions from %s probes  -  dropped %s with no gene region, %s with missing values, and %s region(s) under %d probes",
      format(nrow(collapsed$X), big.mark = ","),
      format(collapsed$n_probes_in, big.mark = ","),
      format(collapsed$n_unassigned, big.mark = ","),
      format(collapsed$n_incomplete, big.mark = ","),
      format(collapsed$n_small_features, big.mark = ","),
      collapsed$min_probes),
    x = if (continuous) stats::setNames(dat$.x, rownames(dat)) else NULL,
    strata = strata, samples = rownames(dat),
    dropped_covars = dropped, masked = length(masked),
    covars = used
  )
}
